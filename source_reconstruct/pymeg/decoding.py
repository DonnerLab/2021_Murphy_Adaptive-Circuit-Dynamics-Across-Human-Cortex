"""
Implememt multivariate pattern classification in source space across all vertices


The idea is to couple source reconstruction with a decoding approach:

  0. To use both hemispheres at once one needs to construct special labels that contain both areas.
  1. Use lcmv.py to perform source reconstruction for one ROI.
  2. Here: Provide a custom accumulate function that performs decoding.
"""
import numpy as np
import mne
from functools import partial
from itertools import product
from scipy.stats import uniform
from sklearn import svm
from sklearn.model_selection import (
    cross_validate,
    cross_val_predict,
    RandomizedSearchCV,
)
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

#from imblearn.over_sampling import RandomOverSampler

#from conf_analysis.behavior import metadata
#from conf_analysis.meg import preprocessing

#from imblearn.pipeline import Pipeline
from sklearn.metrics.scorer import make_scorer
from sklearn.metrics import roc_auc_score, mean_squared_error
from sklearn.utils.multiclass import type_of_target
from sklearn.feature_selection import SelectFromModel, SelectKBest

from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.linear_model import LinearRegression

from sklearn.decomposition import PCA
from sklearn.utils import shuffle
import pandas as pd

from pymeg.lcmv import complex_tfr, tfr2power_estimator

njobs = 1


def get_lcmv(tfr_params, epochs, filters, njobs=njobs):
    times = epochs[0].times
    return multi_apply_lcmv(epochs, times, filters, tfr_params)


class Decoder(object):
    def __init__(self, target, classifier=None):
        """

        Args:
            classifier: (name, Sklearn classifier object)
                A tuple that contains a name at first position and a
                sklearn classifier at second position.
            target: pandas Series
                This Series needs to be indexed by trial numbers
                that are also given to lcmv.py for source
                reconstruction.

        """
        self.classifier = classifier
        self.target = target

        if classifier is None:
            self.clf_name = "SVClin"
            self.clf = None
        else:
            self.clf_name, self.clf = classifier

    def __call__(self, *args, **kwargs):
        return self.classify(*args, **kwargs)

    def classify(self, data, time, freq, trial, roi,
        average_vertices=False, use_phase=True):
        """Perform decoding on source reconstructed data.

        Decoding is carried out for each time point separately.

        Args:
            data: ndarray
                4d: ntrials x freq x vertices x time
            time: ndarray
                time points that match last dimension of data
            freq: value
                estimation value for this set of data
            trial: ndarray
                Needs to match first dim of data
            roi: str
                Name of the roi that this comes from
            average_vertices: False or int
                Average vertices per hemisphere?
                If not false must be index
                array to indicate which vertices belong
                to which hemisphere.
            use_phase: bool
                Include phase as feature?
        Returns:
            A pandas DataFrame that contains decoding accuracies
        """
        n_trials = data.shape[0]
        # This can be moved into the loop to save memory?

        target_vals = self.target    # PM: original was target.loc[trial]
        idnan = np.isnan(target_vals)
        scores = []
        
        average_vertices = np.array(average_vertices,ndmin=1)  # PM: ensures average_vertices is 1d array even if set to default of False

        for idx_t, tp in enumerate(time):
            if use_phase:
                d = data[:, :, :, idx_t].copy()
                phase = np.angle(d)
                del d
                if average_vertices.any()==True:

                    phase = np.stack(
                            (phase[:, :, average_vertices].mean(2),
                             phase[:, :, ~average_vertices].mean(2)),
                            2
                        )
            power = ((data[:, :, :, idx_t] * data[:, :, :, idx_t].conj()).real)
            if average_vertices.any()==True:
                power = np.stack(
                            (power[:, :, average_vertices].mean(2),
                             power[:, :, ~average_vertices].mean(2)),
                            2
                        )
            # Build prediction matrix:
            if use_phase:
                X = np.hstack(
                    [
                        power[:, :, :].reshape((n_trials, -1)),
                        phase[:, :, :].reshape((n_trials, -1)),
                    ]
                )
            else:
                X = power[:, :, :].reshape((n_trials, -1))
            print('Time:', tp, 'Size:', X.shape)
            if self.clf is None:
                C = 10/X.shape[1]
                print('C=', C)
                clf = Pipeline(
                    [
                       ("Scaling", StandardScaler()),
                       ("PCA", PCA(n_components=0.95, svd_solver='full')),
                       ("Upsampler", RandomOverSampler(sampling_strategy="minority")),
                       ("FeatureSelection", SelectFromModel(svm.LinearSVC(C=C, penalty="l1", dual=False, max_iter=50000))),
                       ("SVClin", svm.LinearSVC(max_iter=5000, dual=False, penalty="l2", C=1/2)),
                    ]
                )
            else:
                clf = self.clf
            s = categorize(clf, target_vals[~idnan], X[~idnan, :])
            s["latency"] = tp
            s["roi"] = roi
            scores.append(s)
        return pd.DataFrame(scores)


def categorize(clf, target, data, njobs=1):
    """
    Expects a pandas series and a pandas data frame.
    Both need to be indexed with the same index.
    """
    from sklearn.metrics import make_scorer

    corr_scorer = make_scorer(lambda x, y: np.corrcoef(x, y)[0, 1])
    metrics = {"correlation":corr_scorer}

    score = cross_validate(
        clf, data, target, cv=10, scoring=metrics, return_train_score=False, n_jobs=njobs
    )
    del score["fit_time"]
    del score["score_time"]
    score = {k: np.mean(v) for k, v in list(score.items())}
    print(score)

    return score


def multi_apply_lcmv(tfrdata, times, filters, tfr_params, max_ori_out="signed"):
    """Apply Linearly Constrained Minimum Variance (LCMV) beamformer weights.


    Args:
        tfrdata: list of ndarray
            Data to be reconstructed. Each element in this list is one set
            of epochs, which will be reconstructed by the corresponding
            filter.
            Each element should be either n_trials x n_sensors x Y x n_time
            or trials x sensors x time. Reconstruction treats epochs and
            dim Y as independent dimensions.
        est_key: value
            A key to identify this reconstruction (e.g. F for power)
        est_vals: sequence
            Values that identify different reconstructions along dimension Y
            for a single epoch, e.g. the frequency for power reconstructions.
            Needs to be length Y.
        events: list of ndarray
            Identifiers for different epochs. Needs to be of length n_trials.
        times: array
            Time of entries in last dimension of input data.
        infos: list of mne info structure
            Info structure of the epochs which are to be reconstructed.
            One element per element of tfrdata
        filters: list of filter dicts
            List of filter dicts, one for each element in tfrdata


    Returns:
        ndarray of source reconstructed epochs, events, times, est_vals
    """
    from pymeg.lcmv import _apply_lcmv

    results = []
    evs = []
    for epochs, flt in zip(tfrdata, filters):
        assert len(flt) == 1
        evs.append(epochs.events[:, 2])
        info = epochs.info

        epochs, tfrtimes, est_val, est_key = complex_tfr(
            epochs._data[:,:,:], times, **tfr_params
        )

        nfreqs = epochs.shape[2]
        info["sfreq"] = 1.0 / np.diff(tfrtimes)[0]
        eres = []
        for freq in range(nfreqs):
            filter = flt[list(flt.keys())[0]]
            mne.set_log_level("ERROR")
            data = np.stack(
                [
                    x._data
                    for x in _apply_lcmv(
                        data=epochs[:, :, freq, :],
                        filters=filter,
                        info=info,
                        tmin=tfrtimes.min(),
                        max_ori_out=max_ori_out,
                    )
                ]
            )

            eres.append(data)
        results.append(np.stack(eres, 1))
    return np.vstack(results), np.concatenate(evs), est_val, tfrtimes


def multiclass_roc(y_true, y_predict, **kwargs):
    from sklearn.preprocessing import label_binarize
    from sklearn.metrics import roc_auc_score

    return roc_auc_score(
        label_binarize(y_true, classes=[-2, -1, 1, 2]),
        label_binarize(y_predict, classes=[-2, -1, 1, 2]),
        **kwargs
    )
