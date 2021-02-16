from __future__ import division
from __future__ import print_function

import logging
import os

from itertools import product

import numpy as np
import pandas as pd

from joblib import Memory
from joblib import Parallel, delayed

from mne import compute_covariance
from mne.beamformer import make_lcmv
from mne.time_frequency.tfr import _compute_tfr


memory = Memory(cachedir=os.environ['PYMEG_CACHE_DIR'], verbose=0)

fois = np.arange(10, 150, 5)
default_tfr = {'foi': fois, 'cycles': fois * 0.1, 'time_bandwidth': 2,
               'n_jobs': 1, 'est_val': fois, 'est_key': 'F'}


def complex_tfr(x, time, est_val=None, est_key=None, sf=600., foi=None,
                cycles=None, time_bandwidth=None, n_jobs=1, decim=10):
    """Estimate power of epochs in array x."""
    if len(x.shape) == 2:
        x = x[np.newaxis, :, :]
    y = _compute_tfr(
        x, foi, sfreq=sf, method='multitaper', decim=decim, n_cycles=cycles,
        zero_mean=True, time_bandwidth=time_bandwidth, n_jobs=n_jobs,
        use_fft=True, output='complex')
    return y, time[::decim], est_val, est_key


def broadband_est(x, time, est_val=[-1], est_key='BB', **kwargs):
    """Estimate broadband from source reconstructed time course"""
    return x, time, est_val, est_key


def tfr2power_estimator(x):
    """Compute power on source reconstructed FFT results"""
    return ((x * x.conj()).real).mean(1)


def accumulate(data, time, est_key, est_val, roi, trial):
    '''Transform source reconstruction results to a DataFrane.

    Args:
        data: ndarray
            If ntrials x vertices x time in which case the
            function will average across vertices.
            If ntrials x time will be directly converted to df.
        time: ndarray
            time points that match last dimension of data
        est_key: value
            estimation key for this value
        est_val: value
            estimation value for this set of data
        roi: str
            Name of the roi that this comes from
        trial: ndarray
            Needs to match first dim of data
    Returns:
        A pandas DataFrame that contains source reconstructed data
        with hierarchical index to describe each data point.
    '''
    if data.ndim == 3:
        data = flip_and_avg_vertices(data).mean(1)

    # Now ntrials x time
    df = pd.DataFrame(data, index=trial, columns=time)
    df.columns.name = 'time'
    df.index.name = 'trial'
    df = df.stack().reset_index()
    df.loc[:, 'est_key'] = est_key
    df.loc[:, 'est_val'] = est_val
    df.set_index(['trial', 'time', 'est_key', 'est_val'], inplace=True)
    df.columns = [roi]
    return df


def flip_and_avg_vertices(data):
    """Correct random sign flips in reconstructed vertices.

    Average over vertices but correct for random flips first
    Correction is done by ensuring positive correlations
    between vertices

    Args:
        data: ndarray
            A 2D array with vertices in the first dimension and time
            in the second.

    Returns:
        A single time series constructed by averaging over vertices.
    """
    for i in range(data.shape[1]):
        if np.corrcoef(data[:, i, :].ravel(),
                       data[:, 0, :].ravel())[0, 1] < 0:
            data[:, i, :] = -data[:, i, :]
    if all(data.mean((0, 2)) < 0):
        data = -data
    return data


@memory.cache
def setup_filters(info, forward, data_cov, noise_cov, labels,
                  reg=0.05, pick_ori='max-power', njobs=4):
    """Construct LCMV filters for ROIs.

    Args:
        info: MNE info structure
        forward: MNE forward solution
        data_cov: MNE Data covariance matrix
        noise_cov: MNE Noise covariance matrix
        labels: list
            A list of MNE label objects for which to construct
            filters.
        reg: float
            Regularization filter passed on to mne.make_filter
        pick_ori: str
            mne pick_ori argument of mne.make_filter
        njobs: int
            Number of cores to use.

    Returns:
        A dictionary that maps ROIs to filters
    """
    logging.info('Getting filter')
    tasks = []
    for l in labels:
        tasks.append(delayed(get_filter)(
            info, forward, data_cov,
            noise_cov, label=l, reg=reg,
            pick_ori='max-power'))

    filters = Parallel(n_jobs=njobs, verbose=1)(tasks)
    return {name: f for name, f in filters}


def reconstruct_tfr(
        filters, info, epochs,  events, times,
        estimator=complex_tfr, est_args=default_tfr,
        post_func=tfr2power_estimator, accumulate_func=accumulate,
        njobs=4):
    '''Reconstruct time frequency representation of epochs.

    Parallelization is applied across filters, i.e. regions of interest. This
    function calls par_reconstruct with appropriate default settings for TFR
    reconstruction. Change the est_args argument to specify parameters
    for the time frequency conversion.

    Args:
        filters: Dictionary returned by setup_filters
        info: MNE info object
        epochs: ndarray
            Data array of MNE epochs object
        events: ndarray
            Vector that assigns unique identifiers to different epochs
            in data array.
        times: ndarray
            Vector that assigns time points to time dimension in data
            array
        estimator: function
            A function that is applied to the sensor space data
            before  source reconstruction. Use 'complex_tfr' to project
            the  sensor space TFR representation into source space.
        est_args: dict
            A dict that is **pased to estimator, e.g. parameters for the
            TFR transformation.
        post_func: function
            A function that is applied to the source reconstructed data.
            To get TFR in source space you can pass 'complex_tfr' as
            'estimator' and then use 'tfr2power_estimator' here to
            compute power from complex FFT output.
        accumulate_func: function
            A function that takes the output of post func, the estimation keys,
            estimation values, time points, the region of interest and trial
            identifiers as inputs and returns a pandas DataFrame.
        njobs: int
            Number of cores to parallelize over.
    Returns:
        Concatenated outputs across regions of interest.
    '''
    M = par_reconstruct(
        pre_estimator=estimator, pre_est_args=est_args, epochs=epochs,
        events=events, times=times, info=info, filters=filters,
        post_func=post_func, accumulate_func=accumulate_func, njobs=njobs)
    return pd.concat([pd.concat(m, axis=0) for m in M if len(m) > 0], axis=1)
    # return pd.concat(M, axis=1)


def reconstruct_broadband(
        filters, info, epochs,  events, times,
        estimator=broadband_est, est_args={}, post_func=None,
        accumulate_func=accumulate, njobs=4):
    '''Reconstruct broadband activity from a set of regions of interest.

    Parallelization is applied across filters, i.e. regions of interest. This
    function calls par_reconstruct with appropriate default settings for
    broadband reconstruction.

    See reconstruct_tfr for description of arguments
    '''
    M = par_reconstruct(
        pre_estimator=estimator, pre_est_args=est_args, epochs=epochs,
        events=events, times=times, info=info, filters=filters,
        post_func=None, accumulate_func=accumulate_func, njobs=njobs)

    return pd.concat([pd.concat(m, axis=0) for m in M if len(m) > 0], axis=1)


def par_reconstruct(pre_estimator, pre_est_args, epochs, events, times,
                    info, filters, post_func=tfr2power_estimator,
                    accumulate_func=accumulate, njobs=4):
    '''Source reconstruct epochs with flexible transform before and after.

    This function performs source reconstruction and can transform the
    input data before reconstruction (e.g. for TFR) and after
    reconstruction (e.g. to compute power from complex FFT output). Output data
    can be  passed through yet another function to shape into the desired
    output.

    The flow of data through this function therefore is:

        epochs -> pre_estimator -> to source space ->
            post_func -> accumulate_func

    Args:
        pre_estimator: function
            A function that is applied to the sensor space data
            before  source reconstruction. Use 'complex_tfr' to project
            the  sensor space TFR representation into source space.
        pre_est_args: dict
            A dict that is **pased to pre_estimator, e.g. additional arguments
            to customize behavior of this function.
        epochs: ndarray
            Data array of MNE epochs object
        events: ndarray
            Vector that assigns unique identifiers to different epochs
            in data array.
        times: ndarray
            Vector that assigns time points to time dimension in data
            array
        info: MNE info object
        filters: Dictionary returned by setup_filters
        post_func: function
            A function that is applied to the source reconstructed data.
            To get TFR in source space you can pass 'complex_tfr' as
            'estimator' and then use 'tfr2power_estimator' here to
            compute power from complex FFT output.
        accumulate_func: function
            A function that takes the output of post func, the estimation keys,
            estimation values, time points, the region of interest and trial
            identifiers as inputs and returns a pandas DataFrame.
        njobs: int
            Number of cores to parallelize over.
    Returns:
        List that contains output for each ROI.
    '''
    pre_est_args['n_jobs'] = njobs
    logging.info('Applying pre-estimator function, params: ' +
                 str(pre_est_args))
    tfrdata, times, est_val, est_key = pre_estimator(
        epochs, times, **pre_est_args)
    logging.info('Done with pre-estimator. Data has shape ' +
                 str(tfrdata.shape) + ' now')
    tasks = []

    for filter in filters.keys():
        tasks.append(delayed(apply_lcmv)(
            tfrdata, est_key, est_val, events,
            times, info,
            {filter: filters[filter]}, post_func=post_func,
            accumulate_func=accumulate_func))
    logging.info(
        'Prepared %i tasks for parallel execution with %i jobs' %
        (len(tasks), njobs)
    )
    return Parallel(n_jobs=njobs, verbose=1)(tasks)


def apply_lcmv(tfrdata, est_key, est_vals, events, times, info,
               filters, post_func=None, accumulate_func=None,
               max_ori_out='signed'):
    """Apply Linearly Constrained Minimum Variance (LCMV) beamformer weights.


    Args:    
        tfrdata: ndarray
            Data to be reconstructed.
            Should be either n_trials x n_sensors x Y x n_time
            or trials x sensors x time. Reconstruction treats epochs and
            dim Y as independent dimensions.
        est_key: value
            A key to identify this reconstruction (e.g. F for power)
        est_vals: sequence
            Values that identify different reconstructions along dimension Y
            for a single epoch, e.g. the frequency for power reconstructions.
            Needs to be length Y.
        events: array
            Identifiers for different epochs. Needs to be of length n_trials.
        times: array
            Time of entries in last dimension of input data.
        info: mne info structure
            Info structure of the epochs which are to be reconstructed
        filters: dict
            Contains ROI names as keys and MNE filter dicts as values.
        post_func: function
            This function is applied to the reconstructed epochs, useful
            to convert complex TFR estimates into power values.
        accumulate_func: function
            Function that is applied after post_func has been applied.
            Can for example be used to transform the output into a dataframe.
        max_ori_out: str, default 'signed'
            This is passed to the MNE LCMV function which at the moment
            requires this to be 'signed'

    Returns:
        List of source reconstructed epochs transformed by post_func.
    """

    if accumulate_func is None:
        accumulate_func = lambda x: x
    if tfrdata.ndim == 3:
        # Must be trials x n_sensors x t_time
        tfrdata = tfrdata[:, :, np.newaxis, :]
    nfreqs = tfrdata.shape[2]
    assert(len(est_vals) == nfreqs)
    # ntrials = tfrdata.shape[0]
    info['sfreq'] = 1. / np.diff(times)[0]
    results = []
    for freq, (roi, filter) in product(range(nfreqs), filters.items()):
        if filter['weights'].size > 0:
            data = np.stack([x._data for x in
                             _apply_lcmv(data=tfrdata[:, :, freq, :],
                                         filters=filter,
                                         info=info, tmin=times.min(),
                                         max_ori_out=max_ori_out)])
            if post_func is None:
                results.append(accumulate_func(
                    data, est_key=est_key, time=times, est_val=est_vals[freq],
                    roi=roi, trial=events))
            else:
                data = post_func(data)
                results.append(accumulate_func(
                    data, est_key=est_key, time=times,
                    est_val=est_vals[freq], roi=roi, trial=events))
    return results


@memory.cache
def get_cov(epochs, tmin=0, tmax=1):
    """Compute a covariance matrix with default settings.

    This is mainly a helper function to cache computation of covariance
    matrices.
    """
    return compute_covariance(epochs, tmin=tmin, tmax=tmax, method='shrunk')


@memory.cache
def get_noise_cov(epochs, tmin=-0.5, tmax=0):
    """Compute a noise covariance matrix with default settings.

    This is mainly a helper function to cache computation of covariance
    matrices.
    """
    return compute_covariance(epochs, tmin=tmin, tmax=tmax, method='shrunk')


def get_filter(info, forward, data_cov, noise_cov, label=None, reg=0.05,
               pick_ori='max-power'):
    """Comput LCMV filter for one region of interest."""
    filter = make_lcmv(info=info,
                       forward=forward,
                       data_cov=data_cov,
                       noise_cov=noise_cov,
                       reg=0.05,
                       pick_ori='max-power',
                       label=label)
    del filter['data_cov']
    del filter['noise_cov']
    del filter['src']
    return label.name, filter


def get_filters(estimator, epochs, forward, source, noise_cov, data_cov,
                labels):
    """Compute LCMV filters for a list of regions of interest."""
    return {l.name: get_filter(epochs.info, forward, data_cov, noise_cov,
                               label=l, reg=0.05, pick_ori='max-power')
            for l in labels}


def _apply_lcmv(data, filters, info, tmin, max_ori_out):
    """Apply LCMV spatial filter to data for source reconstruction.
    
    Copied directly from MNE to remove dependence on source space in
    filter. This makes the filter much smaller and easier to use 
    multiprocessing here.

    Original authors: 
    Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
              Roman Goj <roman.goj@gmail.com>
              Britta Westner <britta.wstnr@gmail.com>

    Original License: BSD (3-clause)
    """
    from mne.source_estimate import _make_stc
    from mne.minimum_norm.inverse import combine_xyz
    if max_ori_out != 'signed':
        raise ValueError('max_ori_out must be "signed", got %s'
                         % (max_ori_out,))

    if isinstance(data, np.ndarray) and data.ndim == 2:
        data = [data]
        return_single = True
    else:
        return_single = False

    W = filters['weights']

    #subject = _subject_from_forward(filters)
    for i, M in enumerate(data):
        if len(M) != len(filters['ch_names']):
            raise ValueError('data and picks must have the same length')

        if filters['is_ssp']:
            raise RuntimeError('SSP not supported here')

        if filters['whitener'] is not None:
            M = np.dot(filters['whitener'], M)

        # project to source space using beamformer weights
        vector = False
        if filters['is_free_ori']:
            sol = np.dot(W, M)
            if filters['pick_ori'] == 'vector':
                vector = True
            else:
                sol = combine_xyz(sol)
        else:
            # Linear inverse: do computation here or delayed
            if (M.shape[0] < W.shape[0] and
                    filters['pick_ori'] != 'max-power'):
                sol = (W, M)
            else:
                sol = np.dot(W, M)
            if filters['pick_ori'] == 'max-power' and max_ori_out == 'abs':
                sol = np.abs(sol)

        tstep = 1.0 / info['sfreq']
        yield _make_stc(sol, vertices=filters['vertices'], tmin=tmin,
                        tstep=tstep, subject='NN', vector=vector,
                        source_nn=filters['source_nn'])
