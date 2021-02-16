import os
import pandas as pd
from glob import glob
import numpy as np
from joblib import Parallel, delayed
from joblib import Memory
import logging

from pymeg import atlas_glasser

memory = Memory(cachedir=os.environ['PYMEG_CACHE_DIR'], verbose=0)

backend = 'loky'


class Cache(object):
    """A cache that can prevent reloading from disk.

    Can be used as a context manager.
    """

    def __init__(self, cache=True):
        self.store = {}
        self.cache = cache

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.clear()

    def get(self, globstring):
        if self.cache:
            if globstring not in self.store:
                self.store[globstring] = self._load_tfr_data(globstring)
            else:
                logging.info('Returning cached object: %s' % globstring)
            return self.store[globstring]
        else:
            return self._load_tfr_data(globstring)

    def clear(self):
        self.cache = {}

    def _load_tfr_data(self, globstring):
        """Load all files identified by glob string"""
        logging.info('Loading data for: %s' % globstring)
        tfr_data_filenames = glob(globstring)
        logging.info('This is these filenames: %s' % str(tfr_data_filenames))
        tfrs = []
        for f in tfr_data_filenames:
            logging.info('Now working on: %s' % str(f))
            tfr = pd.read_hdf(f)
            logging.info('Done loading, now pivoting.')
            tfr = pd.pivot_table(tfr.reset_index(), values=tfr.columns, index=[
                                 'trial', 'est_val'], columns='time').stack(-2)
            tfr.index.names = ['trial', 'freq', 'area']
            tfrs.append(tfr)
        logging.info('Concate now.')
        tfr = pd.concat(tfrs)
        logging.info('Done _load_tfr_data.')
        return tfr


def baseline_per_sensor_get(tfr, baseline_time=(-0.25, -0)):
    '''
    Get average baseline
    '''
    time = tfr.columns.get_level_values('time').values.astype(float)
    id_base = (time >= baseline_time[0]) & (time <= baseline_time[1])
    base = tfr.loc[:, id_base].groupby(['freq', 'area']).mean().mean(
        axis=1)  # This should be len(nr_freqs * nr_hannels)
    return base


def baseline_per_sensor_apply(tfr, baseline):
    '''
    Baseline correction by dividing by average baseline
    '''
    def div(x):
        freqs = x.index.get_level_values('freq').values[0]
        areas = x.index.get_level_values('area').values[0]
        bval = float(baseline
                     .loc[
                         baseline.index.isin([freqs], level='freq')
                         & baseline.index.isin([areas], level='area')])
        return (x - bval) / bval * 100
    return tfr.groupby(['freq', 'area']).apply(div)


@memory.cache(ignore=['cache'])
def load_tfr_contrast(data_globstring, base_globstring, meta_data, conditions,
                      baseline_time, n_jobs=1, baseline_per_condition=True,
                      cache=Cache(cache=False)):
    """Load a set of data files and turn them into contrasts.
    """
    tfrs = []
    # load data:
    tfr_data = cache.get(data_globstring)
    # Make sure that meta_data and tfr_data overlap in trials
    tfr_trials = np.unique(tfr_data.index.get_level_values('trial').values)
    meta_trials = np.unique(meta_data.reset_index().loc[:, 'hash'].values)
    assert(any([t in meta_trials for t in tfr_trials]))

    # data to baseline:
    if not (data_globstring == base_globstring):
        tfr_data_to_baseline = cache.get(base_globstring)
    else:
        tfr_data_to_baseline = tfr_data

    if baseline_per_condition:
        # apply condition ind, collapse across trials, and get baseline::
        tfr_data_to_baseline = tfr_data_to_baseline.groupby(
            ['freq', 'area']).mean()
    # compute contrasts
    tasks = []
    for condition in conditions:
        tasks.append((tfr_data, tfr_data_to_baseline, meta_data,
                      condition, baseline_time, baseline_per_condition))

    tfr_conditions = Parallel(n_jobs=n_jobs, verbose=1, backend=backend)(
        delayed(make_tfr_contrasts)(*task) for task in tasks)
    weight_dicts = [t[1] for t in tfr_conditions]
    weights = weight_dicts.pop()
    [weights.update(w) for w in weight_dicts]
    # weights = {(k, v) for k, v in [t[1] for t in tfr_conditions]}
    tfrs.append(pd.concat([t[0] for t in tfr_conditions if t[0] is not None]))
    tfrs = pd.concat(tfrs)
    return tfrs, weights


def make_tfr_contrasts(tfr_data, tfr_data_to_baseline, meta_data,
                       condition, baseline_time, baseline_per_condition=True):

    # unpack:
    condition_ind = meta_data.loc[meta_data[condition] == 1, "hash"]

    if baseline_per_condition:
        # apply condition ind, collapse across trials, and get baseline::
        tfr_data_to_baseline = (tfr_data_to_baseline.loc[
            tfr_data_to_baseline.index.isin(condition_ind, level='trial'), :]
            .groupby(['freq', 'area']).mean())

    baseline = baseline_per_sensor_get(
        tfr_data_to_baseline, baseline_time=baseline_time)

    # apply condition ind, and collapse across trials:
    tfr_data_condition = (tfr_data.loc[
        tfr_data.index.isin(condition_ind, level='trial'), :])
    num_trials_in_condition = len(np.unique(
        tfr_data_condition.index.get_level_values('trial')))
    if num_trials_in_condition == 0:
        return None, {condition: num_trials_in_condition}
    tfr_data_condition = tfr_data_condition.groupby(['freq', 'area']).mean()

    # apply baseline, and collapse across sensors:
    tfr_data_condition = baseline_per_sensor_apply(
        tfr_data_condition, baseline=baseline).groupby(['freq', 'area']).mean()

    tfr_data_condition['condition'] = condition
    tfr_data_condition = tfr_data_condition.set_index(
        ['condition', ], append=True, inplace=False)
    tfr_data_condition = tfr_data_condition.reorder_levels(
        ['area', 'condition', 'freq'])
    return tfr_data_condition, {condition: num_trials_in_condition}


@memory.cache(ignore=['cache'])
def single_conditions(conditions, data_glob, base_glob, meta_data,
                      baseline_time, baseline_per_condition=False,
                      n_jobs=1, cache=Cache(cache=False)):

    tfr_condition, weights = load_tfr_contrast(
        data_glob, base_glob, meta_data,
        list(conditions), baseline_time, n_jobs=n_jobs,
        baseline_per_condition=baseline_per_condition,
        cache=cache)
    return tfr_condition.groupby(
        ['area', 'condition', 'freq']).mean(), weights


@memory.cache(ignore=['cache'])
def pool_conditions(conditions, data_globs, base_globs, meta_data,
                    baseline_time, baseline_per_condition=False,
                    n_jobs=1, cache=Cache(cache=False)):
    weights = {}
    tfrs = {}
    for i, (data_glob, base_glob) in enumerate(
            zip(ensure_iter(data_globs), ensure_iter(base_globs))):
        # tfr, weight = single_conditions(
        #    conditions, data_glob, base_glob, meta_data, baseline_time,
        #    n_jobs=n_jobs,
        #    cache=cache)
        tfr, weight = load_tfr_contrast(
            data_glob, base_glob, meta_data,
            list(conditions), baseline_time, n_jobs=n_jobs,
            baseline_per_condition=baseline_per_condition,
            cache=cache)
        tfrs[i] = tfr
        weights[i] = weight
        # Compute total trials per condition
    total_weights = {}
    for i, w in weights.items():
        for k, v in w.items():
            if k not in total_weights:
                total_weights[k] = v
            else:
                total_weights[k] += v
    # Apply weights to each tfr
    ind_weights = {}
    for k in total_weights.keys():
        ind_weights[k] = []
    for key in tfrs.keys():
        tfr = tfrs[key]
        for condition in total_weights.keys():
            condition_ind = tfr.index.get_level_values(
                'condition') == condition
            w = weights[key][condition] / total_weights[condition]
            tfr.loc[condition_ind, :] *= w
            ind_weights[condition].append(w)
        tfrs[key] = tfr
    for condition, weights in ind_weights.items():
        logging.info("weights for %s -> %s, sum=%f" %
                     (condition, str(weights), sum(weights)))
    tfrs = pd.concat(tfrs.values()).groupby(
        ['area', 'condition', 'freq']).sum()
    return tfrs


@memory.cache(ignore=['cache'])
def compute_contrast(contrasts, hemis, data_globstring, base_globstring,
                     meta_data, baseline_time, baseline_per_condition=True,
                     n_jobs=1, cache=Cache(cache=False)):
    """Compute a single contrast from tfr data
    Args:
        contrast: dict
            Contains contrast names as keys and len==2 tuples as values. The
            tuples contain a list of condition names first and then a set of
            weights for each condition. Condition names identify columns in
            the meta data that are one for each trial that belongs to
            this condition.
        hemi: str
            Can be:
                'lh_is_ipsi' if contrast is ipsi-contra hemi and left hemi is
                    ipsi.
                'rh_is_ipsi' if contrast is ipis-contra and right hemi is ipsi
                'avg' if contrast should be averaged across hemispheres
        data_globstring: list
            Each string in data_globstring selects a set of filenames if
            passed through glob. Condition averages and baselines are then
            computed for each group of filenames identified by one entry
            in data_globstring. This is useful for, e.g. computing
            conditions per session first, then averaging them and then
            computing contrasts across sessions.
        base_globstring: string or list
            Same as data_globstring but selects data to use for baselining
        meta_data: data frame
            Meta data DataFrame with as many rows as trials.
        baseline_time: tuple

    """

    # load for all subjects:
    tfr_condition = []
    from functools import reduce
    from itertools import product

    conditions = set(
        reduce(lambda x, y: x + y, [x[0] for x in contrasts.values()]))

    tfr_condition = pool_conditions(conditions, data_globstring,
                                    base_globstring, meta_data,
                                    baseline_time, n_jobs=n_jobs,
                                    cache=cache)

    # Lower case all area names
    # FIXME: Set all area names to lower case!
    all_clusters, _, _, _ = atlas_glasser.get_clusters()
    tfr_areas = np.array([a for a in tfr_condition.index.levels[
        np.where(np.array(tfr_condition.index.names) == 'area')[0][0]]])
    tfr_areas_lower = np.array([area.lower() for area in tfr_areas])
    for cluster, areas in all_clusters.items():
        new_areas = []
        for area in areas:
            idx = np.where(tfr_areas_lower == area.lower())[0]
            if len(idx) == 1:
                new_areas.append(tfr_areas[idx[0]])
        all_clusters[cluster] = new_areas
    # mean across sessions:
    tfr_condition = tfr_condition.groupby(
        ['area', 'condition', 'freq']).mean()
    cluster_contrasts = []
    for cur_contrast, hemi, cluster in product(contrasts.items(), hemis,
                                               all_clusters.keys()):
        contrast, (conditions, weights) = cur_contrast
        logging.info('Start computing contrast %s for cluster %s' %
                     (contrast, cluster))
        right = []
        left = []
        for condition in conditions:
            tfrs_rh = []
            tfrs_lh = []
            for area in all_clusters[cluster]:
                area_idx = tfr_condition.index.isin([area], level='area')
                condition_idx = tfr_condition.index.isin(
                    [condition], level='condition')
                subset = tfr_condition.loc[area_idx & condition_idx].groupby(
                    ['freq']).mean()
                if 'rh' in area:
                    tfrs_rh.append(subset)
                else:
                    tfrs_lh.append(subset)
            # What happens when an area is not defined for both hemis?
            if (len(tfrs_lh) == 0) and (len(tfrs_rh) == 0):
                logging.warn('Skipping condition %s in cluster %s' %
                             (condition, cluster))
                continue
            try:
                left.append(pd.concat(tfrs_lh))
            except ValueError:
                pass
            try:
                right.append(pd.concat(tfrs_rh))
            except ValueError:
                pass

        if (len(left) == 0) and (len(right) == 0):
            logging.warn('Skipping cluster %s' % (cluster))
            continue
        if hemi == 'rh_is_ipsi':
            left, right = right, left
        if 'is_ipsi' in hemi:
            if not len(left) == len(right):
                logging.warn('Skipping cluster %s: does not have the same number of lh/rh rois' %
                             (cluster))
                continue
            tfrs = [left[i] - right[i]
                    for i in range(len(left))]
        else:
            if (len(right) == 0) and (len(left) == len(weights)):
                tfrs = left
            elif (len(left) == 0) and (len(right) == len(weights)):
                tfrs = right
            else:
                tfrs = [(right[i] + left[i]) / 2
                        for i in range(len(left))]
        assert(len(tfrs) == len(weights))
        tfrs = [tfr * weight for tfr, weight in zip(tfrs, weights)]
        tfrs = reduce(lambda x, y: x + y, tfrs)
        tfrs = tfrs.groupby('freq').mean()
        tfrs.loc[:, 'cluster'] = cluster
        tfrs.loc[:, 'contrast'] = contrast
        tfrs.loc[:, 'hemi'] = hemi
        cluster_contrasts.append(tfrs)
    logging.info('Done compute contrast')
    return pd.concat(cluster_contrasts)


def augment_data(meta, response_left, stimulus):
    """Augment meta data with fields for specific cases

    Args:
        meta: DataFrame
        response_left: ndarray
            1 if subject made a left_response / yes response
        stimulus: ndarray
            1 if a left_response is correct
    """
    # add columns:
    meta["all"] = 1

    meta["left"] = response_left.astype(int)
    meta["right"] = (~response_left).astype(int)

    meta["hit"] = ((response_left == 1) & (stimulus == 1)).astype(int)
    meta["fa"] = ((response_left == 1) & (stimulus == 0)).astype(int)
    meta["miss"] = ((response_left == 0) & (stimulus == 1)).astype(int)
    meta["cr"] = ((response_left == 0) & (stimulus == 0)).astype(int)
    return meta


def set_jw_style():
    import matplotlib
    import seaborn as sns
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    sns.set(style='ticks', font='Arial', font_scale=1, rc={
        'axes.linewidth': 0.05,
        'axes.labelsize': 7,
        'axes.titlesize': 7,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'legend.fontsize': 6,
        'xtick.major.width': 0.25,
        'xtick.minor.width': 0.25,
        'ytick.major.width': 0.25,
        'text.color': 'Black',
        'axes.labelcolor': 'Black',
        'xtick.color': 'Black',
        'ytick.color': 'Black', })
    sns.plotting_context()


def plot_mosaic(tfr_data, vmin=-25, vmax=25, cmap='RdBu_r',
                ncols=4, epoch='stimulus', stats=False,
                threshold=0.05):
    from mne.viz.utils import _plot_masked_image as pmi
    if epoch == "stimulus":
        time_cutoff = (-0.5, 1.35)
        xticks = [0, 0.25, 0.5, 0.75, 1]
        xticklabels = ['0\nStim on', '', '.5', '', '1\nStim off']
        yticks = [25, 50, 75, 100, 125]
        yticklabels = ['25', '', '75', '', '125']
        xmarker = [0, 1]
        baseline = (-0.25, 0)
    else:
        time_cutoff = (-1, .5)
        xticks = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5]
        xticklabels = ['-1', '', '-0.5', '', '0\nResponse', '', '0.5']
        yticks = [1, 25, 50, 75, 100, 125]
        yticklabels = ['1', '25', '', '75', '', '125']
        xmarker = [0, 1]
        baseline = None
    from matplotlib import gridspec
    import pylab as plt
    import seaborn as sns
    set_jw_style()
    sns.set_style('ticks')
    nrows = (len(atlas_glasser.areas) // ncols) + 1
    gs = gridspec.GridSpec(nrows, ncols)
    gs.update(wspace=0.01, hspace=0.01)

    for i, (name, area) in enumerate(atlas_glasser.areas.items()):
        try:
            column = np.mod(i, ncols)
            row = i // ncols
            plt.subplot(gs[row, column])
            times, freqs, tfr = get_tfr(tfr_data.query('cluster=="%s"' %
                                                       area), time_cutoff)
            # cax = plt.gca().pcolormesh(times, freqs, np.nanmean(
            #    tfr, 0), vmin=vmin, vmax=vmax, cmap=cmap, zorder=-2)
            mask = None
            if stats:
                _, _, cluster_p_values, _ = get_tfr_stats(
                    times, freqs, tfr, threshold)
                sig = cluster_p_values.reshape((tfr.shape[1], tfr.shape[2]))
                mask = sig < threshold
            cax = pmi(plt.gca(),  np.nanmean(tfr, 0), times, yvals=freqs,
                      yscale='linear', vmin=vmin, vmax=vmax,
                      mask=mask, mask_alpha=1,
                      mask_cmap=cmap, cmap=cmap)

            # plt.grid(True, alpha=0.5)
            for xmark in xmarker:
                plt.axvline(xmark, color='k', lw=1, zorder=-1, alpha=0.5)

            plt.yticks(yticks, [''] * len(yticks))
            plt.xticks(xticks, [''] * len(xticks))
            set_title(name, times, freqs, plt.gca())
            plt.tick_params(direction='inout', length=2, zorder=100)
            plt.xlim(time_cutoff)
            plt.ylim([1, 147.5])
            plt.axhline(10, color='k', lw=1, alpha=0.5, linestyle='--')
        except ValueError as e:
            print(name, area, e)
    plt.subplot(gs[nrows - 2, 0])

    sns.despine(left=True, bottom=True)
    plt.subplot(gs[nrows - 1, 0])

    pmi(plt.gca(),  np.nanmean(tfr, 0) * 0, times, yvals=freqs,
        yscale='linear', vmin=vmin, vmax=vmax,
        mask=None, mask_alpha=1,
        mask_cmap=cmap, cmap=cmap)
    plt.xticks(xticks, xticklabels)
    plt.yticks(yticks, yticklabels)
    for xmark in xmarker:
        plt.axvline(xmark, color='k', lw=1, zorder=-1, alpha=0.5)
    if baseline is not None:
        plt.fill_between(baseline, y1=[1, 1],
                         y2=[150, 150], color='k', alpha=0.5)
    plt.tick_params(direction='in', length=3)
    plt.xlim(time_cutoff)
    plt.ylim([1, 147.5])
    plt.xlabel('time [s]')
    plt.ylabel('Freq [Hz]')
    sns.despine(ax=plt.gca())


def plot_tfr(df, vmin=-5, vmax=5, cmap='RdBu_r', threshold=0.05):
    import pylab as plt
    from mne.viz.utils import _plot_masked_image as pmi
    times, freqs, tfr = get_tfr(df, (-np.inf, np.inf))
    T_obs, clusters, cluster_p_values, h0 = get_tfr_stats(
        times, freqs, tfr, 0.05)
    sig = cluster_p_values.reshape((tfr.shape[1], tfr.shape[2]))

    cax = pmi(plt.gca(), tfr.mean(0), times, yvals=freqs,
              yscale='linear', vmin=vmin, vmax=vmax, mask=sig < threshold,
              mask_alpha=1, mask_cmap=cmap, cmap=cmap, linewidth=1)

    return cax, times, freqs, tfr


@memory.cache()
def get_tfr_stats(times, freqs, tfr, threshold=0.05):
    from mne.stats import permutation_cluster_1samp_test as cluster_test
    return cluster_test(
        tfr, threshold={'start': 0, 'step': 0.2},
        connectivity=None, tail=0, n_permutations=1000, n_jobs=2)


def plot_tfr_stats(times, freqs, tfr, threshold=0.05):
    import pylab as plt
    from matplotlib.colors import LinearSegmentedColormap
    T_obs, clusters, cluster_p_values, h0 = get_tfr_stats(
        times, freqs, tfr, threshold)
    vmax = 1
    cmap = LinearSegmentedColormap.from_list('Pvals', [(0 / vmax, (1, 1, 1, 0)),
                                                       (0.04999 / vmax,
                                                        (1, 1, 1, 0)),
                                                       (0.05 / vmax,
                                                        (1, 1, 1, 0.5)),
                                                       (1 / vmax, (1, 1, 1, 0.5))]
                                             )
    sig = cluster_p_values.reshape((tfr.shape[1], tfr.shape[2]))

    # df = np.array(list(np.diff(freqs) / 2) + [freqs[-1] - freqs[-2]])
    # dt = np.array(list(np.diff(times) / 2) + [times[-1] - times[-2]])
    # from scipy.interpolate import interp2d
    # print(np.unique((sig < threshold).astype(float)))
    # i = interp2d(times, freqs, sig.astype(float))
    # X, Y = (np.linspace(times[0], times[-1], len(times) * 25),
    #        np.linspace(freqs[0], freqs[-1], len(freqs) * 25))
    # Z = i(X.ravel(), Y.ravel())

    # plt.gca().pcolormesh(times, freqs, sig, vmin=0, vmax=1, cmap=cmap)

    plt.gca().contour(
        times, freqs,
        sig, (threshold),
        linewidths=0.5, colors=('black'))
    return X, Y, Z


def set_title(text, times, freqs, axes):
    import pylab as plt
    x = np.mean(times)
    y = np.max(freqs)
    plt.text(x, y, text, fontsize=8,
             verticalalignment='top', horizontalalignment='center')


def get_tfr(tfr, time_cutoff):
    # variables:
    times = np.array(tfr.columns, dtype=float)
    freqs = np.array(
        np.unique(tfr.index.get_level_values('freq')), dtype=float)
    time_ind = (times > time_cutoff[0]) & (times < time_cutoff[1])
    time_ind = (times > time_cutoff[0]) & (times < time_cutoff[1])

    tfrs = [tfr.loc[tfr.index.isin([subj], level='subject'), time_ind].values
            for subj in np.unique(tfr.index.get_level_values('subject'))]
    # data:
    X = np.stack(tfrs)
    return times[time_ind], freqs, X


def plot_cluster(names, view):
    from pymeg import atlas_glasser
    all_clusters, _, _, _ = atlas_glasser.get_clusters()
    label_names = []
    for name in names:
        cluster_name = atlas_glasser.areas[name]
        label_names.extend(all_clusters[cluster_name])

    plot_roi('lh', label_names, 'r')


#@memory.cache
def plot_roi(hemi, labels, color, annotation='HCPMMP1',
             view='parietal',
             fs_dir=os.environ['SUBJECTS_DIR'],
             subject_id='S04', surf='inflated'):
    import matplotlib
    import os
    import glob
    from surfer import Brain
    from mne import Label
    color = np.array(matplotlib.colors.to_rgba(color))

    brain = Brain(subject_id, hemi, surf, offscreen=False)
    labels = [label.replace('-rh', '').replace('-lh', '') for label in labels]
    # First select all label files

    label_names = glob.glob(os.path.join(
        fs_dir, subject_id, 'label', 'lh*.label'))
    label_names = [label for label in label_names if any(
        [l in label for l in labels])]

    for label in label_names:
        brain.add_label(label, color=color)

    # Now go for annotations
    from nibabel.freesurfer import io
    ids, colors, annot_names = io.read_annot(os.path.join(
        fs_dir, subject_id, 'label', 'lh.%s.annot' % annotation),
        orig_ids=True)

    for i, alabel in enumerate(annot_names):
        if any([label in alabel.decode('utf-8') for label in labels]):
            label_id = colors[i, -1]
            vertices = np.where(ids == label_id)[0]
            l = Label(np.sort(vertices), hemi='lh')
            brain.add_label(l, color=color)
    brain.show_view(view)
    return brain.screenshot()


def ensure_iter(input):
    if isinstance(input, str):
        yield input
    else:
        try:
            for item in input:
                yield item
        except TypeError:
            yield input
