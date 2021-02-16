"""
Aggregate source reconstructed files across chunks and sessions into
one unified format.

Aggregation will combine individual ROIs into larger clusters more
suitable for MEG data.

The source reconstruction process outputs data frames that are
indexed by time, est_val, est_type, epoch etc. in the rows
and ROIs in the columns:

                                 lh.wang2015atlas.V1v-lh ...
trial time      est_key est_val
506   -1.500000 F       10                  2.206130e-25 ...
      -1.483333 F       10                  3.374152e-25 ...
      -1.466667 F       10                  4.967676e-25 ...
      -1.450000 F       10                  6.293999e-25 ...
      -1.433333 F       10                  6.862688e-25 ...


After aggregation all chunks from one session will be combined. ROIs
will be aggregated into clusters. baseline corrected and converted to
percent signal change. Cluster can be averaged across hemis,
be lateralized or kept as is. The resulting structure is:

time                                               -0.750000 ...
hemi     cluster                       trial freq
Averaged HCPMMP1_audiotory_association 503   10   -73.221365 ...
                                       504   10   -66.933821 ...
                                       505   10   -64.982795 ...
                                       506   10   -69.250634 ...
                                       507   10   -35.822782 ...

Aggregated files can be saved in matlab readable HDF files. These
files are organized as follows:

    /Averaged
        /Cluster 1
            - freq A: 2D dataset (trial x time)
            - freq B ...
        /Cluster 2
            ...
    /Lateralized
        ...

Indices into the rows and columns of individual datasets are stored
in their attributes.

"""
import pandas as pd
import numpy as np


def aggregate_files(data_globstring, base_globstring, baseline_time,
                    hemis=['Averaged', 'Pair', 'Lateralized'],
                    cache=None, to_decibels=False):
    """Read source reconstructed files, baseline correct and aggregate into
    area clusters.

    Args:
        data_globstring: globstring that selects data files
        base_globstring: globstring that selects baseline files
        baseline_time: 2-tuple
            Defines time to use for baselining
        hemis: List of strings
            Can contain the following:
                'Averaged': Mean over hemispheres per area
                'Pair': Return areas as they are per hemisphere
                'Lateralized': Subtract left from right hemisphere
        to_decibels: bool
            Convert data to decibels - implies subtracting log baseline
            but no conversion to percent signal change.
    Returns:
        DataFrame that contains time points as columns and is indexed
        by time, frequency and cluster in the row index.
    """
    from pymeg.contrast_tfr import Cache
    tfr_baseline = None
    if not (data_globstring == base_globstring):
        with Cache() as base_cache:
            tfr_baseline = base_cache.get(base_globstring)
            if to_decibels:
                tfr_baseline = 10 * np.log10(tfr_baseline)
            tfr_baseline = tfr_baseline.groupby(['freq', 'area']).mean()

    if cache is None:
        cache = Cache()
    tfr_data = cache.get(data_globstring)
    if to_decibels:
        tfr_data = 10 * np.log10(tfr_data)
    if tfr_baseline is None:
        tfr_baseline = tfr_data.groupby(['freq', 'area']).mean()

    baseline = tfr_baseline.loc[:, slice(*baseline_time)].mean(1)
    baseline.name = 'baseline'
    cols = tfr_data.columns
    tfr_data = tfr_data.join(baseline, on=['freq', 'area'])
    baseline = tfr_data.loc[:, 'baseline']
    if to_decibels:
        tfr_data = tfr_data.loc[:, cols].sub(baseline, axis=0)
    else:
        tfr_data = ((tfr_data.loc[:, cols].sub(
            baseline, axis=0)).div(baseline, axis=0)) * 100
    aggs = aggregate(tfr_data, hemis)
    return aggs


def agg2hdf(agg, filename):
    """Convert an aggregate into a HDF file.

    The resulting HDF file encodes the hemi and cluster
    index hierarchically as groups. Each frequency is a dataset
    that is itself 2D: trials x time.

    Indices into the 2D datasets are saved in the datasets attrs.
    Row indices are prefixed with 'rows_' and column indices with
    'cols_'.

    Args:
        agg: DataFrame
        filename: str
            Path to file.
    """
    import h5py
    with h5py.File(filename, mode='w') as store:
        for (hemi, cluster, freq), data in agg.groupby(['hemi', 'cluster', 'freq']):
            try:
                grp = store.create_group(hemi + '/' + cluster)
            except ValueError:
                grp = store[hemi + '/' + cluster]
            dset = grp.create_dataset(
                str(freq), data=data.values, compression="gzip", compression_opts=7)
            dset.attrs['cols_time'] = data.columns.values.astype(float)
            for index in data.index.names:
                index_vals = data.index.get_level_values(index).values
                if index_vals.dtype == object:
                    index_vals = [str(i).encode('utf-8') for i in index_vals]
                dset.attrs[
                    'rows_' + index] = index_vals


def delayed_agg(filename, hemi=None, cluster=None, freq=None):
    from functools import partial
    return partial(hdf2agg, filename, hemi=hemi, cluster=cluster, freq=freq)


def hdf2agg(filenames, hemi=None, cluster=None, freq=None):
    return pd.concat([_hdf2agg(f, hemi, cluster, freq)
                      for f in ensure_iter(filenames)])


def _hdf2agg(filename, hemi=None, cluster=None, freq=None):
    """Convert HDF file back to aggregate DataFrame.

    Args:
        filename: Path to aggregate file
        hemi: str, default None
            Restrict return to these hemi combination strategy.
        cluster: str, default None
            Restrict return to this cluster
        freq: int, default None
            Restrict return to this frequency
    Returns:
        DataFrame indexed by hemi, cluster, trial and freq in the rows
        and time in the columns.
    """
    import h5py
    dfs = []
    with h5py.File(filename, mode='r') as store:
        for fhemi, hd in store.items():
            if hemi is not None and not (str(hemi) == fhemi):
                continue
            for fcluster, cd in hd.items():
                if (cluster is not None) and (
                        not any([str(c) == fcluster for c in
                                 ensure_iter(cluster)])):
                    continue
                for fF, Fd in cd.items():
                    if freq is not None and not (str(freq) == fF):
                        continue
                    dfs.append(get_df_from_hdf(Fd))

    return pd.concat(dfs)


def get_df_from_hdf(dataset):
    """Convert HDF dataset to pandas Dataframe
    """
    cols = []
    rows = []
    row_names = []
    col_names = []
    for key, values in dataset.attrs.items():
        dname = values.dtype.name

        if 'bytes' in dname:
            values = values.astype('U')
            try:
                values = [float(i) for i in values]
            except ValueError:
                pass
        if key.startswith('cols_'):
            cols.append(values)
            col_names.append(key.replace('cols_', ''))
        if key.startswith('rows_'):
            rows.append(values)
            row_names.append(key.replace('rows_', ''))
    index = pd.MultiIndex.from_arrays(rows)
    index.names = row_names

    cols = pd.MultiIndex.from_arrays(cols)
    cols.names = col_names
    return pd.DataFrame(dataset[:], index=index, columns=cols)


def aggregate(tfr_data, hemis):
    """Aggregate individual areas into clusters.
    """
    from itertools import product
    from pymeg import atlas_glasser
    all_clusters, _, _, _ = atlas_glasser.get_clusters()
    clusters = []
    tfr_areas = np.unique(tfr_data.index.get_level_values('area'))
    for hemi, cluster in product(hemis, all_clusters.keys()):
        print('Working on %s, %s' % (hemi, cluster))
        tfrs_rh = [area for area in all_clusters[cluster] if 'rh' in area]
        tfrs_lh = [area for area in all_clusters[cluster] if 'lh' in area]
        tfrs_rh = [t for t in tfr_areas if any(
            [a.lower() in t.lower() for a in tfrs_rh])]
        tfrs_lh = [t for t in tfr_areas if any(
            [a.lower() in t.lower() for a in tfrs_lh])]
        lh_idx = tfr_data.index.isin(tfrs_lh, level='area')
        rh_idx = tfr_data.index.isin(tfrs_rh, level='area')
        left = tfr_data.loc[lh_idx, :].groupby(
            ['freq', 'trial']).mean()
        right = tfr_data.loc[rh_idx, :].groupby(
            ['freq', 'trial']).mean()

        if hemi == 'Pair':
            left.loc[:, 'cluster'] = cluster + '_LH'
            left.loc[:, 'hemi'] = 'Pair'
            right.loc[:, 'cluster'] = cluster + '_RH'
            right.loc[:, 'hemi'] = 'Pair'
            clusters.append(left)
            clusters.append(right)
        else:
            if hemi == 'Lateralized':
                tfrs = left - right
            elif hemi == 'Averaged':
                tfrs = (right + left) / 2
            tfrs.loc[:, 'cluster'] = cluster
            tfrs.loc[:, 'hemi'] = hemi
            clusters.append(tfrs)
    df = pd.concat(clusters)
    df.set_index(['cluster', 'hemi'], append=True, inplace=True)
    return df.reorder_levels(
        ['hemi', 'cluster', 'trial', 'freq'])


def ensure_iter(input):
    if isinstance(input, str):
        yield input
    else:
        try:
            for item in input:
                yield item
        except TypeError:
            yield input
