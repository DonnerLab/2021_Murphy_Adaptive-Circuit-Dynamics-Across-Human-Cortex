import logging
import mne
import numpy as np
import os

from joblib import Memory

from os import makedirs
from os.path import join
from glob import glob

from pymeg import lcmv as pymeglcmv
from pymeg import source_reconstruction as pymegsr


memory = Memory(cachedir=os.environ['PYMEG_CACHE_DIR'])
path = '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/'


def set_n_threads(n):
    import os
    os.environ['OPENBLAS_NUM_THREADS'] = str(n)
    os.environ['MKL_NUM_THREADS'] = str(n)
    os.environ['OMP_NUM_THREADS'] = str(n)

subjects = {'DCB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'DHB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1)],
            'ECB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'EMB': [(1, 1), (1, 2), (2, 1), (2, 3), (3, 1), (3, 2)],
            'EXF': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'EXG': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'GSB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'HBC': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'JTB': [(1, 1), (1, 2), (2, 1), (2, 3), (3, 1), (3, 2)],
            'KSV': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'NIF': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'OMF': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'PDP': [(1, 1), (1, 2), (2, 2), (2, 3), (3, 1), (3, 2)],
            'QNV': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2)],
            'TFD': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'TNB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)],
            'TSJ': [(1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2)]}


def submit():
    from pymeg import parallel
    for subject, tasks in subjects.items():
        for session, recording in tasks:
            for signal in ['BB']: # for signal in ['BB', 'HF', 'LF']:
                parallel.pmap(
                    extract, [(subject, session, recording, signal)],
                    walltime='15:00:00', memory=50, nodes=1, tasks=5,
                    name='sr' + str(subject) + '_' + str(session) + str(recording),
                    ssh_to=None, env='mne')

def lcmvfilename(subject, session, signal, recording, chunk=None):
    try:
        makedirs(path)
    except:
        pass
    if chunk is None:
        filename = '%s-SESS%i-%i-%s-lcmv.hdf' % (
            subject, session, recording, signal)
    else:
        filename = '%s-SESS%i-%i-%s-chunk%i-lcmv.hdf' % (
            subject, session, recording, signal, chunk)
    return join(path, filename)


def get_stim_epoch(subject, session, recording, hipass=None, lopass=None):
    from pymeg import preprocessing as pymegprepr
    globstring = '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/%s-%i_0%i_*fif.gz' % (
        subject, session, recording)
    filenames = glob(globstring)[0]
    epochs = mne.read_epochs(filenames)
    epochs.times = epochs.times - 1  # PM: this was somehow necessary for initial Conv2mne pipeline, but *NOT* for induced
    epochs = epochs.pick_channels(
        [x for x in epochs.ch_names if x.startswith('M')])
    if (hipass is not None) or (lopass is not None):  # filter epoched data if desired
        epochs = epochs.filter(hipass,lopass)
    id_time = (-0.25 <= epochs.times) & (epochs.times <= 0)
    means = epochs._data[:, :, id_time].mean(-1)
    epochs._data -= means[:, :, np.newaxis]
    min_time, max_time = epochs.times.min() + 0.75, epochs.times.max() - 0.75
    data_cov = pymeglcmv.get_cov(epochs, tmin=min_time, tmax=max_time)
    return data_cov, epochs, filenames


def extract(subject, session, recording, signal_type='BB',
            BEM='three_layer', debug=False, chunks=100, njobs=4):
    mne.set_log_level('WARNING')
    pymeglcmv.logging.getLogger().setLevel(logging.INFO)
    set_n_threads(1)

    logging.info('Reading stimulus data')
    data_cov, epochs, epochs_filename = get_stim_epoch(
        subject, session, recording)

    raw_filename = glob('/home/pmurphy/meg_data/surprise/%s-%i*_0%i.ds' %
                        (subject, session, recording))
    assert len(raw_filename) == 1
    raw_filename = raw_filename[0]

    trans_filename = glob('/home/pmurphy/meg_data/surprise/MRIs/trans_mats/%s_%i_0%i*fif' % (
        subject, session, recording))[0]
    logging.info('Setting up source space and forward model')

    forward, bem, source = pymegsr.get_leadfield(
        subject, raw_filename, epochs_filename, trans_filename, bem_sub_path='bem_ft')
    labels = pymegsr.get_labels(subject)
    labels = pymegsr.labels_exclude(labels,
                                    exclude_filters=['wang2015atlas.IPS4',
                                                     'wang2015atlas.IPS5',
                                                     'wang2015atlas.SPL',
                                                     'JWDG_lat_Unknown'])
    labels = pymegsr.labels_remove_overlap(
        labels, priority_filters=['wang', 'JWDG'],)

    # Now chunk Reconstruction into blocks of ~100 trials to save Memory
    fois_h = np.arange(36, 162, 4)
    fois_l = np.arange(1, 36, 1)
    tfr_params = {
        'HF': {'foi': fois_h, 'cycles': fois_h * 0.25, 'time_bandwidth': 6 + 1,
               'n_jobs': njobs, 'est_val': fois_h, 'est_key': 'HF', 'sf': 400,
               'decim': 20},
        'LF': {'foi': fois_l, 'cycles': fois_l * 0.4, 'time_bandwidth': 1 + 1,
               'n_jobs': njobs, 'est_val': fois_l, 'est_key': 'LF', 'sf': 400,
               'decim': 20}
    }

    events = epochs.events[:, 2]
    filters = pymeglcmv.setup_filters(epochs.info, forward, data_cov, None, labels)

    set_n_threads(1)

    for i in range(0, len(events), chunks):
        filename = lcmvfilename(
            subject, session, signal_type, recording, chunk=i)
        if os.path.isfile(filename):
            continue
        if signal_type == 'BB':
            logging.info('Starting reconstruction of BB signal')
            M = pymeglcmv.reconstruct_broadband(
                filters, epochs.info, epochs._data[i:i + chunks],
                events[i:i + chunks],
                epochs.times, njobs=1)
        else:
            logging.info('Starting reconstruction of TFR signal')
            M = pymeglcmv.reconstruct_tfr(
                filters, epochs.info, epochs._data[i:i + chunks],
                events[i:i + chunks], epochs.times,
                est_args=tfr_params[signal_type],
                njobs=4)
        M.to_hdf(filename, 'epochs')
    set_n_threads(njobs)
