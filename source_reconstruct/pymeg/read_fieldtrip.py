
import h5py
import numpy as np
import os
import glob

from joblib import Memory

import mne

from mne import create_info
from mne.epochs import EpochsArray


memory = Memory(cachedir=os.environ['PYMEG_CACHE_DIR'], verbose=0)


def fix_chs(rawinfo, einfo):
    ch_names = []
    for k in range(len(einfo['chs'])):
        name = einfo['chs'][k]['ch_name']
        newchan = [x for x in rawinfo['chs']
                   if name in x['ch_name']][0]
        einfo['chs'][k] = newchan
        ch_names.append(newchan['ch_name'])
    einfo['ch_names'] = ch_names
    return einfo


@memory.cache
def get_info_for_epochs(rawname):
    raw = mne.io.ctf.read_raw_ctf(rawname)
    return raw.info


def read_ft_epochs(fname, rawinfo, cachedir=os.environ['PYMEG_CACHE_DIR'],
                   trialinfo_col=-1):
    '''
    Read and cache the output of fieldtrip epochs.

    This function reads a matlab file that contains a 'data' struct with the
    following fields:

        trialinfo: matrix
            Dim is ntrials x nr_meta, the columns contain meta information
            about each trial.
        label: list of strings
            Channel names
        sampleinfo: matrix
            Dim is ntrials x 2, the first column contains the start sample
            of each epoch in the raw data.
        time: array
            Contains time points for epochs.
        trial: array
            Dim is time x channels x trials, contains the actial data

    This data is parsed into an MNE Epochs object. To correctly assign
    channel locations, types etc. the info structure from the raw data
    that generated the fieldtrip epochs is used. The channel names in
    the fieldtrip structure should still be relatable to the raw
    channel names, relatable here means that a fieldtrip channel name
    must be contained in the raw channel name.

    Args
        fname: str
            Path to .mat file to load
        rawinfo: mne info structure
            Info structure with correct channel locations etc. This
            should be obtained by reading the raw data corresponding
            to the epochs with MNE.
        cachedir: str
            Path where the epochs are saved on disk. If this is
            None the epochs are returned.
        trialinfo_col: int
            Column in trialinfo which contains trial identifier.

    Output
        Returns path to saved epochs if cachedir is not None, else
        it returns the epochs
    '''
    if cachedir is None:
        return _load_ft_epochs(fname, rawinfo, trialinfo_col=trialinfo_col)
    epochs_path = os.path.join(cachedir, fname + '-epo.fif.gz')
    if not os.path.exists(epochs_path):
        epochs = _load_ft_epochs(fname, rawinfo, trialinfo_col=trialinfo_col)
        epochs.save(epochs_path)
    return epochs_path


def _load_ft_epochs(fname, rawinfo, trialinfo_col=-1):
    # load Matlab/Fieldtrip data
    f = h5py.File(fname)
    list(f.keys())
    ft_data = f['data']
    ft_data.keys()

    trialinfo = ft_data['trialinfo']
    channels = ft_data['label']
    sampleinfo = ft_data['sampleinfo']
    time = ft_data['time']
    sfreq = np.around(1 / np.diff(time[:].ravel()), 2)
    assert(len(np.unique(sfreq)) == 1)
    n_time, n_chans, n_trial = ft_data['trial'].shape

    data = np.zeros((n_trial, n_chans, n_time))
    transposed_data = np.transpose(ft_data['trial'])
    for trial in range(n_trial):
        data[trial, :, :] = transposed_data[trial]

    data = data[:, range(n_chans), :]

    chan_names = []
    for i in range(n_chans):
        st = channels[0][i]
        obj = ft_data[st]
        chan_names.append(''.join(chr(j) for j in obj[:]))
    #ch_names = [x + '-3705' for x in chan_names]

    info = create_info(chan_names, sfreq[0])
    events = np.zeros((n_trial, 3), int)
    events[:, 2] = trialinfo[trialinfo_col]
    events[:, 0] = sampleinfo[0]

    epochs = EpochsArray(data, info, events=events, tmin = min(time)[0], verbose=False)
    epochs.info = fix_chs(rawinfo, epochs.info)
    return epochs

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

ftdir = '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed4mne/'
megdir = '/home/pmurphy/meg_data/surprise/'
savedir = '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne_induced/'

for subj, tasks in subjects.items():
    for sess, rec in tasks:
        ftname = ftdir + str(subj) + '-' + str(sess) + '_0' + str(rec) + '_preproc4mne_induced.mat';
        print(ftname)
        megname = glob.glob(megdir + str(subj) + '-' + str(sess) + '*_0' + str(rec) + '.ds')
        print(megname[0])
        
        rawinfo = get_info_for_epochs(megname[0])
        epochs = _load_ft_epochs(ftname, rawinfo)
        epochs.save(savedir + str(subj) + '-' + str(sess) + '_0' + str(rec) + '_preproc4mne_induced.mat' + '-epo.fif.gz')
        del epochs
        



