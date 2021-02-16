'''
Preprocess an MEG data set.

The idea for preprocessing MEG data is modelled around a few aspects of the
confidence data set:
    1. Each MEG dataset is accompanied by a DataFrame that contains metadata for
       each trial.
    2. Trial meta data can be matched to the MEG data by appropriate triggers.
    3. Each MEG datafile contains several blocks of data that can be processed
       independently.

This leads to the following design:
    1. MEG is cut into recording blocks and artifact detection is carried out
    2. Each processed block is matched to meta data and timing (event on- and offsets)
       data is extracted from the MEG and aligned with the behavioral data.
    3. Data is epoched. Sync with meta data is guaranteed by a unique key for
       each trial that is stored along with the epoched data.
'''
from __future__ import division, print_function
import mne
import numpy as np
import pandas as pd
from pymeg.tools import hilbert
from pymeg import artifacts
import logging
from joblib import Memory
import os

memory = Memory(cachedir=os.environ['PYMEG_CACHE_DIR'], verbose=0)


def get_trial_periods(events, trial_start, trial_end):
    '''
    Parse trial start and end times from events.
    '''
    start = np.where(events[:, 2] == trial_start)[0]
    end = np.where(events[:, 2] == trial_end)[0]
    if not len(start) == len(end):
        start_times = events[start, 0]
        start = []
        end_times = events[end, 0]
        for i, e in enumerate(end_times):
            d = start_times - e
            d[d > 0] = -np.inf
            start_index = np.where(events[:, 0] == start_times[
                                   np.argmax(d)])[0][0]
            start.append(start_index)
    return np.array(start), end


def get_meta(raw, mapping, trial_pins, trial_start, trial_end, other_pins=None):
    '''
    Parse block structure from events in MEG files.

    Aggresively tries to fix introduced by recording crashes and late recording
    starts.

    mapping =
    '''

    def pins2num(pins):
        if len(pins) == 0:
            trial = 1
        else:
            # Convert pins to numbers
            trial = sum([2**(8 - pin) for pin in pins])
        return trial

    events, _ = get_events(raw)
    events = events.astype(float)

    if trial_start == trial_end:
        start = np.where(events[:, 2] == trial_start)[0]
        end = np.where(events[:, 2] == trial_end)[0]
        end = np.concatenate((end[1:] - 1, np.array([events.shape[0]])))

    else:
        start, end = get_trial_periods(events, trial_start, trial_end)

    trials = []
    for i, (ts, te) in enumerate(zip(start, end)):
        current_trial = {}
        trial_nums = events[ts:te + 1, 2]
        trial_times = events[ts:te + 1, 0]
        if trial_pins:
            # Find any pins that need special treatment, parse them and remove
            # triggers from trial_nums
            for key, value in trial_pins.items():
                if key in trial_nums:
                    pstart = np.where(trial_nums == key)[0][0] + 1
                    pend = pstart + np.where(trial_nums[pstart:] > 8)[0][0] + 1
                    pvals = trial_nums[pstart:pend]
                    current_trial[value] = pins2num(pvals)
                    trial_nums = np.concatenate(
                        (trial_nums[:pstart], trial_nums[pend:]))
                    trial_times = np.concatenate(
                        (trial_times[:pstart], trial_times[pend:]))

        for trigger, time in zip(trial_nums, trial_times):
            if trigger in mapping.keys():
                key = mapping[trigger][0]
                val = mapping[trigger][1]
            else:
                key = trigger
                val = time
            if key in current_trial.keys():
                try:
                    current_trial[key].append(current_trial[key][-1] + 1)
                    current_trial[key + '_time'].append(time)
                except AttributeError:
                    current_trial[str(key)] = [current_trial[
                        key], current_trial[key] + 1]
                    current_trial[
                        str(key) + '_time'] = [current_trial[str(key) + '_time'], time]
            else:
                current_trial[key] = val
                current_trial[str(key) + '_time'] = time
        trials.append(current_trial)

    meta = pd.DataFrame(trials)

    # Find other pins that are not trial related
    if other_pins:
        nums = events[:, 2]
        for key, value in other_pins.items():
            pstarts = np.where(nums == key)[0] + 1
            for pstart in pstarts:
                t = events[pstart, 0]
                pend = pstart + np.where(nums[pstart:] > 8)[0][0] + 1
                pvals = nums[pstart:pend]
                idx = meta.trial_start_time > t
                meta.loc[idx, value] = pins2num(pvals)

    time_fields = [c for c in meta if str(c).endswith('_time')]
    meta_fields = [c for c in meta if not str(c).endswith('_time')]
    return meta.loc[:, meta_fields], meta.loc[:, time_fields]


def preprocess_block(raw, blinks=True):
    '''
    Apply artifact detection to a block of data.
    '''
    ab = None
    artdef = {}
    if blinks:
        ab = artifacts.annotate_blinks(raw)
        artdef['blinks'] = ab
    am, zm = artifacts.annotate_muscle(raw)
    artdef['muscle'] = zm
    ac, zc, d = artifacts.annotate_cars(raw)
    artdef['cars'] = [zc, d]
    raw, ar, zj, jumps = artifacts.annotate_jumps(raw)
    artdef['jumps'] = zj
    ants = artifacts.combine_annotations(
        [x for x in [ab, am, ac, ar] if x is not None])
    #ants.onset += raw.first_samp/raw.info['sfreq']
    raw.annotations = ants
    artdef.update({'muscle': zm, 'cars': (zc, d), 'jumps': (zj, jumps)})
    return raw, ants, artdef


def mne_events(data, time_field, event_val):
    return np.vstack([
        data[time_field].values,
        0 * data[time_field].values,
        data[event_val].values]).astype(int).T


def ant2time_window(r, ant, onsets, epoch_time=(0, 1)):
    '''
    Create an annotation object that only contains events around time window.

    onsets are event onsets given in samples as defined in timing structure.
    '''
    onsets = (onsets - r.first_samp) / r.info['sfreq']
    event_starts = onsets + epoch_time[0]
    event_ends = onsets + epoch_time[1]
    new_onset, new_duration, new_description = [], [], []

    for ant_onset, ant_duration, description in zip(ant.onset, ant.duration,
                                                    ant.description):
        ant_end = ant_onset + ant_duration
        if all((ant_end < event_starts) | (event_ends < ant_onset)):
            pass
        else:
            new_onset.append(ant_onset)
            new_duration.append(ant_duration)
            new_description.append(description)
    return mne.annotations.Annotations(new_onset, new_duration, new_description)


def get_epoch(raw, meta, timing,
              event='stim_onset_t', epoch_time=(-.2, 1.5),
              epoch_label='hash', reject_time=(None, None)):
    '''
    Cut out epochs from raw data 

    Parameters
    ----------
    raw : raw data
    meta, timing : Dataframes that contain meta and timing information
    event : Column in timing that contains event onsets in sample time
    epoch_time : (start, end) in sec. relative to event onsets defined by 'event'    
    epoch_label : Column in meta that contains epoch labels.
    reject_time : time window for rejection.
    '''
    if reject_time[0] is None:
        reject_time = epoch_time
    print(reject_time)
    fields = set((event, epoch_label))
    joined_meta = (pd.concat([meta, timing], axis=1)
                   .loc[:, fields]
                   .dropna())

    ev = mne_events(joined_meta, event, epoch_label)

    annotations = raw.annotations
    new_ants = ant2time_window(
        raw, annotations, ev[:, 0], epoch_time=reject_time)
    print('Overlapping w/bad events:', len(new_ants.onset))
    raw.annotations = new_ants
    stim_period = mne.Epochs(raw, ev, tmin=epoch_time[0], tmax=epoch_time[1],
                             baseline=None,
                             reject_by_annotation=True,
                             reject_tmin=reject_time[0],
                             reject_tmax=reject_time[1])

    stim_period.load_data()
    if len(stim_period.events) == 0:
        raise RuntimeError('No trials left')
    stim_period = stim_period[[str(i) for i in stim_period.events[:, 2]]]
    raw.annotations = annotations
    return meta, stim_period


def concat(raws, metas, timings):
    '''
    Concatenate a set of raw objects and apply offset to meta to
    keep everything in sync. Should allow to load all sessions of
    a subject. Can then crop to parallelize.
    '''
    raws = [r.copy() for r in raws]
    offsets = np.cumsum([0] + [len(raw) for raw in raws])
    raw = raws[::-1].pop()
    raw.append(raws, preload=False)
    timings = [timing + offset for timing, offset in zip(timings, offsets)]
    timings = pd.concat(timings)
    metas = pd.concat(metas)
    return raw, metas, timings


def apply_baseline(epochs, baseline):
    '''
    Apply baseline correction to M/EEG channels
    '''
    drop_list = []
    chidx = np.array(
        [(x.startswith('M') or x.startswith('E'))
         for x in epochs.ch_names]).astype(bool)
    for epoch, orig in enumerate(epochs.selection):
        # Find baseline epoch for this.
        base = np.where(baseline.selection == orig)[0]
        if len(base) == 0:
            # Reject this one.
            drop_list.append(epoch)
        else:
            base_val = np.squeeze(baseline._data[base, chidx, :]).mean(1)
            epochs._data[epoch, chidx, :] -= base_val[:, np.newaxis]

    return epochs.drop(drop_list), drop_list


@memory.cache
def get_events_from_file(filename):
    raw = mne.io.read_raw_ctf(filename, system_clock='ignore')
    buttons = mne.find_events(raw, 'UPPT002', shortest_event=1)
    triggers = mne.find_events(raw, 'UPPT001', shortest_event=1)
    return triggers, buttons


def get_events(raw):
    buttons = mne.find_events(raw, 'UPPT002', shortest_event=1)
    triggers = mne.find_events(raw, 'UPPT001', shortest_event=1)
    return triggers, buttons


def load_epochs(filenames):
    epochs = []
    for f in filenames:
        e = mne.read_epochs(f)
        epochs.append(e)
    return epochs


def load_meta(filenames):
    return [pd.read_hdf(f, 'meta') for f in filenames]


def concatenate_epochs(epochs, metas):
    '''
    Concatenate a list of epoch and meta objects and set their dev_head_t projection to
    that of the first epoch.
    '''
    dev_head_t = epochs[0].info['dev_head_t']
    epoch_arrays = []
    processed_metas = []
    for e in ensure_iter(epochs):
        e.info['dev_head_t'] = dev_head_t
        e = mne.epochs.EpochsArray(e._data, e.info,
                                   events=e.events, tmin=e.tmin)
        epoch_arrays.append(e)

    if metas is not None:
        for m in ensure_iter(metas):
            processed_metas.append(m)

        return mne.concatenate_epochs(epoch_arrays), pd.concat(processed_metas)
    else:
        return mne.concatenate_epochs(epoch_arrays)


def combine_annotations(annotations, first_samples, last_samples, sfreq):
    '''
    Concatenate a list of annotations objects such that annotations
    stay in sync with the output of mne.concatenate_raws.

    This function assumes that annotations objects come from different raw objects
    that are to be concatenated. In this case the concatenated raw object retains
    the first sample of the first raw object and then treats the data as
    continuous. In contrast, the annotation onsets already shifted by each individual
    raw object's first sample to be in sync. When concatenting annotations this
    needs to be taken into account.

    Parameters
    ----------
    annotations : list of annotations objects, shape (n_objects,)
    first_samples : list of ints, shape (n_objects,)
        First sample of each annotations' raw object.
    last_samples : list of ints, shape (n_objects,)
        Last sample of each annotations' raw object.
    sfreq : int
        Sampling frequency of data in raw objects.
    '''
    if all([ann is None for ann in annotations]):
        return None
    durations = [(1 + l - f) / sfreq for f,
                 l in zip(first_samples, last_samples)]
    offsets = np.cumsum([0] + durations[:-1])

    onsets = [(ann.onset - (fs / sfreq)) + offset
              for ann, fs, offset in zip(annotations, first_samples, offsets) if ann is not None]
    if len(onsets) == 0:
        return mne.annotations.Annotations(onset=[], duration=[], description=None)
    onsets = np.concatenate(onsets) + (first_samples[0] / sfreq)
    return mne.annotations.Annotations(onset=onsets,
                                       duration=np.concatenate(
                                           [ann.duration for ann in annotations]),
                                       description=np.concatenate([ann.description for ann in annotations]))


def ensure_iter(input):
    if isinstance(input, str):
        yield input
    else:
        try:
            for item in input:
                yield item
        except TypeError:
            yield input
