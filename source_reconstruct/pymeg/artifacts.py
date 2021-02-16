'''
Artifact detection on MEG data.

TODO: Update doc.

For now
'''
import glob
import logging

import mne
import numpy as np
import pandas as pd

from .tools import hilbert


def annotate_blinks(raw, ch_mapping={'x': 'UADC002-3705', 'y': 'UADC003-3705', 'p': 'UADC004-3705'}):
    '''
    Detect blinks and annotate as bad blinks
    '''
    logging.info('Annotating blinks artifacts')
    x, y, p = eye_voltage2gaze(raw, ch_mapping=ch_mapping)
    xpos, ypos = x.ravel() / ppd(), y.ravel() / ppd()
    sc = saccade_detection(xpos, ypos, threshold=10, acc_thresh=2000, Hz=1200)
    blinks, scfilt = blink_detection(xpos, ypos, sc)
    if len(blinks) == 0:
        return None
    blink_onsets = raw.times[blinks[:, 0]]
    blink_durations = raw.times[blinks[:, 1]] - raw.times[blinks[:, 0]]
    return mne.Annotations(blink_onsets, blink_durations, 'bad blinks')


def annotate_muscle(raw, cutoff=10):
    logging.info('Annotating muscle artifacts')
    try:
        arts, z = detect_muscle(raw.copy(), cutoff=cutoff)
    except MemoryError:
        print('Memory Error detected:', raw)
        raise RuntimeError(
            'Memory error detected in annotate muscle: ' + raw.info['filename'])

    annotations = None
    if len(arts) > 0:
        annotations = mne.Annotations(arts[:, 0],
                                      arts[:, 1], 'bad muscle')
    return annotations, z


def annotate_cars(raw, cutoff=4.0, der_cutoff=7.):
    logging.info('Annotating car artifacts')
    arts, z, d = detect_cars(raw.copy(), cutoff=cutoff, der_cutoff=der_cutoff)
    annotations = None
    if len(arts) > 0:
        annotations = mne.Annotations(arts[:, 0], arts[:, 1], 'bad car')
    return annotations, z, d


def annotate_jumps(raw, cutoff=25, allowed_before_bad=np.inf):
    logging.info('Annotating jump artifacts')
    arts, z, jumps_per_channel = detect_jumps(raw.copy(), cutoff=cutoff)
    # Need to check for jumps_per_channel
    bads = [k for k, v in jumps_per_channel.items() if v >
            allowed_before_bad]
    arts = [arts[k] for k, v in jumps_per_channel.items() if v <=
            allowed_before_bad]

    if len(bads) > 0:
        if 'bads' in raw.info.keys():
            raw.info['bads'].extend(bads)
        else:
            raw.info['bads'] = bads
    a = []
    for k in arts:
        if len(k) is not 0:
            a.extend(k)
    arts = np.array(a)
    annotations = None
    try:
        if len(arts) > 0:
            annotations = mne.Annotations(arts[:, 0], arts[:, 1], 'bad jump')
    except IndexError:
        pass
    return raw, annotations, z, jumps_per_channel


def detect_cars(raw, cutoff=3.5, der_cutoff=5.0, frequency_band=(None, 1)):
    '''
    Detect cars artifacts on blocks and ignore intermediate data.

    This works analagously to the fieldtrip detect artifact routine.

    Data epochs that are already marked as bad by the annotations in raw will be
    excluded for computation of mean and std.
    '''

    logging.info('Detecting car events with cutoff %i' % cutoff)
    if not hasattr(raw, '_data'):
        logging.info('Loading data for car artifact detection')
        raw.load_data()
    raw.pick_channels([x for x in raw.ch_names if x.startswith('M')])
    raw.filter(l_freq=None, h_freq=1)
    hilb = hilbert(raw._data)
    del raw._data
    hilb = np.abs(hilb).astype(float)

    # Compute IGR and median
    Qs = np.percentile(hilb, [10, 50, 90], axis=1)
    IQR = Qs[2, :] - Qs[0, :]
    m = Qs[1, :]
    zh = ((((hilb - m[:, np.newaxis]) / IQR[:, np.newaxis]))**2).mean(0)

    # Normalize zh to have 80 between 0 an 1
    q80 = np.percentile(zh, [80])
    zh = zh / q80

    # Compute derivative of zh
    d = np.concatenate(([0], np.diff(zh)))
    # Normalize to have 80 between -1 and 1
    d = d / np.diff(np.percentile(d, [10, 80]))
    d[:int(raw.info['sfreq'])] = 0
    d[-int(raw.info['sfreq']):] = 0

    # Compute artifact borders
    art_borders = np.where(np.diff(np.concatenate([[0], zh > cutoff, [0]])))[0]
    artifacts = []
    for start, end in zip(art_borders[0::2], art_borders[1::2]):
        onset_t = max((start - 1) - int(2.5 * raw.info['sfreq']), 0)
        end_t = min((end) + int(2.5 * raw.info['sfreq']), len(zh))
        # Check for derivative
        if d[onset_t:end_t].min() < -der_cutoff and d[onset_t:end_t].max() > der_cutoff:
            onset_t /= raw.info['sfreq']
            end_t /= raw.info['sfreq']
            duration = end_t - onset_t
            artifacts.append((onset_t, duration))
    return np.array(artifacts), zh, d


def detect_muscle(raw, cutoff=10, frequency_band=(110, 140)):
    '''
    Detect muscle artifacts on blocks and ignore intermediate data.

    This works analagously to the fieldtrip detect artifact routine.

    Data epochs that are already marked as bad by the annotations in raw will be
    excluded for computation of mean and std.
    '''

    logging.info('Detecting muscle events with cutoff %i' % cutoff)
    if not hasattr(raw, '_data'):
        logging.info('Loading data for muscle artifact detection')
        raw.load_data()
    raw.pick_channels([x for x in raw.ch_names if x.startswith('M')])

    # filt = mne.filter.band_pass_filter(raw._data, raw.info['sfreq'],
    #                                   frequency_band[0], frequency_band[1],  method='iir',
    #                                   iir_params = dict(order=9, ftype='butter'),
    #                                   copy=False)
    filt = mne.filter.filter_data(raw._data, raw.info['sfreq'],
                                  l_freq=frequency_band[
                                      0], h_freq=frequency_band[1],  method='iir',
                                  iir_params=dict(order=9, ftype='butter'),
                                  copy=False)

    hilb = abs(hilbert(filt)).astype(float)
    del filt
    # Compute IGR and median
    Qs = np.percentile(hilb, [10, 50, 90], axis=1)
    IQR = Qs[2, :] - Qs[0, :]
    m = Qs[1, :]
    zh = ((((hilb - m[:, np.newaxis]) / IQR[:, np.newaxis]))**2).mean(0)

    # Normalize zh to have 80 between 0 an 1
    q80 = np.percentile(zh, [80])
    zh = zh / q80

    art_borders = np.where(np.diff(np.concatenate([[0], zh > cutoff, [0]])))[0]

    artifacts = []
    for start, end in zip(art_borders[0::2], art_borders[1::2]):
        artifacts.append((
                         ((start - 1) / raw.info['sfreq']) - 0.2,
                         ((end - start) / raw.info['sfreq']) + 0.2,
                         ))
    return np.array(artifacts), zh


def detect_jumps(raw, cutoff=25):
    '''
    Detect jumps by convolving with a jump detection filter.
    '''
    logging.info('Detecting muscle events with cutoff %i' % cutoff)
    if not hasattr(raw, '_data'):
        logging.info('Loading data for muscle artifact detection')
        raw.load_data()
    raw.pick_channels([x for x in raw.ch_names if x.startswith('M')])
    #filt = mne.filter.low_pass_filter(raw._data, raw.info['sfreq'], 1)
    filt = mne.filter.filter_data(
        raw._data, raw.info['sfreq'], l_freq=None, h_freq=1)

    jump_kernel = (np.array([1] * 50),
                   [0, 0],
                   np.array([-1] * 50))
    jump_kernel = np.concatenate(jump_kernel)
    filt = 0 * raw._data.copy()
    for i in range(filt.shape[0]):
        filt[i, :] = np.convolve(
            raw._data[i, :], jump_kernel, mode='same') - filt[i, :]
        filt[i, :int(len(jump_kernel) / 2)] = 0
        filt[i, -int(len(jump_kernel) / 2):] = 0

    # Compute IGR and median
    Qs = np.percentile(filt, [10, 50, 90], axis=1)
    IQR = Qs[2, :] - Qs[0, :]
    m = Qs[1, :]
    filt = (((filt - m[:, np.newaxis]) / IQR[:, np.newaxis]))**2
    # Need to keep information about channels here.

    artifacts = {}
    channel_count = {}
    for i, zh in enumerate(filt):
        artifacts[raw.ch_names[i]] = []
        art_borders = np.where(
            np.diff(np.concatenate([[0], zh > cutoff, [0]])))[0]
        channel_count[raw.ch_names[i]] = 0
        for start, end in zip(art_borders[0::2], art_borders[1::2]):
            artifacts[raw.ch_names[i]].append(
                ((start - 1) / raw.info['sfreq'], (end - start) / raw.info['sfreq']))
            channel_count[raw.ch_names[i]] += 1
    return artifacts, zh, channel_count


def eye_voltage2gaze(raw, ranges=(-5, 5), screen_x=(0, 1920),
                     screen_y=(0, 1080),
                     ch_mapping={'x': 'UADC002-3705', 'y': 'UADC003-3705', 'p': 'UADC004-3705'}):
    '''
    Convert analog output of EyeLink 1000+ to gaze coordinates.
    '''
    minvoltage, maxvoltage = ranges
    maxrange, minrange = 1., 0.
    screenright, screenleft = screen_x
    screenbottom, screentop = screen_y

    idx = np.where(np.array(raw.ch_names) == ch_mapping['x'])[0][0]
    R = (raw[idx, :][0] - minvoltage) / (maxvoltage - minvoltage)
    S = R * (maxrange - minrange) + minrange
    x = S * (screenright - screenleft + 1) + screenleft

    idy = np.where(np.array(raw.ch_names) == ch_mapping['y'])[0][0]
    R = (raw[idy, :][0] - minvoltage) / (maxvoltage - minvoltage)
    S = R * (maxrange - minrange) + minrange
    y = S * (screenbottom - screentop + 1) + screentop

    idp = np.where(np.array(raw.ch_names) == ch_mapping['p'])[0][0]
    p = raw[idp, :][0]
    return x, y, p


def eye_voltage2gaze_epochs(epochs, ranges=(-5, 5), screen_x=(0, 1920),
                            screen_y=(0, 1080),
                            ch_mapping={'x': 'UADC002-3705', 'y': 'UADC003-3705', 'p': 'UADC004-3705'}):
    '''
    Convert analog output of EyeLink 1000+ to gaze coordinates.
    '''
    minvoltage, maxvoltage = ranges
    maxrange, minrange = 1., 0.
    screenright, screenleft = screen_x
    screenbottom, screentop = screen_y

    idx = np.where(np.array(epochs.ch_names) == ch_mapping['x'])[0][0]
    R = (epochs._data[:, idx, :].squeeze() -
         minvoltage) / (maxvoltage - minvoltage)
    S = R * (maxrange - minrange) + minrange
    x = S * (screenright - screenleft + 1) + screenleft

    idy = np.where(np.array(epochs.ch_names) == ch_mapping['y'])[0][0]
    R = (epochs._data[:, idy, :].squeeze() -
         minvoltage) / (maxvoltage - minvoltage)
    S = R * (maxrange - minrange) + minrange
    y = S * (screenbottom - screentop + 1) + screentop

    idp = np.where(np.array(epochs.ch_names) == ch_mapping['p'])[0][0]
    p = epochs._data[:, idp:].squeeze()
    return x, y, p

velocity_window_size = 3


def get_velocity(x, y, Hz):
    '''
    Compute velocity of eye-movements.

    'x' and 'y' specify the x,y coordinates of gaze location. The function
    assumes that the values in x,y are sampled continously at a rate specified
    by 'Hz'.
    '''
    Hz = float(Hz)
    distance = ((np.diff(x) ** 2) +
                (np.diff(y) ** 2)) ** .5
    distance = np.hstack(([distance[0]], distance))
    win = np.ones((velocity_window_size)) / float(velocity_window_size)
    velocity = np.convolve(distance, win, mode='same')
    velocity = velocity / (velocity_window_size / Hz)
    acceleration = np.diff(velocity) / (1. / Hz)
    acceleration = np.abs(np.hstack(([acceleration[0]], acceleration)))
    return velocity, acceleration


def saccade_detection(x, y, Hz=1200, threshold=30,
                      acc_thresh=2000):
    '''
    Detect saccades in a stream of gaze location samples.

    Coordinates of x,y are assumed to be in degrees.

    Saccades are detect by a velocity/acceleration threshold approach.
    A saccade starts when a) the velocity is above threshold, b) the
    acceleration is above acc_thresh at least once during the interval
    defined by the velocity threshold.
    '''

    velocity, acceleration = get_velocity(x, y, float(Hz))
    saccades = (velocity > threshold)

    borders = np.where(np.diff(saccades.astype(int)))[0] + 1
    if velocity[1] > threshold:
        borders = np.hstack(([0], borders))

    saccade = 0 * np.ones(x.shape)

    saccade_times = []
    # Only count saccades when acceleration also surpasses threshold
    for i, (start, end) in enumerate(zip(borders[0::2], borders[1::2])):
        if np.sum(acceleration[start:end] > acc_thresh) >= 1:
            saccade[start:end] = 1
            saccade_times.append((start, end))

    return np.array(saccade_times)


def microssacade_detection(x, y, VFAC):
    '''
    Microsaccade detection a la Engbert et al.
    '''
    if len(x) < 5:
        return None
    dt = 1 / 1200.
    kernel = np.array([1., 1., 0., -1., -1.])
    vx = np.convolve(x, kernel, mode='same') / (6 * dt)
    vy = np.convolve(y, kernel, mode='same') / (6 * dt)
    msdx = np.sqrt(np.median((vx - np.median(vx))**2))
    msdy = np.sqrt(np.median((vy - np.median(vy))**2))
    radiusx = VFAC * msdx
    radiusy = VFAC * msdy
    test = (vx / radiusx)**2 + (vy / radiusy)**2
    borders = np.where(np.diff((test > 1).astype(int)))[0] + 1
    if test[0] > 1:
        borders = np.hstack(([0], borders))
    if test[-1] > 1:
        borders = np.hstack((borders, [len(x)]))

    borders = borders.reshape(len(borders) / 2, 2)
    return borders


def blink_detection(x, y, saccades):
    '''
    A blink is everything that is surrounded by two saccades and period in
    between where the eye is off screen.
    '''
    rm_sac = (saccades[:, 0] * 0).astype(bool)
    blinks = []
    skipnext = False
    for i, ((pss, pse), (nss, nse)) in enumerate(zip(saccades[:-1], saccades[1:])):
        if skipnext:
            skipnext = False
            continue
        xavg = x[pse:nss].mean()
        yavg = y[pse:nss].mean()

        if (xavg > 40) and (yavg > 20):
            rm_sac[i:i + 2] = True
            blinks.append((pss, nse))
            skip_next = True

    return np.array(blinks), saccades[~rm_sac, :]


def nan_bad_epochs(data, raw):
    '''
    Overwrite data with NANs for all epochs in the Nx3 artifacts matrix.
    '''
    Hz = raw.info['sfreq']
    if raw.annotations is not None:
        for start, duration, desc in zip(raw.annotations.onset,
                                         raw.annotations.duration,
                                         raw.annotations.description):
            if not 'bad' in desc:
                continue
            start = int(start * Hz)
            data[start:start + int(duration * Hz)] = nan
    return data


def ppd(vieweing_distance=62.5, screen_width=38.0, x_resolution=1450):
    '''
    Compute pixels per degree for the current setup.
    '''
    o = np.tan(0.5 * np.pi / 180) * vieweing_distance
    return 2 * o * x_resolution / screen_width


def combine_annotations(annotations):
    '''
    Add annotations to a raw object. Makes sure that old annotations are kept.
    '''
    if len(annotations) == 0:
        return mne.Annotations(np.array([0]), np.array([0]), 'dummy')
    elif len(annotations) == 1:
        return annotations[0]
    else:
        old = annotations[0]
        for new in annotations[1:]:
            if new is None:
                continue
            orig_time = None
            onsets = np.concatenate((old.onset, new.onset))
            duration = np.concatenate((old.duration, new.duration))
            descr = np.concatenate((old.description, new.description))
            if old.orig_time is not None:
                orig_time = np.concatenate((old.orig_time, new.orig_time))
            old = mne.Annotations(onsets, duration, descr, orig_time)
        return old


def concatenate_annotations(annotations, durations):
    '''
    Concatenat annotations for several consecutive raw objects that are to be
    merged. Each annotation must come with one duration for the raw object that
    it belongs to.

    durations is in s.

    TODO: This was because of a bug in pymne - which should have been fixed upstream by now.
    '''
    durations[0] = 0
    onsets = np.concatenate([a.onset + (li)
                             for a, li in zip(annotations, cumsum(durations))])
    durations = np.concatenate([a.duration
                                for a in annotations])
    description = np.concatenate([a.description
                                  for a in annotations])
    return mne.Annotations(onsets, durations, description)
