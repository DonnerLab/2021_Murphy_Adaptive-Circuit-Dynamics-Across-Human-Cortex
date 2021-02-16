import os
from os.path import join
import mne
import pandas as pd
import numpy as np
from joblib import Memory

from . import preprocessing


memory = Memory(cachedir=os.environ['PYMEG_CACHE_DIR'], verbose=0)
subjects_dir = '/home/pmurphy/meg_data/surprise/MRIs/fs_converted/'  # set to SUBJECTS_DIR for this study


def set_fs_subjects_dir(directory):
    """Set freesurfer subjectdir environment variable"""
    global subjects_dir
    os.environ['SUBJECTS_DIR'] = directory
    subjects_dir = directory


def check_bems(subjects):
    """Create a plot of all BEM segmentations."""
    for sub in subjects:
        fig = mne.viz.plot_bem(subject=sub,
                               subjects_dir=subjects_dir,
                               brain_surfaces='white',
                               orientation='coronal')


@memory.cache
def get_source_space(subject):
    """Return source space.

    Mainly a helper function to provide caching of source space
    computation.
    """
    return mne.setup_source_space(subject, spacing='oct6',   # default spacing is 'oct6', but this may mean 
                                  subjects_dir=subjects_dir, # FEF is not found - in such cases, switch to 'oct7'
                                  add_dist=False)


@memory.cache
def get_info(raw_filename, epochs_filename):
    """Return an info dict for a set of epochs.

    Args:
        raw_filename : string
            Path of raw data that was the basis for creating
            the epochs.
        epochs_filename: string
            Epochs file name.

    Returns:
        MNE Infor structure.
    """
    trans, fiducials, info = get_head_correct_info(
        raw_filename, epochs_filename)
    return info


@memory.cache
def get_leadfield(subject, raw_filename, epochs_filename, trans_filename,
                  conductivity=(0.3, 0.006, 0.3), njobs=4, bem_sub_path='bem'):
    """Compute leadfield with presets for this subject

    Args:    
        subject : str
            Name of freesurfer subject
        raw_filename : str
            Filename that points to the raw data for this lead field.
            This file will be used to extract a CTF transformation matrix
            and an info struct.
        epochs_filename : str
            Filename from which fiducial locations will be extracted.
        trans_filename : str
            Points to transformation file between fiducials and MRI.
        conductivity : 3-tuple of floats
            Conductivities for BEM model
        njobs: int
            Number of cores to paralellize over
        bem_sub_path: str
            Sub-path of freesurfer subject path where to read bem
            surfaces from

    Returns:   
        Tuple of (forward model, BEM model, source space)
    """
    src = get_source_space(subject)
    model = make_bem_model(
        subject=subject,
        ico=None,
        conductivity=conductivity,
        subjects_dir=subjects_dir,
        bem_sub_path=bem_sub_path)
    bem = mne.make_bem_solution(model)
    info = get_info(raw_filename, epochs_filename)
    fwd = mne.make_forward_solution(
        info,
        trans=trans_filename,
        src=src,
        bem=bem,
        meg=True,
        eeg=False,
        mindist=2.5,
        n_jobs=njobs)
    return fwd, bem, fwd['src']


def make_bem_model(subject, ico=4, conductivity=(0.3, 0.006, 0.3),
                   subjects_dir=None, verbose=None, bem_sub_path='bem'):
    """Create a BEM model for a subject.

    Copied from MNE python, adapted to read surface from fieldtrip / spm
    segmentation.
    """
    import os.path as op
    from mne.io.constants import FIFF
    conductivity = np.array(conductivity, float)
    if conductivity.ndim != 1 or conductivity.size not in (1, 3):
        raise ValueError('conductivity must be 1D array-like with 1 or 3 '
                         'elements')
    subjects_dir = mne.utils.get_subjects_dir(subjects_dir, raise_error=True)
    subject_dir = op.join(subjects_dir, subject)
    bem_dir = op.join(subject_dir, bem_sub_path)
    inner_skull = op.join(bem_dir, 'inner_skull.surf')
    outer_skull = op.join(bem_dir, 'outer_skull.surf')
    outer_skin = op.join(bem_dir, 'outer_skin.surf')
    surfaces = [inner_skull, outer_skull, outer_skin]
    ids = [FIFF.FIFFV_BEM_SURF_ID_BRAIN,
           FIFF.FIFFV_BEM_SURF_ID_SKULL,
           FIFF.FIFFV_BEM_SURF_ID_HEAD]
    if len(conductivity) == 1:  # use only 1-layer BEM model in desired cases
        surfaces = surfaces[:1]
        ids = ids[:1]
    surfaces = mne.bem._surfaces_to_bem(surfaces, ids, conductivity, ico)
    mne.bem._check_bem_size(surfaces)
    return surfaces


@memory.cache
def get_labels(subject, filters=['*wang2015atlas*', '*JWDG.lr*'],
               annotations=['HCPMMP1'], sdir=None):
    """Read ROI labels from annotations and label files.

    This defines the ROIs that can be used for source reconstruction.
    ROIs originate from either freesurfer annotation files or label 
    files. 

    Args:
        subject: str
            Name of freesurfer subject
        filters: list of glob strings.
            A list of strings that select label files by globbing them,
            e.g. shell based selection of strings ('*wang*' selects all
            labels that contain the string 'wang').
        annotations: list of str
            Name of annotation files to load
        sdir: str, default None
            Overwrites freesurfer subject dir.

    Returns:
        List of MNE label objects.
    """
    global subjects_dir
    import glob
    if sdir is not None:
        subject_dir = sdir
    else:
        subject_dir = subjects_dir

    labels = []
    for filter in filters:
        labels += glob.glob(join(subject_dir, subject, 'label', filter))
    labels = [mne.read_label(label, subject) for label in labels]
    for annotation in annotations:
        annot = mne.read_labels_from_annot(
            subject, parc=annotation, subjects_dir=subject_dir)
        annot = [a for a in annot if not '???' in a.name]
        labels.extend(annot)
    return labels

def labels_exclude(labels, exclude_filters=['wang2015atlas.IPS4', 'wang2015atlas.IPS5', 
                                        'wang2015atlas.SPL', 'JWG_lat_Unknown']):
    labels_to_exclude = []
    for l in labels:
        is_exclude = np.array([(f in l.name) for f in exclude_filters]).sum() > 0
        if is_exclude:
            labels_to_exclude.append(l)
    for l in labels_to_exclude:
        labels.remove(l)
    return labels

def labels_remove_overlap(labels, priority_filters=['wang', 'JWG'],):
    
    
    vertices_nr_ori = sum([l.vertices.shape[0] for l in labels])

    labels_no_overlap = []

    for hemi in ['lh', 'rh']:

        labels_hemi = [l for l in labels if l.hemi == hemi]

        # get all category 1 vertices:
        cat1_vertices = []
        for l in labels_hemi:
            is_priority = np.array([(f in l.name) for f in priority_filters]).sum() > 0
            if is_priority:
                cat1_vertices.append(l.vertices)
        cat1_vertices = np.unique(np.concatenate(cat1_vertices))

        # remove category1 vertices from all other labels:
        for i, l in enumerate(labels_hemi):
            is_priority = np.array([(f in l.name) for f in priority_filters]).sum() > 0
            if not is_priority:
                cat1_indices = np.isin(labels_hemi[i].vertices, cat1_vertices)
                labels_hemi[i].vertices = labels_hemi[i].vertices[~cat1_indices]
        
        labels_no_overlap.extend(labels_hemi)

    vertices_nr_new = sum([l.vertices.shape[0] for l in labels_no_overlap])
    
    print('vertices ori: {}'.format(vertices_nr_ori))
    print('vertices new:  {}'.format(vertices_nr_new))
    print('excluded: {}'.format(vertices_nr_ori-vertices_nr_new))
    
    return labels_no_overlap
    
    # vertices_nr_ori = sum([l.vertices.shape[0] for l in labels])

    ## get all category 1 vertices:
    # cat1_vertices = []
    # for l in labels:
    #     is_priority = np.array([(f in l.name) for f in priority_filters]).sum() > 0
    #     if is_priority:
    #         cat1_vertices.append(l.vertices)
    # cat1_vertices = np.unique(np.concatenate(cat1_vertices))

    ## remove category1 vertices from all other labels:
    # for i, l in enumerate(labels):
    #     is_priority = np.array([(f in l.name) for f in priority_filters]).sum() > 0
    #     if not is_priority:
    #         cat1_indices = np.isin(labels[i].vertices, cat1_vertices)
    #         labels[i].vertices = labels[i].vertices[~cat1_indices]
    
    # vertices_nr_new = sum([l.vertices.shape[0] for l in labels])
    
    # print('vertices ori: {}'.format(vertices_nr_ori))
    # print('vertices new:  {}'.format(vertices_nr_new))
    # print('excluded: {}'.format(vertices_nr_ori-vertices_nr_new))
    
    # return labels


'''
Transformation matrix MEG<>T1 space.
'''


@memory.cache
def get_head_correct_info(raw_filename, epoch_filename, N=-1):
    """Get transformation matrix, fiducial positions and infor structure.

    The returned info structure contains fiducial locations computed from 
    the epoch data. 
    """
    trans = get_ctf_trans(raw_filename)
    fiducials = get_ref_head_pos(epoch_filename, trans, N=N)
    raw = mne.io.ctf.read_raw_ctf(raw_filename)
    info = replace_fiducials(raw.info, fiducials)
    return trans, fiducials, info



def get_trans_epoch(raw_filename, epoch_filename):
    from os.path import join
    save_path=os.environ['PYMEG_CACHE_DIR']
    save_file = join(save_path, 
        epoch_filename.split("/")[-1].replace(".fif","").replace(".gz", "") + "-trans.fif")
    if os.path.isfile(save_file):
        return save_file
    trans, fiducials, info = get_head_correct_info(
        raw_filename, epoch_filename)

    mne.io.meas_info.write_info(save_file, info)
    return save_file


def make_trans(subject, raw_filename, epoch_filename, trans_name):
    """Create coregistration between MRI and MEG space.

    Call MNE gui to create a MEG<>MRI transformation matrix
    """
    import os
    import time
    fid_epochs = get_trans_epoch(raw_filename, epoch_filename)
    cmd = 'mne coreg --high-res-head -d %s -s %s -f %s' % (
        subjects_dir, subject, fid_epochs)
    print(cmd)
    os.system(cmd)
    mne.gui.coregistration(subject, inst=hs_ref.name,
                           subjects_dir=subjects_dir)
    print('--------------------------------')
    print('Please save trans file as:')
    print(trans_name)
    while not os.path.isfile(trans_name):
        #print('Waiting for transformation matrix to appear')
        time.sleep(5)



@memory.cache
def get_ref_head_pos(filename,  trans, N=-1):
    """Compute average head position from epochs.

    Args:
        filename: str
            Epochs file to load
        trans: dict
            A dictionary that contains t_ctf_dev_dev
            transformation matrix, e.g. output of
            get_ctf_trans
    Returns:
        Dictionary that contains average fiducial positions.
    """
    from mne.transforms import apply_trans
    data = preprocessing.load_epochs([filename])[0]
    cc = head_loc(data.decimate(10))
    nasion = np.stack([c[0] for c in cc[:N]]).mean(0)
    lpa = np.stack([c[1] for c in cc[:N]]).mean(0)
    rpa = np.stack([c[2] for c in cc[:N]]).mean(0)
    nasion, lpa, rpa = nasion.mean(-1), lpa.mean(-1), rpa.mean(-1)

    return {'nasion': apply_trans(trans['t_ctf_dev_dev'], np.array(nasion)),
            'lpa': apply_trans(trans['t_ctf_dev_dev'], np.array(lpa)),
            'rpa': apply_trans(trans['t_ctf_dev_dev'], np.array(rpa))}


def replace_fiducials(info, fiducials):
    """Replace initial fiducial measuremnt with new estimates

    CTF systems measure fiducial location at the beginning of the measurement.
    When used with online head loc over multiple sessions these measurements
    are not accurate. This is because subjects are guided to the head position
    of previous sessions.

    Args:
        info: MNE info structure
        fiducials: dict
            Dictionary that contains fiducial positions, e.g.
            see output of get_ref_head_pos.
    Returns:
        Info structure with updated head position.
    """
    from mne.io import meas_info
    fids = meas_info._make_dig_points(**fiducials)
    info = info.copy()
    dig = info['dig']
    for i, d in enumerate(dig):
        if d['kind'] == 3:
            if d['ident'] == 3:
                dig[i]['r'] = fids[2]['r']
            elif d['ident'] == 2:
                dig[i]['r'] = fids[1]['r']
            elif d['ident'] == 1:
                dig[i]['r'] = fids[0]['r']
    info['dig'] = dig
    return info


def head_movement(epochs):
    """Compute head movement from epochs.

    Returns the circumcenter of the three fiducials for each time point.
    """
    ch_names = np.array(epochs.ch_names)
    channels = {'x': ['HLC0011', 'HLC0012', 'HLC0013'],
                'y': ['HLC0021', 'HLC0022', 'HLC0023'],
                'z': ['HLC0031', 'HLC0032', 'HLC0033']}
    channel_ids = {}
    for key, names in channels.items():
        ids = [np.where([n in ch for ch in ch_names])[0][0] for n in names]
        channel_ids[key] = ids

    data = epochs._data
    ccs = []
    for e in range(epochs._data.shape[0]):
        x = np.stack([data[e, i, :] for i in channel_ids['x']])
        y = np.stack([data[e, i, :] for i in channel_ids['y']])
        z = np.stack([data[e, i, :] for i in channel_ids['z']])
        cc = circumcenter(x, y, z)
        ccs.append(cc)
    return np.stack(ccs)


@memory.cache
def get_head_loc(epochs):
    cc = head_loc(epochs)
    trans, fiducials, info = get_head_correct_info(subject, session)
    nose_coil = np.concatenate([c[0] for c in cc], -1)
    left_coil = np.concatenate([c[1] for c in cc], -1)
    right_coil = np.concatenate([c[2] for c in cc], -1)
    nose_coil = apply_trans(trans['t_ctf_dev_dev'], nose_coil.T)
    left_coil = apply_trans(trans['t_ctf_dev_dev'], left_coil.T)
    right_coil = apply_trans(trans['t_ctf_dev_dev'], right_coil.T)

    nose_coil = (nose_coil**2).sum(1)**.5
    left_coil = (left_coil**2).sum(1)**.5
    right_coil = (right_coil**2).sum(1)**.5
    return nose_coil, left_coil, right_coil


def head_loc(epochs):
    ch_names = np.array(epochs.ch_names)
    channels = {'x': ['HLC0011', 'HLC0012', 'HLC0013'],
                'y': ['HLC0021', 'HLC0022', 'HLC0023'],
                'z': ['HLC0031', 'HLC0032', 'HLC0033']}
    channel_ids = {}
    for key, names in channels.items():
        ids = [np.where([n in ch for ch in ch_names])[0][0] for n in names]
        channel_ids[key] = ids

    data = epochs._data
    ccs = []
    if len(epochs._data.shape) > 2:
        for e in range(epochs._data.shape[0]):
            x = np.stack([data[e, i, :] for i in channel_ids['x']])
            y = np.stack([data[e, i, :] for i in channel_ids['y']])
            z = np.stack([data[e, i, :] for i in channel_ids['z']])
            ccs.append((x, y, z))
    else:
        x = np.stack([data[i, :] for i in channel_ids['x']])
        y = np.stack([data[i, :] for i in channel_ids['y']])
        z = np.stack([data[i, :] for i in channel_ids['z']])
        ccs.append((x, y, z))
    return ccs


def get_ctf_trans(filename):
    """Get transformation matrix between sensors and head space."""
    from mne.io.ctf.res4 import _read_res4
    from mne.io.ctf.hc import _read_hc
    from mne.io.ctf.trans import _make_ctf_coord_trans_set

    res4 = _read_res4(filename)  # Read the magical res4 file
    coils = _read_hc(filename)  # Read the coil locations

    # Investigate the coil location data to get the coordinate trans
    coord_trans = _make_ctf_coord_trans_set(res4, coils)
    return coord_trans


def circumcenter(coil1, coil2, coil3):
    """Determines position and orientation of the circumcenter of fiducials.
    Adapted from:    
    http://www.fieldtriptoolbox.org/example/how_to_incorporate_head_movements_in_meg_analysis
    CIRCUMCENTER determines the position and orientation of the circumcenter
    of the three fiducial markers (MEG headposition coils).

    Args:
        coil1-3: 3xN array
            X,y,z-coordinates of the 3 coils [3 X N],[3 X N],[3 X N] where N
            is timesamples/trials.
    Returns:
        X,y,z-coordinates of the circumcenter [1-3 X N], and the orientations
        to the x,y,z-axes [4-6 X N].
    A. Stolk, 2012
    """
    N = coil1.shape[1]
    cc = np.zeros((6, N)) * np.nan
    # x-, y-, and z-coordinates of the circumcenter
    # use coordinates relative to point `a' of the triangle
    xba = coil2[0, :] - coil1[0, :]
    yba = coil2[1, :] - coil1[1, :]
    zba = coil2[2, :] - coil1[2, :]
    xca = coil3[0, :] - coil1[0, :]
    yca = coil3[1, :] - coil1[1, :]
    zca = coil3[2, :] - coil1[2, :]

    # squares of lengths of the edges incident to `a'
    balength = xba * xba + yba * yba + zba * zba
    calength = xca * xca + yca * yca + zca * zca

    # cross product of these edges
    xcrossbc = yba * zca - yca * zba
    ycrossbc = zba * xca - zca * xba
    zcrossbc = xba * yca - xca * yba

    # calculate the denominator of the formulae
    denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc
                         + zcrossbc * zcrossbc)

    # calculate offset (from `a') of circumcenter
    xcirca = ((balength * yca - calength * yba) * zcrossbc -
              (balength * zca - calength * zba) * ycrossbc) * denominator
    ycirca = ((balength * zca - calength * zba) * xcrossbc -
              (balength * xca - calength * xba) * zcrossbc) * denominator
    zcirca = ((balength * xca - calength * xba) * ycrossbc -
              (balength * yca - calength * yba) * xcrossbc) * denominator

    cc[0, :] = xcirca + coil1[0, :]
    cc[1, :] = ycirca + coil1[1, :]
    cc[2, :] = zcirca + coil1[2, :]
    # orientation of the circumcenter with respect to the x-, y-, and z-axis
    # coordinates
    v = np.stack([cc[0, :].T, cc[1, :].T, cc[2, :].T]).T
    vx = np.stack([np.zeros((N,)).T, cc[1, :].T, cc[2, :].T]).T
    # on the x - axis
    vy = np.stack([cc[0, :].T, np.zeros((N,)).T, cc[2, :].T]).T
    # on the y - axis
    vz = np.stack([cc[0, :].T, cc[1, :].T, np.zeros((N,)).T]).T
    # on the z - axis
    thetax, thetay = np.zeros((N,)) * np.nan, np.zeros((N,)) * np.nan
    thetaz = np.zeros((N,)) * np.nan
    for j in range(N):

        # find the angles of two vectors opposing the axes
        thetax[j] = np.arccos(np.dot(v[j, :], vx[j, :]) /
                              (np.linalg.norm(v[j, :]) * np.linalg.norm(vx[j, :])))
        thetay[j] = np.arccos(np.dot(v[j, :], vy[j, :]) /
                              (np.linalg.norm(v[j, :]) * np.linalg.norm(vy[j, :])))
        thetaz[j] = np.arccos(np.dot(v[j, :], vz[j, :]) /
                              (np.linalg.norm(v[j, :]) * np.linalg.norm(vz[j, :])))

        # convert to degrees
        cc[3, j] = (thetax[j] * (180 / np.pi))
        cc[4, j] = (thetay[j] * (180 / np.pi))
        cc[5, j] = (thetaz[j] * (180 / np.pi))
    return cc


def ensure_iter(input):
    if isinstance(input, basestring):
        yield input
    else:
        try:
            for item in input:
                yield item
        except TypeError:
            yield input


def clear_cache():
    memory.clear()


def add_volume_info(subject, surface, subjects_dir, volume='T1'):
    """Add volume info from MGZ volume
    """
    import os.path as op
    from mne.bem import _extract_volume_info
    from mne.surface import (read_surface, write_surface)
    subject_dir = op.join(subjects_dir, subject)
    mri_dir = op.join(subject_dir, 'mri')
    T1_mgz = op.join(mri_dir, volume + '.mgz')
    new_info = _extract_volume_info(T1_mgz)
    print(new_info.keys())
    rr, tris, volume_info = read_surface(surface,
                                         read_metadata=True)

    # volume_info.update(new_info)  # replace volume info, 'head' stays
    print(volume_info.keys())
    import numpy as np
    if 'head' not in volume_info.keys():
        volume_info['head'] = np.array([2,  0, 20], dtype=np.int32)
    write_surface(surface, rr, tris, volume_info=volume_info)

