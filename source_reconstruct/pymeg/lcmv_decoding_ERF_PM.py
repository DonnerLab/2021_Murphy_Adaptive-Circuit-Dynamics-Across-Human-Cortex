# Some imports:
import logging
import mne
import numpy as np
import os
import scipy.io as sio

from joblib import Memory # Provides caching of results

from os import makedirs
from os.path import join
from glob import glob

from pymeg import lcmv as pymeglcmv
from pymeg import source_reconstruction as pymegsr
from pymeg import decoding_ERF

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn import linear_model

import datetime
import pandas as pd

from pymeg.lcmv_peter import get_stim_epoch
from pymeg.source_reconstruction import (
    get_leadfield,
    make_trans,
)

# Setup some paths:
memory = Memory(cachedir='/mnt/homes/home024/pmurphy/tmp/')
subjects_dir = "/home/pmurphy/meg_data/surprise/MRIs/fs_converted" # freesurfer subject dirs
trans_dir = "/home/pmurphy/meg_data/surprise/MRIs/trans_mats" # transofrmation matrices

# We need to restrict the number of threads that we use for the cluster
def set_n_threads(n):
    os.environ["OPENBLAS_NUM_THREADS"] = str(n)
    os.environ["MKL_NUM_THREADS"] = str(n)
    os.environ["OMP_NUM_THREADS"] = str(n)

# make dict of subjects/sessions/recordings
# subjects = {'DCB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'DHB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1)],
#             'ECB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'EMB': [(1, 1), (1, 2), (2, 1), (2, 3), (3, 1), (3, 2)],
#             'EXF': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'EXG': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'GSB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'HBC': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'JTB': [(1, 1), (1, 2), (2, 1), (2, 3), (3, 1), (3, 2)],
#             'KSV': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'NIF': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'OMF': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'PDP': [(1, 1), (1, 2), (2, 2), (2, 3), (3, 1), (3, 2)],
#             'QNV': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2)],
#             'TFD': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
#             'TNB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)],
#             'TSJ': [(1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2)]}
subjects = {'DCB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'DHB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1)],
            'ECB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'EMB': [(1, 1), (1, 2), (2, 1), (2, 3), (3, 1), (3, 2)],
            'EXF': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'EXG': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'GSB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)]}

# typical processing demands (in # cores) per subject
mem_demand = {'DCB': 4, 'DHB': 3, 'ECB': 3, 'EMB': 3, 'EXF': 3, 'EXG': 3,
            'GSB': 3, 'HBC': 4, 'JTB': 3, 'KSV': 4, 'NIF': 3, 'OMF': 3,
            'PDP': 3, 'QNV': 4, 'TFD': 3, 'TNB': 3, 'TSJ': 3}

# typical processing demands (in # vertices) per areas
mem_area = {'vfcPrimary': 3650, 'vfcEarly': 9610, 'vfcV3ab': 3280,
            'vfcIPS01': 4200, 'vfcIPS23': 2200, 'JWG_aIPS': 870,
            'JWG_IPS_PCeS': 4600, 'JWG_M1': 2900, 'HCPMMP1_premotor': 9900}

# Mapping areas to labels
areas_to_labels = {
    "vfcPrimary": [
        u"lh.wang2015atlas.V1d-lh", u"rh.wang2015atlas.V1d-rh",
        u"lh.wang2015atlas.V1v-lh", u"rh.wang2015atlas.V1v-rh",
    ],
    "vfcEarly": [
        u"lh.wang2015atlas.V2d-lh", u"rh.wang2015atlas.V2d-rh",
        u"lh.wang2015atlas.V2v-lh", u"rh.wang2015atlas.V2v-rh",
        u"lh.wang2015atlas.V3d-lh", u"rh.wang2015atlas.V3d-rh",
        u"lh.wang2015atlas.V3v-lh", u"rh.wang2015atlas.V3v-rh",
        u"lh.wang2015atlas.hV4-lh", u"rh.wang2015atlas.hV4-rh",
    ],
#     "vfcVO": [
#         u"lh.wang2015atlas.VO1-lh", u"rh.wang2015atlas.VO1-rh",
#         u"lh.wang2015atlas.VO2-lh", u"rh.wang2015atlas.VO2-rh",
#     ],
#     "vfcPHC": [
#         u"lh.wang2015atlas.PHC1-lh", u"rh.wang2015atlas.PHC1-rh",
#         u"lh.wang2015atlas.PHC2-lh", u"rh.wang2015atlas.PHC2-rh",
#     ],
    "vfcV3ab": [
        u"lh.wang2015atlas.V3A-lh", u"rh.wang2015atlas.V3A-rh",
        u"lh.wang2015atlas.V3B-lh", u"rh.wang2015atlas.V3B-rh",
    ],
#     "vfcTO": [
#         u"lh.wang2015atlas.TO1-lh", u"rh.wang2015atlas.TO1-rh",
#         u"lh.wang2015atlas.TO2-lh", u"rh.wang2015atlas.TO2-rh",
#     ],
#     "vfcLO": [
#         u"lh.wang2015atlas.LO1-lh", u"rh.wang2015atlas.LO1-rh",
#         u"lh.wang2015atlas.LO2-lh", u"rh.wang2015atlas.LO2-rh",
#     ],
    "vfcIPS01": [
        u"lh.wang2015atlas.IPS0-lh", u"rh.wang2015atlas.IPS0-rh",
        u"lh.wang2015atlas.IPS1-lh", u"rh.wang2015atlas.IPS1-rh",
    ],
    "vfcIPS23": [
        u"lh.wang2015atlas.IPS2-lh", u"rh.wang2015atlas.IPS2-rh",
        u"lh.wang2015atlas.IPS3-lh", u"rh.wang2015atlas.IPS3-rh",
    ],
    "JWG_aIPS": ["lh.JWDG.lr_aIPS1-lh", "rh.JWDG.lr_aIPS1-rh", ],
    "JWG_IPS_PCeS": ["lh.JWDG.lr_IPS_PCes-lh", "rh.JWDG.lr_IPS_PCes-rh", ],
    "JWG_M1": ["lh.JWDG.lr_M1-lh", "rh.JWDG.lr_M1-rh", ],
    "HCPMMP1_premotor": (
        ["L_{}_ROI-lh".format(area) for area in [
            "55b", "6d", "6a", "FEF", "6v", "6r", "PEF"]] +
        ["R_{}_ROI-rh".format(area) for area in [
            "55b", "6d", "6a", "FEF", "6v", "6r", "PEF"]]
    ),
}


# Submit to cluster. One job per ROI and subject
def submit(only_glasser=False):
    from pymeg import parallel
    from itertools import product

    for area in list(areas_to_labels.keys()):
        for subject in subjects.keys():
            print("Submitting %s -> %s" % (subject, area))
            parallel.pmap(
                decode,
                [(subject, area, subjects[subject])],
                walltime="80:00:00",
                memory=mem_demand[subject]*10 -10,
                nodes=1,
                tasks=mem_demand[subject] -1,
                env="mne",
                name="dcd_" + area + subject,
            )


# This function returns labels for one subject
def get_labels(subject, only_glasser):
    if not only_glasser:
        labels = pymegsr.get_labels(
            subject=subject,
            filters=["*wang*.label", "*JWDG*.label"],
            annotations=["HCPMMP1"],
            sdir=subjects_dir,
        )
        labels = pymegsr.labels_exclude(
            labels=labels,
            exclude_filters=[
                "wang2015atlas.IPS4",
                "wang2015atlas.IPS5",
                "wang2015atlas.SPL",
                "JWDG_lat_Unknown",
            ],
        )
        labels = pymegsr.labels_remove_overlap(
            labels=labels, priority_filters=["wang", "JWDG"]
        )
    else:
        labels = pymegsr.get_labels(
            subject=subject,
            filters=["select_nothing"],
            annotations=["HCPMMP1"],
            sdir=subjects_dir,
        )
    return labels


#@memory.cache
def decode(
    subject,
    area,
    sessinfo,
    epoch_type="stimulus",
    only_glasser=False,
    BEM="three_layer",
    debug=False,
    target="response",
):
    # Only show warning and above:
    mne.set_log_level("WARNING")
    pymeglcmv.logging.getLogger().setLevel(logging.INFO)

    set_n_threads(1)

    # Get all labels for this subject
    labels = get_labels(subject, only_glasser)
    # And keep only the ones that we are interested in (area parameter of this function) -> needs to be in areas_to_labels dict defined above
    labels = [x for x in labels if any([cl for cl in areas_to_labels[area] if cl in x.name])]
    print(labels)

    if len(labels) < 1:
        raise RuntimeError('Expecting at least two labels')

    # Turn individual labels into one big label that contains all vertices that
    # belong to one ROI
    label = labels.pop()
    for l in labels:
        label += l

    print('Selecting this label for area %s:'%area, label)

    ###### PM ATTEMPT TAKING SESS/REC INFO & EXISTING FUNCTIONS INTO ACCOUNT
    data=[]; fwds=[]; bems=[]; sources=[];
    low_pass_fs = None   # low-pass filter cutoff; set to None is no filter required
    for sess, rec in sessinfo:
        logging.info("Reading data for %s, sess %i, rec %i "% (subject,sess,rec))
        data_cov,epoch,epoch_filename = get_stim_epoch(subject, sess, rec, lopass=low_pass_fs)
        data.append((data_cov, epoch));

        logging.info("Setting up source space and forward model")
        raw_filename = glob('/home/pmurphy/meg_data/surprise/%s-%i*_0%i.ds' %
                            (subject, sess, rec))[0]
        trans_filename = glob('/home/pmurphy/meg_data/surprise/MRIs/trans_mats/%s_%i_0%i*fif' % (
            subject, sess, rec))[0]
        forward, bem, source = get_leadfield(
            subject, raw_filename, epoch_filename, trans_filename, bem_sub_path='bem_ft')
        fwds.append(forward); # bems.append(bem); sources.append(source);

    # pull trial IDs
    events = [d[1].events[:, 2] for d in data]
    events = np.hstack(events)

    # Compute LCMV filters for each session
    filters = []
    for (data_cov, epochs), forward in zip(data, fwds):
        filters.append(
            pymeglcmv.setup_filters(epochs.info, forward, data_cov, None, [label])
        )
    set_n_threads(1)

    # specify decoding settings
    clf = Pipeline(
        [
           ("Scaling", StandardScaler()),
           ("PCA", PCA(n_components=0.95, svd_solver='full')),
           ("RidgeReg", linear_model.Ridge(alpha=10)),
        ]
    )

    # specify sample onsets and window for decoding
    smpon = np.arange(0.4,5,0.4)   # vector of sample onsets (s)
    smpwin = [-0.10001, 1.40001]   # window for decoding, relative to sample onset (s) - going marginally outside desired bnds important to catch all times
    basewin = [-0.075, 0]   # window relative to sample onset for sample-specific baseline (set to None if not required)
    subsmp = 200    # frequency (Hz) at which to subsample (set to None if no subsampling required - original fs = 400 Hz); NB: should be factor of 400

    # load to-be-decoded variables & check that trials are appropiately aligned
    matname = ('/home/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed4mne/BehavPupil/%s_4decode.mat' % (subject))
    mat = sio.loadmat(matname)

    mat_events = np.int64(np.concatenate(mat["tIDs"]))  # convert matlab events to same type as python events
    assert np.array_equal(events,mat_events[:len(events)])

    # Perform source reconstruction, using for each session the appropriate filter
    # Iterates over sample positions to mitigate memory demands
    all_smp = []  # inialize DataFrame containing all decoding results, across sample positions/variables
    for smp in range(len(smpon)):
        fname = "/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/decodeERF/SSbase_nofilt/%s_%s_%s_finegrainERF.hdf" % (subject, area, str(smp+1))
        # the try: except: block implements caching, if output is already there, don't do it again.
        try:
            all_s = pd.read_hdf(fname)
        except FileNotFoundError:
            # perform source reconstruction of ERF data
            erfdata, events, times = decoding_ERF.get_lcmv(
                [d[1].copy().crop(smpon[smp]+smpwin[0],smpon[smp]+smpwin[1]) for d in data], filters, njobs=6    # erfdata is ntrials*nvertices*ntimes
            )
            times = times-smpon[smp]  # making times vector relative to sample onset

            # sample-specific baseline subtraction
            if basewin is not None:
                id_time = (basewin[0] <= times) & (times <= basewin[1])
                means = erfdata[:, :, id_time].mean(-1)
                erfdata -= means[:, :, np.newaxis]

            # subsample
            if subsmp is not None:
                erfdata = erfdata[:,:,0::round(400/subsmp)]
                times = times[0::round(400/subsmp)]

            # loop through variables to be decoded
            all_s = []
            for target in ["LLR", "vPE"]:
                # pull variable to be decoded
                target_vals = mat[target]   #  target_vals will be a numpy ndarray, ntrials*nsamples
                target_vals = target_vals[:len(events),:]

                # perform decoding
                dcd = decoding_ERF.Decoder(target_vals[:,smp],("RidgeReg",clf))
                k = dcd.classify(
                    erfdata, times, events, area,   # feeding in times aligned to smp onset
                    average_vertices=False
                )
                k.loc[:, "target"] = target  # include target_val label in dataframe
                all_s.append(k)   #  append decoding results for this target_val combo

            all_s = pd.concat(all_s)         # concatenate all target_vals
            all_s.loc[:, 'ROI'] = area       # and include ROI/sample position labels
            all_s.loc[:, "sample"] = str(smp+1)
            all_s.to_hdf(fname, "df")  # save once all target_vals have been iterated over

        all_smp.append(all_s)  # append decoding results for this sample position

    all_smp = pd.concat(all_smp)  # concatenate all sample positions
    all_smp.to_csv("/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/decodeERF/SSbase_nofilt/%s_%s_full_finegrainERF.csv" % (subject, area))

    return k


