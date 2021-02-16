# Some imports:
import logging
import mne
import numpy as np
import os
import scipy.io as sio
import time

from joblib import Memory # Provides caching of results

from os import makedirs
from os.path import join
from glob import glob

from pymeg import lcmv as pymeglcmv
from pymeg import source_reconstruction as pymegsr
from pymeg import decoding

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn import linear_model

import datetime
import pandas as pd

from pymeg.lcmv_peter import get_stim_epoch
from pymeg.source_reconstruction import (
    get_leadfield,
    make_trans,   # NB: source_reconstruction also has get_trans_epoch, get_ctf_trans - not sure which we want, NW has get_trans
)

# Setup some paths:
memory = Memory(cachedir='/mnt/homes/home024/pmurphy/tmp/')
# path = "/home/nwilming/conf_meg/sr_labeled/"    ###### PM: WHAT IS THIS PATH FOR?????
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
subjects = {'OMF': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'PDP': [(1, 1), (1, 2), (2, 2), (2, 3), (3, 1), (3, 2)],
            'QNV': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2)],
            'TFD': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)],
            'TNB': [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)],
            'TSJ': [(1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2)]}

# typical processing demands (in # cores) per subject
#mem_demand = {'DCB': 6, 'DHB': 4, 'ECB': 4, 'EMB': 3, 'EXF': 3, 'EXG': 4,  # SETTINGS FOR VERTEX-PRESERVED DECODING
#              'GSB': 3, 'HBC': 3, 'JTB': 3, 'KSV': 3, 'NIF': 3, 'OMF': 3,
#              'PDP': 4, 'QNV': 7, 'TFD': 4, 'TNB': 4, 'TSJ': 7}
mem_demand = {'DCB': 3, 'DHB': 2, 'ECB': 2, 'EMB': 2, 'EXF': 2, 'EXG': 2,   # SETTINGS FOR VERTEX-AVERAGED DECODING
              'GSB': 2, 'HBC': 2, 'JTB': 2, 'KSV': 2, 'NIF': 2, 'OMF': 2,
              'PDP': 2, 'QNV': 3, 'TFD': 2, 'TNB': 2, 'TSJ': 3}

# typical processing demands (in # vertices) per area (not used but useful reference fo setting ntasks at submit)
mem_area = {'vfcPrimary': 3650, 'vfcEarly': 9610, 'vfcV3ab': 3280,
            'vfcIPS01': 4200, 'vfcIPS23': 2200, 'JWG_aIPS': 870,
            'JWG_IPS_PCeS': 4600, 'JWG_M1': 2900, 'HCPMMP1_premotor': 9900}

# Mapping areas to labels
areas_to_labels = {
    "vfcPrimary": [
        u"lh.wang2015atlas.V1d-lh", u"rh.wang2015atlas.V1d-rh",
        u"lh.wang2015atlas.V1v-lh", u"rh.wang2015atlas.V1v-rh",
    ],
    #"vfcEarly": [
    #    u"lh.wang2015atlas.V2d-lh", u"rh.wang2015atlas.V2d-rh",
    #    u"lh.wang2015atlas.V2v-lh", u"rh.wang2015atlas.V2v-rh",
    #    u"lh.wang2015atlas.V3d-lh", u"rh.wang2015atlas.V3d-rh",
    #    u"lh.wang2015atlas.V3v-lh", u"rh.wang2015atlas.V3v-rh",
    #    u"lh.wang2015atlas.hV4-lh", u"rh.wang2015atlas.hV4-rh",
    #],
    # "vfcVO": [
    #     u"lh.wang2015atlas.VO1-lh", u"rh.wang2015atlas.VO1-rh",
    #     u"lh.wang2015atlas.VO2-lh", u"rh.wang2015atlas.VO2-rh",
    # ],
    # "vfcPHC": [
    #     u"lh.wang2015atlas.PHC1-lh", u"rh.wang2015atlas.PHC1-rh",
    #     u"lh.wang2015atlas.PHC2-lh", u"rh.wang2015atlas.PHC2-rh",
    # ],
    #"vfcV3ab": [
    #    u"lh.wang2015atlas.V3A-lh", u"rh.wang2015atlas.V3A-rh",
    #    u"lh.wang2015atlas.V3B-lh", u"rh.wang2015atlas.V3B-rh",
    #],
    # "vfcTO": [
    #     u"lh.wang2015atlas.TO1-lh", u"rh.wang2015atlas.TO1-rh",
    #     u"lh.wang2015atlas.TO2-lh", u"rh.wang2015atlas.TO2-rh",
    # ],
    # "vfcLO": [
    #     u"lh.wang2015atlas.LO1-lh", u"rh.wang2015atlas.LO1-rh",
    #     u"lh.wang2015atlas.LO2-lh", u"rh.wang2015atlas.LO2-rh",
    # ],
    #"vfcIPS01": [
    #    u"lh.wang2015atlas.IPS0-lh", u"rh.wang2015atlas.IPS0-rh",
    #    u"lh.wang2015atlas.IPS1-lh", u"rh.wang2015atlas.IPS1-rh",
    #],
    #"vfcIPS23": [
    #    u"lh.wang2015atlas.IPS2-lh", u"rh.wang2015atlas.IPS2-rh",
    #    u"lh.wang2015atlas.IPS3-lh", u"rh.wang2015atlas.IPS3-rh",
    #],
    #"JWG_aIPS": ["lh.JWDG.lr_aIPS1-lh", "rh.JWDG.lr_aIPS1-rh", ],
    #"JWG_IPS_PCeS": ["lh.JWDG.lr_IPS_PCes-lh", "rh.JWDG.lr_IPS_PCes-rh", ],
    #"JWG_M1": ["lh.JWDG.lr_M1-lh", "rh.JWDG.lr_M1-rh", ],
    #"HCPMMP1_premotor": (
    #    ["L_{}_ROI-lh".format(area) for area in [
    #        "55b", "6d", "6a", "FEF", "6v", "6r", "PEF"]] +
    #    ["R_{}_ROI-rh".format(area) for area in [
    #        "55b", "6d", "6a", "FEF", "6v", "6r", "PEF"]]
    #),
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
                [(subject, area, subjects[subject])],  ##### PM: added session/recording info as input
                walltime="200:00:00",
                memory= mem_demand[subject]*10 + 10,
                nodes=1,
                tasks= mem_demand[subject] + 1,
                env="mne",   ##### PM: switched from py36 to mne
                name="dcd_" + area + subject,
            )
            time.sleep(1)


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
    for sess, rec in sessinfo:
        logging.info("Reading data for %s, sess %i, rec %i "% (subject,sess,rec))
        data_cov,epoch,epoch_filename = get_stim_epoch(subject, sess, rec)  #### PM: my get_stim_epoch fun may be different to NW's, has 3 outputs...
        data.append((data_cov, epoch)); ##### PM: N.B. data may in fact need to be output currently called 'epochs'

        logging.info("Setting up source space and forward model")
        raw_filename = glob('/home/pmurphy/meg_data/surprise/%s-%i*_0%i.ds' %
                            (subject, sess, rec))[0]
        trans_filename = glob('/home/pmurphy/meg_data/surprise/MRIs/trans_mats/%s_%i_0%i*fif' % (
            subject, sess, rec))[0]
        forward, bem, source = get_leadfield(                                 ########### PM: this is my get_leadfield function - may still need work to be integrated here
            subject, raw_filename, epoch_filename, trans_filename, bem_sub_path='bem_ft')
        fwds.append(forward); # bems.append(bem); sources.append(source);


    # Define TFR parameters
    fois = np.arange(35, 101, 5)   # PM: np.arange(36, 162, 4)
    #lfois = np.hstack([np.arange(1, 11, 1),np.arange(12, 31, 2)])    # PM: np.arange(1, 36, 1)
    lfois = np.hstack([np.arange(1, 7, 1),np.arange(16, 31, 2)])   # NB: THIS SETTING EXCLUDES 7-15 Hz!!!!
    tfr_params = {
        "HF": {               # PM changed from 'F' to 'HF'
            "foi": fois,
            "cycles": fois * 0.25,  # PM: fois * 0.25
            "time_bandwidth": 6 + 1,   # PM: 6 + 1
            "n_jobs": 1,
            "est_val": fois,
            "est_key": "HF",   # PM changed from 'F' to 'HF'
            "sf": 400,         # PM added
            "decim": 20,       # PM added
        },
        "LF": {
            "foi": lfois,
            "cycles": lfois * 0.4,  # PM: fois * 0.4
            "time_bandwidth": 1 + 1,     # PM: 1 + 1
            "n_jobs": 1,
            "est_val": lfois,
            "est_key": "LF",
            "sf": 400,         # PM added
            "decim": 20,       # PM added
        },
    }

    events = [d[1].events[:, 2] for d in data]
    events = np.hstack(events)

    # Compute LCMV filters for each session
    filters = []
    for (data_cov, epochs), forward in zip(data, fwds):   # PM: N.B. data_cov may not be part of data as currently configured
        filters.append(
            pymeglcmv.setup_filters(epochs.info, forward, data_cov, None, [label])
        )
    set_n_threads(1)
    
    # Specify vertex -> hemisphere mapping array --- COMMENT IN IF WANT TO AVERAGE ACROSS VERTICES
    f = filters[0][label.name]
    avg_vertices = np.zeros((len(f['vertices'][0]) + len(f['vertices'][1]))).astype(bool)
    avg_vertices[:len(f['vertices'][0])] = True
    
    # specify decoding settings
    clf = Pipeline(
        [
           ("Scaling", StandardScaler()),
           ("PCA", PCA(n_components=0.95, svd_solver='full')),
           ("RidgeReg", linear_model.Ridge(alpha=1)),
        ]
    )

    # specify sample onsets and window for decoding
    smpon = np.arange(0.4,5,0.4)   # vector of sample onsets (s)
    smpwin = [-0.10001, 1.40001]   # window for decoding, relative to sample onset (s) - going marginally outside desired bnds important to catch all times

    # load to-be-decoded variables & check that trials are appropiately aligned
    matname = ('/home/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed4mne/BehavPupil/%s_4decode.mat' % (subject))
    mat = sio.loadmat(matname)

    mat_events = np.int64(np.concatenate(mat["tIDs"]))  # convert matlab events to same type as python events
    assert np.array_equal(events,mat_events[:len(events)])    ##### PM: remove [:len(events)] after testing single session

    # Perform source reconstruction, using for each session the appropriate filter
    # Iterates over sample positions to mitigate memory demands
    all_smp = []  # inialize DataFrame containing all decoding results, across sample positions/variables
    for smp in range(len(smpon)):
        fname = "/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/decodeAv_nophase_noalpha/%s_%s_%s_avTF.hdf" % (subject, area, str(smp+1))
        # the try: except: block implements caching, if output is already there, don't do it again.
        try:
            all_s = pd.read_hdf(fname)
        except FileNotFoundError:
            # perform source reconstruction of TF data
            HF_tfrdata, events, HF_freq, times = decoding.get_lcmv(   # PM: padding by 0.2s (max TF win / 2) for accurate TF estimation
                tfr_params["HF"], [d[1].copy().crop(smpon[smp]+smpwin[0]-0.2,smpon[smp]+smpwin[1]+0.2) for d in data], filters, njobs=6    #d[1].copy().crop() pulls out sample-aligned data
            )
            LF_tfrdata, events, LF_freq, times = decoding.get_lcmv(
                tfr_params["LF"], [d[1].copy().crop(smpon[smp]+smpwin[0]-0.2,smpon[smp]+smpwin[1]+0.2) for d in data], filters, njobs=6
            )

            # Concatenate data
            tfrdata = np.hstack([HF_tfrdata, LF_tfrdata])
            del LF_tfrdata, HF_tfrdata
            freq = np.concatenate([HF_freq, LF_freq])

            # loop through variables to be decoded
            ctimes = (smpwin[0]+smpon[smp] <= times) & (times <= smpwin[1]+smpon[smp])  # indices of time-points without padding
            all_s = []
            for target in ["LLR", "vPE"]:
                # pull variable to be decoded
                target_vals = mat[target]   #  target_vals will be a numpy ndarray, ntrials*nsamples
                target_vals = target_vals[:len(events),:]    ##### PM: remove this line after testing single session

                # perform decoding
                dcd = decoding.Decoder(target_vals[:,smp],("RidgeReg",clf))
                k = dcd.classify(
                    tfrdata[:,:,:,ctimes], times[ctimes]-smpon[smp], freq, events, area,   # feeding in times aligned to smp onset
                    average_vertices=avg_vertices, use_phase=False            ####### PM: NB, set average_vertices to False if want to preserve vertices as separate features
                )
                k.loc[:, "target"] = target  # include target_val label in dataframe
                all_s.append(k)   #  append decoding results for this target_val combo

            all_s = pd.concat(all_s)         # concatenate all target_vals
            all_s.loc[:, 'ROI'] = area       # and include ROI/sample position labels
            all_s.loc[:, "sample"] = str(smp+1)
            all_s.to_hdf(fname, "df")  # save once all target_vals have been iterated over

        all_smp.append(all_s)  # append decoding results for this sample position

    all_smp = pd.concat(all_smp)  # concatenate all sample positions
    all_smp.to_csv("/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/decodeAv_nophase_noalpha/%s_%s_full_avTF.csv" % (subject, area))

    return k

# # code for reading already created individual sample files to csv
# subject = "ECB"
# area = area = "vfcPrimary"
# for alphaRR in [1]:
#     all_smp = []
#     for smp in [0]:
#         fname = "/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/decode/%s_%s_%s_%s_finegrainTF_nophase.hdf" % (subject, area, str(smp+1), str(alphaRR))
#         all_s = pd.read_hdf(fname)
#         all_smp.append(all_s)
#     all_smp = pd.concat(all_smp)  # concatenate all sample positions
#     all_smp.to_csv("/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/decode/%s_%s_%s_full_finegrainTF_nophase.csv" % (subject, area, str(alphaRR)))
