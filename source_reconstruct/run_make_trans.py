"""
Script for initializing pymeg make_trans function for making transformation matrix
"""

import os

os.environ["PYMEG_CACHE_DIR"] = "/mnt/homes/home024/pmurphy/tmp"
os.environ["SUBJECTS_DIR"] = "/mnt/homes/home024/pmurphy/meg_data/surprise/MRIs/fs_converted/"

subj = input("Subject ID: ")
sess = input("Session #: ")
rec = input("Recording #: ")

epochs = "%s%s%s%s%s%s%s" % ("/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/", subj, "-", sess, "_", rec, "_preproc4mne.mat-epo.fif.gz")

import glob
rawfile = glob.glob("%s%s%s%s%s%s%s" % ("/mnt/homes/home024/pmurphy/meg_data/surprise/", subj, "-", sess, "*", rec, ".ds"))
rawfile = rawfile[0]

from pymeg import source_reconstruction as sr
sr.make_trans(subj, rawfile, epochs, "%s%s%s%s%s" % (subj,"-",sess,"_",rec))
