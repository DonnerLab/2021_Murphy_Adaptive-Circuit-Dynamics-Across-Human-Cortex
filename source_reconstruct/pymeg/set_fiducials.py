"""
Script for initializing pymeg make_trans function for making transformation matrix
"""

import os

os.environ["SUBJECTS_DIR"] = "/mnt/homes/home024/pmurphy/meg_data/surpriseD/MRIs/fs_converted/"

import mne
mne.gui.coregistration(subjects_dir='/mnt/homes/home024/pmurphy/meg_data/surpriseD/MRIs/fs_converted/')
