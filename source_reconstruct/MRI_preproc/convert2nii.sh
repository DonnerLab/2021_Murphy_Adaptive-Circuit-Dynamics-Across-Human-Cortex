#!/bin/sh

# Submits n jobs to the torque queing system

# for i in 13142 14412 14588 14684 15565 16049 16222 16395 16631 16873 17699 18106 18146 18656 18657 18911 20020 20042 20365
for i in 20020 20042 20365
do
  
  echo 'Start file ' $i
  
  /home/pmurphy/bin/dcm2niix -z y -o /home/pmurphy/meg_data/surprise/MRIs/nii_converted /home/pmurphy/meg_data/surprise/MRIs/*$i/*6
  
done
