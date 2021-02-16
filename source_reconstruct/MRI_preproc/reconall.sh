#!/bin/sh

#PBS -q batch
#PBS -l walltime=500:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=14gbs

# -- run in the current working (submission) directory --
cd $PBS_O_WORKDIR

chmod g=wx $PBS_JOBNAME

export SUBJECTS_DIR=/home/pmurphy/meg_data/surprise/MRIs/fs_converted/

IDs=( 'KSV' 'EXJ' 'TSJ' 'JTB' 'EXF' 'ECB' 'EMB' 'TFD' 'GSB' 'EXG' 'OMF' 'NIF' 'DHB' 'HBC' 'QNV' 'DCB' 'TNB' 'PDP')
MRIs=(13142 14412 14588 14684 15565 16049 16395 16631 17699 18106 18146 18656 18657 18911 20020 20042 20365 00000)

recon-all -subject ${IDs[var11]} -i /home/pmurphy/meg_data/surprise/MRIs/nii_converted/*${MRIs[var11]}*.gz -all 1> /home/pmurphy/meg_data/surprise/MRIs/fs_converted/Log_files/"$PBS_JOBID"1.out 2> /home/pmurphy/meg_data/surprise/MRIs/fs_converted/Log_files/"$PBS_JOBID"1.err

