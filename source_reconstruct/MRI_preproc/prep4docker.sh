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

xhemireg --s ${IDs[var11]} 1> /home/pmurphy/meg_data/surprise/MRIs/fs_converted/Log_files/"$PBS_JOBID"1.out 2> /home/pmurphy/meg_data/surprise/MRIs/fs_converted/Log_files/"$PBS_JOBID"1.err
surfreg --s ${IDs[var11]} --t fsaverage_sym --lh 1> /home/pmurphy/meg_data/surprise/MRIs/fs_converted/Log_files/"$PBS_JOBID"1.out 2> /home/pmurphy/meg_data/surprise/MRIs/fs_converted/Log_files/"$PBS_JOBID"1.err
surfreg --s ${IDs[var11]} --t fsaverage_sym --lh --xhemi 1> /home/pmurphy/meg_data/surprise/MRIs/fs_converted/Log_files/"$PBS_JOBID"1.out 2> /home/pmurphy/meg_data/surprise/MRIs/fs_converted/Log_files/"$PBS_JOBID"1.err



