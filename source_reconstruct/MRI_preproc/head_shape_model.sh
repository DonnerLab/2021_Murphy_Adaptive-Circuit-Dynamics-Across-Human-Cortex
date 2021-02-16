# This script makes subject-specific high density head shape models via MNE

IDs=( 'KSV' 'EXJ' 'TSJ' 'JTB' 'EXF' 'ECB' 'EMB' 'TFD' 'GSB' 'EXG' 'OMF' 'NIF' 'DHB' 'HBC' 'DCB' 'TNB' 'PDP' 'QNV')

for s in {1..16}
do
    
export subjid=${IDs[s]}
export SUBJECTS_DIR=/home/pmurphy/meg_data/surprise/MRIs/fs_converted

mne make_scalp_surfaces -s ${subjid} -d ${SUBJECTS_DIR}

done
