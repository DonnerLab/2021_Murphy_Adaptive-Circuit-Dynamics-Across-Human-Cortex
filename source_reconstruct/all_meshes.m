addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts

allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXJ' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};

% make_mesh('S03', '/home/nwilming/fs_subject_dir')

for i = 2:length(allsubj)
    fprintf('\n\n ------------------ \n Starting subject %s \n ------------------ \n\n',allsubj{i})
   make_mesh(allsubj{i}, '/home/pmurphy/meg_data/surprise/MRIs/fs_converted');
end