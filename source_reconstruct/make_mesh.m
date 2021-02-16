function make_mesh(subject, fs_path)
% subject = 'S01';
% fs_path = '/home/nwilming/fs_subject_dir';
addpath('/home/nwilming/spm12/')
% addpath('/home/nwilming/fieldtrip-20170914/')
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'

ft_defaults
mri_path = fullfile(fs_path, subject, 'mri/T1.mgz');
mri = ft_read_mri(mri_path);
mri2 = mri;
mri2.transform = mri.hdr.tkrvox2ras;
mri2.coordsys = 'spm';
cfg = [];
cfg.output = {'brain','skull','scalp'};
segmentedmri = ft_volumesegment(cfg,mri2);
cfg = [];
cfg.tissue={'scalp','skull','brain'};

n = (2 + (81920/2))/8;
cfg.numvertices = [n, n, n];
bnd=ft_prepare_mesh(cfg,segmentedmri);

names = {'outer_skin.surf', 'outer_skull.surf', 'inner_skull.surf'};
bem_dir = fullfile(fs_path, subject, 'bem_ft');
mkdir(bem_dir)
for ii = 1:3
    path = fullfile(bem_dir, names{ii});
    ft_write_headshape(path, bnd(ii), 'format', 'freesurfer');
end
end