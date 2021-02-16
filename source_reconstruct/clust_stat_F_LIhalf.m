function clust_stat_F_LIhalf(nclust)

allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};
clusters = {'vfcPrimary';'vfcEarly';'vfcVO';'vfcPHC';'vfcTO';'vfcLO';'vfcV3ab';'vfcIPS01';'vfcIPS23';'vfcFEF';...   % Wang
    'JWG_aIPS';'JWG_IPS_PCeS';'JWG_M1';...                                                                          % JW
    'HCPMMP1_cingulate_pos';'HCPMMP1_paracentral_midcingulate';'HCPMMP1_insular_front_opercular';'HCPMMP1_premotor';'HCPMMP1_dlpfc';'HCPMMP1_frontal_inferior';'HCPMMP1_frontal_orbital_polar';... % Glasser
    'post_medial_frontal';'ant_medial_frontal';'vent_medial_frontal'};  % PM-defined medial PFC clusters

cnames = {'V1','V2-V4','VO1/2','PHC','MT+','LO1/2','V3A/B','IPS0/1','IPS2/3','FEF',...
    'aIPS','IPS/PCeS','M1',...
    'PCC','MCC','Insula','Premotor','dlPFC','vlPFC','OFC',...
    'pmFC','amPFC','vmPFC'};

induced = 1;
modeltype = 'fitted_linIU';  % switch b/w 'normative', 'fitted', 'fitted_np', 'fitted_lin', 'fitted_npIU' & 'fitted_linIU'
priortype = 'psi';  % switch between 'LPR' (prior term will be posterior from previous sample) & 'psi' (prior term will be Glaze-transformed prior for current sample)
surprisetype = 'pCP';  % switch between 'pCP' and 'surprise'

coeftype = 'beta';   % switch b/w 'beta' and 'tscore'

addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults

if induced, str1 = 'Conv2mne_induced'; else str1 = 'Conv2mne'; end
if strcmp(modeltype,'normative'), str2 = 'output';
elseif strcmp(modeltype,'fitted'), str2 = 'output_fitted';
elseif strcmp(modeltype,'fitted_np'), str2 = 'output_fitted_np';
elseif strcmp(modeltype,'fitted_lin'), str2 = 'output_fitted_lin';
elseif strcmp(modeltype,'fitted_npIU'), str2 = 'output_fitted_npIU';
elseif strcmp(modeltype,'fitted_linIU'), str2 = 'output_fitted_linIU'; end
loadpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/',str1,'/agg/',str2,'/'];
savepath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/',str1,'/agg/',str2,'/av/'];

foi = [1 100];  % bounds on frequencies for plotting/cluster-corrected stats
clustalpha = 0.05;  % % intial threhsold for forming clusters

trlwin = [-0.5 5.8];
smpwin = [0 1.0];

if strcmp(surprisetype,'pCP'), sstr='_pCP'; else sstr=''; end
if strcmp(coeftype,'beta'), bstr='_beta'; else bstr=''; end

% Concatenate data across subjects
samps = 2:12;  % sub-selection of samples for certain analyses (set 2:12 for all samples)
fprintf('Loading data...\nSubject ')
for s = 1:length(allsubj)
    fprintf('%d, ',s)
    if strcmp(priortype,'psi')
        load([loadpath,allsubj{s},'_F_',clusters{nclust},'_psi_LIhalf',sstr,bstr,'.mat'])  % load subject-specific results for this cluster
        xtra = load([loadpath,allsubj{s},'_F_',clusters{nclust},'_psi_LIhalf',sstr,bstr,'.mat'],'freqs','times');  % nonsense couple of lines to deal with freqs.m and times.m function/naming conflicts
    else
        load([loadpath,allsubj{s},'_F_',clusters{nclust},'_LIhalf',bstr,sstr,'.mat'])
        xtra = load([loadpath,allsubj{s},'_F_',clusters{nclust},'_LIhalf',bstr,sstr,'.mat'],'freqs','times');
    end
    
    % Store data
    gaB_priorS_MR_wC(:,:,s) = squeeze(mean(TpriorS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3));
    gaB_llrS_MR_wC(:,:,s) = squeeze(mean(TllrS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3));
    gaB_llrXsurpriseS_MR_wC(:,:,s) = squeeze(mean(TllrXsurpriseS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3));
    gaB_llrXcertaintyS_MR_wC(:,:,s) = squeeze(mean(TllrXcertaintyS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3));
    gaB_llrXpupilS_MR_wC(:,:,s) = squeeze(mean(TllrXpupilS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),samps-1),3));
            
    % Add structures to group-level arrays
    cstruct = struct;
    cstruct.label = {cnames{nclust}};  % arbitrary label for motor lateralization 'channel'
    cstruct.fsample = 20;  % sampling rate
    cstruct.freq = freq(freq>=foi(1) & freq<=foi(2));
    cstruct.time = smptimes(smptimes>=smpwin(1) & smptimes<=smpwin(2));
    cstruct.dimord = 'chan_freq_time';
    temp_size = [1 size(squeeze(mean(TpriorS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3)))];
        
    cstruct.powspctrm = reshape(mean(TpriorS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3),temp_size); allB_priorS_MR_wC{s} = cstruct;
    cstruct.powspctrm = reshape(mean(TllrS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3),temp_size); allB_llrS_MR_wC{s} = cstruct;
    cstruct.powspctrm = reshape(mean(TllrXsurpriseS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3),temp_size); allB_llrXsurpriseS_MR_wC{s} = cstruct;
    cstruct.powspctrm = reshape(mean(TllrXcertaintyS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3),temp_size); allB_llrXcertaintyS_MR_wC{s} = cstruct;
    cstruct.powspctrm = reshape(mean(TllrXpupilS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3),temp_size); allB_llrXpupilS_MR_wC{s} = cstruct;
    cstruct.powspctrm = reshape(mean(TllrXsurpriseS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3)-...
                                   mean(TllrXcertaintyS_MR_wC(freq>=foi(1) & freq<=foi(2),smptimes>=smpwin(1) & smptimes<=smpwin(2),:),3),temp_size); allB_surpMcert_MR_wC{s} = cstruct;
        
    cstruct.powspctrm = zeros(size(cstruct.powspctrm)); allBnull{s} = cstruct;
end
fprintf('done.\n')

% Set up cfg stucture for cluster-based permutation tests
cfg_stat = [];
cfg_stat.latency     = 'all';
cfg_stat.frequency   = 'all';
cfg_stat.avgoverchan = 'yes';
cfg_stat.avgovertime = 'no';
cfg_stat.avgoverfreq = 'no';
cfg_stat.parameter   = 'powspctrm';
cfg_stat.method      = 'montecarlo';
cfg_stat.statistic   = 'depsamplesT';  % also possibly try 'depsamplesregrT'
cfg_stat.alpha       = 0.05;
cfg_stat.correctm    = 'cluster';
cfg_stat.clusteralpha = clustalpha;  % intial threhsold for forming clusters
cfg_stat.clusterstatistic = 'maxsum';  % method for quantifying combined cluster test statistics
cfg_stat.correcttail = 'prob';
cfg_stat.numrandomization = 10000;

Nsub = length(allsubj);
cfg_stat.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg_stat.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg_stat.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg_stat.uvar                = 2; % the 2nd row in cfg.design contains the subject number

% Calculate cluster-corrected stats
stat1 = ft_freqstatistics(cfg_stat,allB_priorS_MR_wC{:},allBnull{:});
stat2 = ft_freqstatistics(cfg_stat,allB_llrS_MR_wC{:},allBnull{:});
stat3 = ft_freqstatistics(cfg_stat,allB_llrXsurpriseS_MR_wC{:},allBnull{:});
stat4 = ft_freqstatistics(cfg_stat,allB_llrXcertaintyS_MR_wC{:},allBnull{:});
stat5 = ft_freqstatistics(cfg_stat,allB_llrXpupilS_MR_wC{:},allBnull{:});
stat6 = ft_freqstatistics(cfg_stat,allB_surpMcert_MR_wC{:},allBnull{:});

% Save
if strcmp(priortype,'psi')
    savename = [savepath,'Clust_stats_F_',clusters{nclust},'_psi_LIhalf',sstr,bstr,'.mat'];
else 
    savename = [savepath,'Clust_stats_F_',clusters{nclust},'_LIhalf',sstr,bstr,'.mat'];
end
save(savename,'stat*','ga*','all_trl_avg','all_smp_avg','cfg_stat')

