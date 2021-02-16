clear, close all

% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %
% --- %% CHOICE LATERALIZATION & SAMPLE-WISE REGRESSIONS %% --- %
% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %
% Specify inits
allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};

basetype = 'dB_common';  % options: 'dB_st', 'dB_common', 'pc_st' or 'pc_common'
modeltype = 'fitted_linIU';  % switch b/w 'normative', 'fitted', 'fitted_np' & 'fitted_lin'
surprisetype = 'pCP';  % switch between 'pCP' and 'surprise'

coeftype = 'beta';   % switch b/w 'beta' and 'tscore'

addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
addpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts/FMINSEARCHBND')
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/agg/output/av/'))

if strcmp(modeltype,'normative'), str2 = 'output';
elseif strcmp(modeltype,'fitted'), str2 = 'output_fitted';
elseif strcmp(modeltype,'fitted_np'), str2 = 'output_fitted_np';
elseif strcmp(modeltype,'fitted_lin'), str2 = 'output_fitted_lin';
elseif strcmp(modeltype,'fitted_npIU'), str2 = 'output_fitted_npIU';
elseif strcmp(modeltype,'fitted_linIU'), str2 = 'output_fitted_linIU'; end
loadpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/TF/wMotorLoc/',basetype,filesep,str2,filesep];
savepath = [loadpath,'av',filesep];

clustalpha = 0.05;  % % intial threhsold for forming clusters
nperms = 10000;

if strcmp(surprisetype,'pCP'), surpstr = '_pCP'; else surpstr = ''; end
if strcmp(coeftype,'beta'), bstr = '_beta'; else bstr = ''; end

%%%%%%%%%%%%%%%%%%
%%% Trial-wise %%%
%%%%%%%%%%%%%%%%%%
fprintf('Loading data...\nSubject ')

for s = 1:length(allsubj)
    fprintf('%d, ',s)
    % Load data
    load([loadpath,allsubj{s},'_trialwise_output_appML',surpstr,bstr,'.mat'],'freqs','tf_avg','TposteriorTF')
    freqsL = freqs;
        
    trl_avgL = tf_avg;
    TposteriorL = TposteriorTF;
    
    load([loadpath,allsubj{s},'_trialwise_output_appML_HiFreq',bstr,'.mat'],'freqs','tf_avg','TposteriorTF','trltimes')
    allfreqs = [freqsL freqs];
        
    % Store data for grand-average, sensor-resolved TFR
    trl_avg(:,:,s) = [trl_avgL; tf_avg];
    trl_Tposterior(:,:,s) = [TposteriorL; TposteriorTF];

    % Create structure for this subject to mimic output from FT_FREQANALYSIS
    cstruct = struct;
    cstruct.label = {'LI'};  % arbitrary label for motor lateralization 'channel'
    cstruct.fsample = 20;  % sampling rate
    cstruct.freq = allfreqs;
    cstruct.time = trltimes;
    cstruct.dimord = 'chan_freq_time';
    
    % Add structures to group-level arrays
    cstruct.powspctrm = reshape([trl_avgL; tf_avg],[1,size(trl_avg,1),size(trl_avg,2)]); all_trl_avg{s} = cstruct;
    cstruct.powspctrm = reshape([TposteriorL; TposteriorTF],[1,size(trl_Tposterior,1),size(trl_Tposterior,2)]); all_trl_Tposterior{s} = cstruct;
    cstruct.powspctrm = zeros(size(cstruct.powspctrm));
    all_trl_null{s} = cstruct;
end
fprintf('done.\nRunning stats.\n')

% Run cluster-based permutation tests
cfg_stat = [];
cfg_stat.latency     = 'all';
cfg_stat.frequency   = 'all';
cfg_stat.avgoverchan = 'no';
cfg_stat.avgovertime = 'no';
cfg_stat.avgoverfreq = 'no';
cfg_stat.parameter   = 'powspctrm';
cfg_stat.method      = 'montecarlo';
cfg_stat.statistic   = 'depsamplesT';  % also possibly try 'depsamplesregrT'
cfg_stat.alpha       = 0.05;
cfg_stat.correctm    = 'cluster';
% cfg_stat.correctm    = 'no';
cfg_stat.clusteralpha = clustalpha;  % intial threhsold for forming clusters
cfg_stat.clusterstatistic = 'maxsum';  % method for quantifying combined cluster test statistics
cfg_stat.correcttail = 'prob';
cfg_stat.numrandomization = nperms;

Nsub = length(allsubj);
cfg_stat.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg_stat.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg_stat.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg_stat.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat1 = ft_freqstatistics(cfg_stat,all_trl_avg{:},all_trl_null{:});
stat2 = ft_freqstatistics(cfg_stat,all_trl_Tposterior{:},all_trl_null{:});

%%%%%%%%%%%%%%%%%%%
%%% Sample-wise %%%
%%%%%%%%%%%%%%%%%%%
sess_r = [2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2];  % # session regressors per subject
pred_r = [2 3 3 3 3 3 3 3 3 3 2];  % # other predictors per sample postion

fprintf('Loading sample-wise data...\nSubject ')

samps = 2:12;  % samples to include for pupil*LLR analyses
DVpure_nLL_sum=[]; EVpure_nLL_sum=[]; DVpure_nLL_sum1=[]; EVpure_nLL_sum1=[]; Lpure_nLL_sum1=[];
for s = 1:length(allsubj)
    fprintf('%d, ',s)
    % Load data
    load([loadpath,allsubj{s},'_samplewise_output_appML',surpstr,bstr,'.mat'])
    freqsL = freqs;
    
    smp_avgL = smp_tf_avg;
    TpriorS_MRL = TpriorS_tf; TllrS_MRL = TllrS_tf; TllrXsurpriseS_MRL = TllrXsurpriseS_tf; TllrXuncertS_MRL = TllrXuncertS_tf; TllrXpupilS_MRL = TllrXpupilS_tf;
    
    for smp = 1:size(Rsq_DVpure_tf_BIC,3)
        DVpure_nLL_L(:,:,smp) = (Rsq_DVpure_tf_BIC(:,:,smp)-((1+pred_r(smp)+sess_r(s))*log(size(LLR_full,1))))./2;  % retrieve -LL for selected regression models, used later for superBIC
        EVpure_nLL_L(:,:,smp) = (Rsq_evidencepure_tf_BIC(:,:,smp)-((1+pred_r(smp)+sess_r(s))*log(size(LLR_full,1))))./2;
        
        DVpure_nLL_L1(:,:,smp) = (Rsq_DVpure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;  % retrieve -LL for selected regression models, used later for superBIC
        EVpure_nLL_L1(:,:,smp) = (Rsq_evidencepure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;
        Lpure_nLL_L1(:,:,smp) = (Rsq_Lpure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;
    end
    Rsq_DVdiff_L = Rsq_DV_tf-Rsq_evidence_tf; Rsq_pure_diff_L = Rsq_DVpure_tf-Rsq_evidencepure_tf;
    BIC_pure_diff_L = Rsq_DVpure_tf_BIC-Rsq_evidencepure_tf_BIC; BIC_pure_diff_L1 = Rsq_DVpure_tf_BIC1-Rsq_evidencepure_tf_BIC1; BIC_pure_PsiLdiff_L1 = Rsq_DVpure_tf_BIC1-Rsq_Lpure_tf_BIC1;
        
    load([loadpath,allsubj{s},'_samplewise_output_appML_HiFreq',surpstr,bstr,'.mat'])
    allfreqs = [freqsL freqs];
    
    % Store data for grand-average, sensor-resolved TFR
    smp_avg(:,:,s) = [squeeze(mean(smp_avgL(:,:,samps-1),3)); squeeze(mean(smp_tf_avg(:,:,samps-1),3))];

    gaB_priorS_MR(:,:,s) = [squeeze(mean(TpriorS_MRL(:,:,samps-1),3)); squeeze(mean(TpriorS_tf(:,:,samps-1),3))];
    gaB_llrS_MR(:,:,s) = [squeeze(mean(TllrS_MRL(:,:,samps-1),3)); squeeze(mean(TllrS_tf(:,:,samps-1),3))];
    gaB_llrXsurpriseS_MR(:,:,s) = [squeeze(mean(TllrXsurpriseS_MRL(:,:,samps-1),3)); squeeze(mean(TllrXsurpriseS_tf(:,:,samps-1),3))];
    gaB_llrXuncertS_MR(:,:,s) = [squeeze(mean(TllrXuncertS_MRL(:,:,samps-1),3)); squeeze(mean(TllrXuncertS_tf(:,:,samps-1),3))];
    gaB_llrXpupilS_MR(:,:,s) = [squeeze(mean(TllrXpupilS_MRL(:,:,samps-1),3)); squeeze(mean(TllrXpupilS_tf(:,:,samps-1),3))];
    
    gaB_surpriseMuncertS_MR(:,:,s) = [squeeze(mean(TllrXsurpriseS_MRL(:,:,samps-1),3)); squeeze(mean(TllrXsurpriseS_tf(:,:,samps-1),3))]-[squeeze(mean(TllrXuncertS_MRL(:,:,samps-1),3)); squeeze(mean(TllrXuncertS_tf(:,:,samps-1),3))];
    
    % gaRsq_DVdiff(:,:,s) = [squeeze(mean(Rsq_DVdiff_L(:,:,samps-1),3)); squeeze(mean(Rsq_DV_tf(:,:,samps-1)-Rsq_evidence_tf(:,:,samps-1),3))];
    gaRsq_pure_diff(:,:,s) = [squeeze(mean(Rsq_pure_diff_L(:,:,samps-1),3)); squeeze(mean(Rsq_DVpure_tf(:,:,samps-1)-Rsq_evidencepure_tf(:,:,samps-1),3))];
    gaRsq_pure_diff_BIC(:,:,s) = [squeeze(mean(BIC_pure_diff_L(:,:,samps-1),3)); squeeze(mean(Rsq_DVpure_tf_BIC(:,:,samps-1)-Rsq_evidencepure_tf_BIC(:,:,samps-1),3))];
    gaRsq_pure_diff_BIC1(:,:,s) = [squeeze(mean(BIC_pure_diff_L1(:,:,samps-1),3)); squeeze(mean(Rsq_DVpure_tf_BIC1(:,:,samps-1)-Rsq_evidencepure_tf_BIC1(:,:,samps-1),3))];
    gaRsq_pure_PsiLdiff_BIC1(:,:,s) = [squeeze(mean(BIC_pure_PsiLdiff_L1(:,:,samps-1),3)); squeeze(mean(Rsq_DVpure_tf_BIC1(:,:,samps-1)-Rsq_Lpure_tf_BIC1(:,:,samps-1),3))];
    
    for smp = 1:size(Rsq_DVpure_tf_BIC,3)
        DVpure_nLL(:,:,smp) = (Rsq_DVpure_tf_BIC(:,:,smp)-((1+pred_r(smp)+sess_r(s))*log(size(LLR_full,1))))./2;  % retrieve -LL for selected regression models, used later for superBIC
        EVpure_nLL(:,:,smp) = (Rsq_evidencepure_tf_BIC(:,:,smp)-((1+pred_r(smp)+sess_r(s))*log(size(LLR_full,1))))./2;
        npreds(smp) = 1+pred_r(smp)+sess_r(s);
        
        DVpure_nLL1(:,:,smp) = (Rsq_DVpure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;  % retrieve -LL for selected regression models, used later for superBIC
        EVpure_nLL1(:,:,smp) = (Rsq_evidencepure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;
        Lpure_nLL1(:,:,smp) = (Rsq_Lpure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;
        npreds1(smp) = 1+1+sess_r(s);
    end
    
    nObs(s) = size(LLR_full,1);  % recording # observations for use in later super BIC
    nPreds(s) = mean(npreds(samps-1)); nPreds1(s) = mean(npreds1(samps-1));  % recording average # predictors for use in later super BIC
    
    m = [squeeze(mean(DVpure_nLL_L(:,:,samps-1),3)); squeeze(mean(DVpure_nLL(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    DVpure_nLL_sum = sum(cat(3,DVpure_nLL_sum,m),3); % sum averaged nLL across subjects
    m = [squeeze(mean(EVpure_nLL_L(:,:,samps-1),3)); squeeze(mean(EVpure_nLL(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    EVpure_nLL_sum = sum(cat(3,EVpure_nLL_sum,m),3);
    
    m = [squeeze(mean(DVpure_nLL_L1(:,:,samps-1),3)); squeeze(mean(DVpure_nLL1(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    DVpure_nLL_sum1 = sum(cat(3,DVpure_nLL_sum1,m),3); % sum averaged nLL across subjects
    m = [squeeze(mean(EVpure_nLL_L1(:,:,samps-1),3)); squeeze(mean(EVpure_nLL1(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    EVpure_nLL_sum1 = sum(cat(3,EVpure_nLL_sum1,m),3);
    m = [squeeze(mean(Lpure_nLL_L1(:,:,samps-1),3)); squeeze(mean(Lpure_nLL1(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    Lpure_nLL_sum1 = sum(cat(3,Lpure_nLL_sum1,m),3);
    
    % Create structure for this subject to mimic output from FT_FREQANALYSIS
    cstruct = struct;
    cstruct.label = {'LI'};  % arbitrary label for motor lateralization 'channel'
    cstruct.fsample = 20;  % sampling rate
    cstruct.freq = allfreqs;
    cstruct.time = smptimes;
    cstruct.dimord = 'chan_freq_time';
    
    % Add structures to group-level arrays
    cstruct.powspctrm = reshape(squeeze(smp_avg(:,:,s)),1,size(smp_avg,1),size(smp_avg,2)); all_smp_avg{s} = cstruct;
    
    cstruct.powspctrm = reshape(squeeze(gaB_priorS_MR(:,:,s)),1,size(gaB_priorS_MR,1),size(gaB_priorS_MR,2)); allB_priorS_MR{s} = cstruct;
    cstruct.powspctrm = reshape(squeeze(gaB_llrS_MR(:,:,s)),1,size(gaB_llrS_MR,1),size(gaB_llrS_MR,2)); allB_llrS_MR{s} = cstruct;
    cstruct.powspctrm = reshape(squeeze(gaB_llrXsurpriseS_MR(:,:,s)),1,size(gaB_llrXsurpriseS_MR,1),size(gaB_llrXsurpriseS_MR,2)); allB_llrXsurpriseS_MR{s} = cstruct;
    cstruct.powspctrm = reshape(squeeze(gaB_llrXuncertS_MR(:,:,s)),1,size(gaB_llrXuncertS_MR,1),size(gaB_llrXuncertS_MR,2)); allB_llrXuncertS_MR{s} = cstruct;
    cstruct.powspctrm = reshape(squeeze(gaB_llrXpupilS_MR(:,:,s)),1,size(gaB_llrXpupilS_MR,1),size(gaB_llrXpupilS_MR,2)); allB_llrXpupilS_MR{s} = cstruct;
    
    cstruct.powspctrm = reshape(squeeze(gaB_surpriseMuncertS_MR(:,:,s)),1,size(gaB_surpriseMuncertS_MR,1),size(gaB_surpriseMuncertS_MR,2)); allB_surpriseMuncertS_MR{s} = cstruct;
    
    % cstruct.powspctrm = reshape(squeeze(gaRsq_DVdiff(:,:,s)),1,size(gaRsq_DVdiff,1),size(gaRsq_DVdiff,2)); allRsq_DVdiff{s} = cstruct;
    cstruct.powspctrm = reshape(squeeze(gaRsq_pure_diff(:,:,s)),1,size(gaRsq_pure_diff,1),size(gaRsq_pure_diff,2)); allRsq_pure_diff{s} = cstruct;
    cstruct.powspctrm = reshape(squeeze(gaRsq_pure_diff_BIC(:,:,s)),1,size(gaRsq_pure_diff_BIC,1),size(gaRsq_pure_diff_BIC,2)); allBIC_pure_diff{s} = cstruct;
    cstruct.powspctrm = reshape(squeeze(gaRsq_pure_diff_BIC1(:,:,s)),1,size(gaRsq_pure_diff_BIC1,1),size(gaRsq_pure_diff_BIC1,2)); allBIC_pure_diff1{s} = cstruct;
    cstruct.powspctrm = reshape(squeeze(gaRsq_pure_PsiLdiff_BIC1(:,:,s)),1,size(gaRsq_pure_PsiLdiff_BIC1,1),size(gaRsq_pure_PsiLdiff_BIC1,2)); allBIC_pure_PsiLdiff1{s} = cstruct;
    
    cstruct.powspctrm = zeros(size(cstruct.powspctrm)); allBnull{s} = cstruct;
end
fprintf('done.\nRunning stats.\n')

% Computing super BIC difference
ga_superBIC_pure_diff = (2*DVpure_nLL_sum)+(sum(nPreds).*log(sum(nObs))) - (2*DVpure_nLL_sum)+(sum(nPreds).*log(sum(nObs)));
ga_superBIC_pure_diff1 = (2*DVpure_nLL_sum1)+(sum(nPreds1).*log(sum(nObs))) - (2*DVpure_nLL_sum1)+(sum(nPreds1).*log(sum(nObs)));
ga_superBIC_pure_PsiLdiff1 = (2*DVpure_nLL_sum1)+(sum(nPreds1).*log(sum(nObs))) - (2*Lpure_nLL_sum1)+(sum(nPreds1).*log(sum(nObs)));

% Run cluster-based permutation test
stat11 = ft_freqstatistics(cfg_stat,all_smp_avg{:},allBnull{:});

stat12 = ft_freqstatistics(cfg_stat,allB_priorS_MR{:},allBnull{:});
stat13 = ft_freqstatistics(cfg_stat,allB_llrS_MR{:},allBnull{:});
stat14 = ft_freqstatistics(cfg_stat,allB_llrXsurpriseS_MR{:},allBnull{:});
stat15 = ft_freqstatistics(cfg_stat,allB_llrXuncertS_MR{:},allBnull{:});
stat16 = ft_freqstatistics(cfg_stat,allB_llrXpupilS_MR{:},allBnull{:});

stat17 = ft_freqstatistics(cfg_stat,allB_surpriseMuncertS_MR{:},allBnull{:});

% stat17 = ft_freqstatistics(cfg_stat,allRsq_DVdiff{:},allBnull{:});
stat18 = ft_freqstatistics(cfg_stat,allRsq_pure_diff{:},allBnull{:});
stat19 = ft_freqstatistics(cfg_stat,allBIC_pure_diff{:},allBnull{:});
stat20 = ft_freqstatistics(cfg_stat,allBIC_pure_diff1{:},allBnull{:});
stat21 = ft_freqstatistics(cfg_stat,allBIC_pure_PsiLdiff1{:},allBnull{:});


% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %
% --- %% BELIEF ENCODING FUNCTIONS %% --- %
% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %
% processing options
priorPLOT = 'LPR';  % switch between 'LPR' (x-axis of encoding profile plots will be posterior from previous sample) & 'psi' (x-axis will be Glaze-transformed prior for current sample)
nbins = 11;  % number of bins for plots
niter = 1000;  % number of iterations for bootsrapped CIs

fs = [8 14];

fprintf('Calculating encoding functions...\n')
allLI=[]; allLIprior=[]; allLIllr=[]; corrLPR=[];
for s = 1:length(allsubj)
    fprintf('Subj %d...\n',s)
    % load lo-freq data
    load([loadpath,allsubj{s},'_samplewise_output_appML',surpstr,bstr,'.mat'],'swLI','swLIlpr','swLIllr','freqs')
    swLIL=swLI; swLIlprL=swLIlpr; swLIllrL=swLIllr; freqsL=freqs;
    % load hi-freq data
    load([loadpath,allsubj{s},'_samplewise_output_appML_HiFreq',surpstr,bstr,'.mat'],'swLI','swLIlpr','swLIllr','LPR_full','psi_full','LLR_full','surprise_full','freqs')
    swLI = cat(2,swLIL,swLI);
    swLIlpr = cat(2,swLIlprL,swLIlpr);
    swLIllr = cat(2,swLIllrL,swLIllr);
    freqs = [freqsL freqs];
    % Average within frequency band of interest
    cLI = squeeze(mean(swLI(:,freqs>=fs(1) & freqs<=fs(2),:),2)).*-1;
    cLIprior = squeeze(mean(swLIlpr(:,freqs>=fs(1) & freqs<=fs(2),:),3)).*-1;
    cLIllr = squeeze(mean(swLIllr(:,freqs>=fs(1) & freqs<=fs(2),:),3)).*-1;
    % normalize computational variables if using model fits
    if ~isempty(strfind(modeltype,'fitted'))
        LPR_full = zscore(LPR_full,0,1);
        psi_full = zscore(psi_full,0,1);
    end
    % Selecting belief term by which to sort LI signal (LPR vs psi)
    if strcmp(priorPLOT,'psi')
        LPR_full(:,1:end-1) = psi_full(:,2:end);
    end
    % bin by computational variables & calculate correlations
    breaks=[];
    fprintf('Calculating actual encoding functions...\n')
    for smp = 2:12
        [allLI(1:nbins,1:2,smp-1,s,1),breaks(smp-1,:)] = bin4plot(LPR_full(:,smp),cLI(:,smp),nbins,0,[]);
        corrLPR(smp-1,s) = corr(allLI(:,1,smp-1,s,1),allLI(:,2,smp-1,s,1),'type','spearman');
    end
    % boostrap
    fprintf('Bootstrap iteration ')
    for i = 1:niter
        if mod(i,50)==0, fprintf('%d,',i), end
        for smp = 2:12
            c = randsample(size(cLI,1),size(cLI,1),true);
            allLI(1:nbins,1:2,smp-1,s,i+1) = bin4plot(LPR_full(c,smp),cLI(c,smp),nbins,0,breaks(smp-1,:));
        end
    end
    fprintf('done.\n')
end


% --- %%%%%%%%%% --- %
% --- %% SAVE %% --- %
% --- %%%%%%%%%% --- %
save([savepath,'Cluster_corr_TF',surpstr,bstr,'.mat'],'stat*','ga*','all_*','allfreqs','cfg_stat','trl_avg','trl_Tposterior','smp_avg',...
    'allLI','breaks','corrLPR','LPR_full','fs','priorPLOT','nbins','niter')

