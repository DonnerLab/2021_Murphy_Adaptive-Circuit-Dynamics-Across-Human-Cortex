clear, close all

% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %
% --- %% CHOICE LATERALIZATION & SAMPLE-WISE REGRESSIONS %% --- %
% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %
% Specify inits
allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};

basetype = 'dB_common';  % options: 'dB_st', 'dB_common', 'pc_st' or 'pc_common'
modeltype = 'fitted_linIU';  % switch b/w 'normative', 'fitted', 'fitted_np' & 'fitted_lin'
surprisetype = 'pCP';  % switch between 'pCP' and 'surprise'

addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
% addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20171231'
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

if strcmp(surprisetype,'pCP'), strs='_pCP'; else strs=''; end

clustalpha = 0.05;  % % intial threhsold for forming clusters
nperms = 10;

%%%%%%%%%%%%%%%%%%
%%% Trial-wise %%%
%%%%%%%%%%%%%%%%%%
fprintf('Loading data...\nSubject ')
all_mEv=[]; all_stEv=[];
for s = 1:length(allsubj)
    fprintf('%d, ',s)
    % Load data
    load([loadpath,allsubj{s},'_output_appML',strs,'_tonly.mat'],'t_avg','t_avg_CP','t_avg_CPnm','t_avg_CProc','t_avg_diffLLR','mEv','stEv','meanEv','n_matches','Tposterior_t','trltimes')
        
    % Store data for grand-average, sensor-resolved TFR
    trl_avg(:,s) = t_avg;
    trl_avg_CP(:,s) = t_avg_CP;
    trl_avg_CPnm(:,s) = t_avg_CPnm;
    trl_avg_CProc(:,s) = t_avg_CProc;
    trl_Tposterior(:,s) = Tposterior_t;
    
    trl_avg_diffLLR(:,s) = t_avg_diffLLR;
    
    all_mEv = [all_mEv; mEv];
    all_stEv = [all_stEv; stEv'];
    ga_meanEv(:,:,s) = meanEv;
    ga_nmatches(:,s) = n_matches';

    % Create structure for this subject to mimic output from FT_FREQANALYSIS
    cstruct = struct;
    cstruct.label = {'LI'};  % arbitrary label for motor lateralization 'channel'
    cstruct.fsample = 20;  % sampling rate
    cstruct.time = trltimes;
    cstruct.dimord = 'chan_time';
    
    % Add structures to group-level arrays
    cstruct.powspctrm = t_avg; all_trl_avg{s} = cstruct;
    cstruct.powspctrm = t_avg_CP; all_trl_avg_CP{s} = cstruct;
    cstruct.powspctrm = t_avg_CPnm; all_trl_avg_CPnm{s} = cstruct;
    cstruct.powspctrm = t_avg_CProc-0.5; all_trl_avg_CProc{s} = cstruct;
    cstruct.powspctrm = Tposterior_t; all_trl_Tposterior{s} = cstruct;
    cstruct.powspctrm = zeros(size(cstruct.powspctrm));
    all_trl_null{s} = cstruct;
end
fprintf('done.\nRunning stats.\n')

% Run cluster-based permutation tests
cfg_stat = [];
cfg_stat.latency     = 'all';
cfg_stat.avgoverchan = 'no';
cfg_stat.avgovertime = 'no';
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

stat1 = ft_timelockstatistics(cfg_stat,all_trl_avg{:},all_trl_null{:});
stat2 = ft_timelockstatistics(cfg_stat,all_trl_Tposterior{:},all_trl_null{:});
stat3 = ft_timelockstatistics(cfg_stat,all_trl_avg_CP{:},all_trl_null{:});
% stat4 = ft_timelockstatistics(cfg_stat,all_trl_avg_CPnm{:},all_trl_null{:});
stat5 = ft_timelockstatistics(cfg_stat,all_trl_avg_CProc{:},all_trl_null{:});

figw = 26; figh = 9;
h = findobj('type','figure'); cfig = length(h)+1;

f = figure; set(cfig,'units','centimeters','pos',[0 0 figw figh],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
onsets = [0.4:0.4:(0.4*12)];
xlims = [-0.3 5.8]; ylims = [-0.68 0.17];
MLcol = [0 0 1]; lw=1; sigw=2.75; sigoffset=0.05; fs=10; fs_label=12;

s1=subplot(1,3,1); hold on
shadedErrorBar(stat1.time,mean(trl_avg.*-1,2),std(trl_avg.*-1,[],2)./sqrt(size(trl_avg,2)),{'color',MLcol,'LineWidth',lw},0)
shadedErrorBar(stat3.time,mean(trl_avg_CP.*-1,2),std(trl_avg_CP.*-1,[],2)./sqrt(size(trl_avg_CP,2)),{'color',[0 0 0],'LineWidth',lw},0)
% shadedErrorBar(stat3.time,mean(trl_avg_CPnm.*-1,2),std(trl_avg_CPnm.*-1,[],2)./sqrt(size(trl_avg_CPnm,2)),{'color',[0.5 0.5 0],'LineWidth',lw},0)
sig = double(stat1.mask); sig(sig==0) = nan; plot(stat1.time,sig.*(ylims(1)+diff(ylims)*sigoffset),'b','LineWidth',sigw)
sig = double(stat3.mask); sig(sig==0) = nan; plot(stat3.time,sig.*(ylims(1)+diff(ylims)*(sigoffset+0.02)),'k','LineWidth',sigw)
% sig = double(stat4.mask); sig(sig==0) = nan; plot(stat4.time,sig.*(ylims(1)+diff(ylims)*(sigoffset+0.04)),'Color',[0.5 0.5 0],'LineWidth',sigw)
line([0 0],ylims,'LineStyle','-','Color',[0 0 0],'LineWidth',1)
line(xlims,[0 0],'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',1)
line([onsets; onsets],[repmat(ylims',1,length(onsets))],'LineStyle','--','Color',[0.8 0.8 0.8]),
set(gca,'FontName','Arial','LineWidth',lw,'FontSize',fs,'TickDir','out','box','off')
xlim(xlims); ylim(ylims)
ylabel('Motor prep. (dB)','FontSize',fs_label), xlabel('Time from trial onset (s)','FontSize',fs_label)

s2=subplot(1,3,2); hold on
binranges = -4:1/3:4;
allbincounts = histc(all_mEv,binranges);
mbincounts = histc(all_stEv,binranges);
b1 = bar(binranges,allbincounts,'histc'); set(b1,'FaceColor','b','EdgeColor','b')
b2 = bar(binranges,mbincounts,'histc'); set(b2,'FaceColor','k','EdgeColor','k')
set(gca,'FontName','Arial','LineWidth',lw,'FontSize',fs-1,'TickDir','out','box','off')
ylabel('Trial count','FontSize',fs_label-2), xlabel('\mu(LLR)','FontSize',fs_label-2), xlim([-3.5 3.5])

ylim_roc = [0.45 0.66];
s3=subplot(1,3,3); hold on
shadedErrorBar(stat5.time,mean(trl_avg_CProc,2),std(trl_avg_CProc,[],2)./sqrt(size(trl_avg_CProc,2)),{'color',[0 0 0],'LineWidth',lw},0)
sig = double(stat5.mask); sig(sig==0) = nan; plot(stat5.time,sig.*(ylim_roc(1)+diff(ylim_roc)*(sigoffset+0.02)),'k','LineWidth',sigw)
line([0 0],ylim_roc,'LineStyle','-','Color',[0 0 0],'LineWidth',1)
line(xlims,[0.5 0.5],'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',1)
line([onsets; onsets],[repmat(ylim_roc',1,length(onsets))],'LineStyle','--','Color',[0.8 0.8 0.8]),
set(gca,'FontName','Arial','LineWidth',lw,'FontSize',fs,'TickDir','out','box','off')
xlim(xlims); ylim(ylim_roc)
ylabel('Area under ROC','FontSize',fs_label), xlabel('Time from trial onset (s)','FontSize',fs_label)

set(s1, 'Position', [0.12, 0.19, 0.29, 0.7])
set(s2, 'Position', [0.195, 0.37, 0.11, 0.17])
set(s3, 'Position', [0.55, 0.19, 0.29, 0.7])


%%%%%%%%%%%%%%%%%%%
%%% Sample-wise %%%
%%%%%%%%%%%%%%%%%%%
fprintf('Loading sample-wise data...\nSubject ')

samps = 2:12;  % samples to include for pupil*LLR analyses

for s = 1:length(allsubj)
    fprintf('%d, ',s)
    % Load data
    load([loadpath,allsubj{s},'_output_appML',strs,'_tonly.mat'])
    
    % Store data for grand-average, sensor-resolved TFR
    smp_avg(:,s) = mean(smp_t_avg(:,samps-1),2);
    
    gaB_priorS_MR(:,s) = mean(TpriorS_t(:,samps-1),2);
    gaB_llrS_MR(:,s) = mean(TllrS_t(:,samps-1),2);
    gaB_llrXsurpriseS_MR(:,s) = mean(TllrXsurpriseS_t(:,samps-1),2);
    gaB_llrXuncertS_MR(:,s) = mean(TllrXuncertS_t(:,samps-1),2);
    gaB_llrXpupilS_MR(:,s) = mean(TllrXpupilS_t(:,samps-1),2);

    gaRsq_DVdiff(:,s) = mean(Rsq_DV_t(:,samps-1)-Rsq_evidence_t(:,samps-1),2);
    gaRsq_pure_diff(:,s) = mean(Rsq_DVpure_t(:,samps-1)-Rsq_evidencepure_t(:,samps-1),2);
    
    % Create structure for this subject to mimic output from FT_FREQANALYSIS
    cstruct = struct;
    cstruct.label = {'LI'};  % arbitrary label for motor lateralization 'channel'
    cstruct.fsample = 20;  % sampling rate
    cstruct.time = smptimes;
    cstruct.dimord = 'chan_time';
    
    % Add structures to group-level arrays
    cstruct.powspctrm = smp_avg(:,s)'; all_smp_avg{s} = cstruct;
    
    cstruct.powspctrm = gaB_priorS_MR(:,s)'; allB_priorS_MR{s} = cstruct;
    cstruct.powspctrm = gaB_llrS_MR(:,s)'; allB_llrS_MR{s} = cstruct;
    cstruct.powspctrm = gaB_llrXsurpriseS_MR(:,s)'; allB_llrXsurpriseS_MR{s} = cstruct;
    cstruct.powspctrm = gaB_llrXuncertS_MR(:,s)'; allB_llrXuncertS_MR{s} = cstruct;
    cstruct.powspctrm = gaB_llrXpupilS_MR(:,s)'; allB_llrXpupilS_MR{s} = cstruct;
    
    cstruct.powspctrm = gaRsq_DVdiff(:,s)'; allRsq_DVdiff{s} = cstruct;
    cstruct.powspctrm = gaRsq_pure_diff(:,s)'; allRsq_pure_diff{s} = cstruct;
    
    cstruct.powspctrm = zeros(size(cstruct.powspctrm)); allBnull{s} = cstruct;
end
fprintf('done.\nRunning stats.\n')

% Run cluster-based permutation test
stat11 = ft_timelockstatistics(cfg_stat,all_smp_avg{:},allBnull{:});

stat12 = ft_timelockstatistics(cfg_stat,allB_priorS_MR{:},allBnull{:});
stat13 = ft_timelockstatistics(cfg_stat,allB_llrS_MR{:},allBnull{:});
stat14 = ft_timelockstatistics(cfg_stat,allB_llrXsurpriseS_MR{:},allBnull{:});
stat15 = ft_timelockstatistics(cfg_stat,allB_llrXuncertS_MR{:},allBnull{:});
stat16 = ft_timelockstatistics(cfg_stat,allB_llrXpupilS_MR{:},allBnull{:});

stat17 = ft_timelockstatistics(cfg_stat,allRsq_DVdiff{:},allBnull{:});
stat18 = ft_timelockstatistics(cfg_stat,allRsq_pure_diff{:},allBnull{:});


% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %
% --- %% BELIEF ENCODING FUNCTIONS %% --- %
% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %
% processing options
priorPLOT = 'LPR';  % switch between 'LPR' (x-axis of encoding profile plots will be posterior from previous sample) & 'psi' (x-axis will be Glaze-transformed prior for current sample)
nbins = 11;  % number of bins for plots
niter = 1000;  % number of iterations for bootsrapped CIs

fprintf('Calculating encoding functions...\n')
allLI=[]; allLIprior=[]; allLIllr=[]; corrLPR=[];
for s = 1:length(allsubj)
    fprintf('Subj %d...\n',s)
    % load lo-freq data
    load([loadpath,allsubj{s},'_output_appML',strs,'_tonly.mat'],'swLI_t','swLIlpr_t','swLIllr_t','LPR_full','psi_full')
    % Average within frequency band of interest
    cLI = swLI_t.*-1;
    cLIprior = swLIlpr_t.*-1;
    cLIllr = swLIllr_t.*-1;
    % normalize computational variables so can plot x-axis both normalized & original
    zLPR_full = zscore(LPR_full,0,1);
    zpsi_full = zscore(psi_full,0,1);
    % Selecting belief term by which to sort LI signal (LPR vs psi)
    if strcmp(priorPLOT,'psi')
        LPR_full(:,1:end-1) = psi_full(:,2:end);
        zLPR_full(:,1:end-1) = zpsi_full(:,2:end);
    end
    % bin by computational variables & calculate correlations
    breaks=[];
    fprintf('Calculating actual encoding functions...\n')
    for smp = 2:12
        [allLI(1:nbins,1:2,smp-1,s,1),breaks(smp-1,:)] = bin4plot(LPR_full(:,smp),cLI(:,smp),nbins,0,[]);
        [allLIz(1:nbins,1:2,smp-1,s,1),breaksz(smp-1,:)] = bin4plot(zLPR_full(:,smp),cLI(:,smp),nbins,0,[]);
        corrLPR(smp-1,s) = corr(allLI(:,1,smp-1,s,1),allLI(:,2,smp-1,s,1),'type','spearman');
    end
    % boostrap
    fprintf('Bootstrap iteration ')
    for i = 1:niter
        if mod(i,50)==0, fprintf('%d,',i), end
        for smp = 2:12
            c = randsample(size(cLI,1),size(cLI,1),true);
            allLI(1:nbins,1:2,smp-1,s,i+1) = bin4plot(LPR_full(c,smp),cLI(c,smp),nbins,0,breaks(smp-1,:));
            allLIz(1:nbins,1:2,smp-1,s,i+1) = bin4plot(zLPR_full(c,smp),cLI(c,smp),nbins,0,breaksz(smp-1,:));
        end
    end
    fprintf('done.\n')
end


% --- %%%%%%%%%% --- %
% --- %% SAVE %% --- %
% --- %%%%%%%%%% --- %
save([savepath,'Cluster_corr_TF',strs,'_tonly.mat'],'stat*','ga*','all_*','cfg_stat','trl_avg','trl_Tposterior','ga_nmatches','smp_avg',...
    'allLI','breaks','allLIz','breaksz','corrLPR','LPR_full','priorPLOT','nbins','niter')
    

