function [] = runFreqAnalysis_appMotorLocalizer_tonly(n)
  % Applies two baselining operations (pre-trial and peri-sample) to MEG
  % TF data and regresses these data (channel*time) onto model-based
  % variables of interest

allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};
subject = allsubj{n};

basetype = 'dB_common';  % options: 'dB_st', 'dB_common', 'pc_st' or 'pc_common'

modeltype = 'fitted_linIU';  % switch b/w 'normative', 'fitted', 'fitted_np', 'fitted_lin', 'fitted_npIU' & 'fitted_linIU'
pupiltype = 'fixed';  % switch between 'fixed' (time-point picked from grand-av surprise encoding) and 'ss' (time-point picked from subject-specific surprise encoding)
priortype = 'psi';  % switch between 'LPR' (prior term will be posterior from previous sample) & 'psi' (prior term will be Glaze-transformed prior for current sample)
surprisetype = 'pCP';  % switch between 'pCP' and 'surprise'

% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
if strcmp(modeltype,'normative'), str2 = 'output';
elseif strcmp(modeltype,'fitted'), str2 = 'output_fitted';
elseif strcmp(modeltype,'fitted_np'), str2 = 'output_fitted_np';
elseif strcmp(modeltype,'fitted_lin'), str2 = 'output_fitted_lin';
elseif strcmp(modeltype,'fitted_npIU'), str2 = 'output_fitted_npIU';
elseif strcmp(modeltype,'fitted_linIU'), str2 = 'output_fitted_linIU'; end

loadpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/TF/wMotorLoc/',basetype,filesep,str2,filesep];

if strcmp(surprisetype,'pCP'), strs = '_pCP'; else strs=''; end

% ==================================================================
% TRIAL-WISE ANALYSIS
% ==================================================================
fprintf('Beginning trial-wise analysis, subject %d of %d...\n',n,length(allsubj))
load([loadpath,subject,'_samplewise_output_appML',strs,'.mat'],'LPR_full','LLR_full','surprise_full','psi_full','prior_full','dil_full','sess_full',...
    'choices_full','pswitch_full','fdist_full','distseq_full','win2');   % load sample-wise computational variables

load([loadpath,subject,'_trialwise_output_appML',strs,'.mat'],'trltimes','freqs','grad','cfg','t_ml');   % load trial-wise data - lo freq
t_mlL = t_ml;
load([loadpath,subject,'_trialwise_output_appML_HiFreq.mat'],'t_ml');   % load trial-wise data - hi freq
t_ml = squeeze(sum(cat(3,t_mlL,t_ml),3));  % sum across freqs (already weighted in previous processing step)

t_avg = mean([t_ml(choices_full==1,:); t_ml(choices_full==0,:).*-1],1);

% ==================================================================
% MAKE CATEGORICAL SESSION REGRESSORS
% ==================================================================
sessions = unique(sess_full);
sess_r = zeros(length(sess_full),length(sessions)-1);
for s = 1:length(sessions)-1
    sess_r(sess_full==sessions(s+1),s) = ones(length(find(sess_full==sessions(s+1))),1);
end

% ==================================================================
% RUN SINGLE-TRIAL REGRESSIONS
% ==================================================================
% Full-trial regressions, time*freq*space
Tposterior_t=[];
fprintf('Running regressions of full-trial motor preparation signal onto final posterior belief...')
for t = 1:size(t_ml,2)  % looping through time-points
    m = regstats(t_ml(:,t),[prior_full(:,end) sess_r],'linear',{'tstat'});  % signed posterior belief
    Tposterior_t(t) = m.tstat.t(2);
end
fprintf(' Done.\n')


% ==================================================================
% SAMPLE-WISE ANALYSIS
% ==================================================================
fprintf('Beginning sample-wise analysis...')
load([loadpath,subject,'_samplewise_output_appML',strs,'.mat'],'smptimes','t_ml');   % load trial-wise data - lo freq
t_mlL = t_ml;
load([loadpath,subject,'_samplewise_output_appML_HiFreq',strs,'.mat'],'t_ml');   % load trial-wise data - hi freq
t_ml = squeeze(sum(cat(4,t_mlL,t_ml),4));  % sum across freqs (already weighted in previous processing step)

smp_t_avg = squeeze(mean(cat(1,t_ml(choices_full==1,:,:),t_ml(choices_full==0,:,:).*-1),1));

% ==================================================================
% RUN SINGLE-TRIAL REGRESSIONS
% ==================================================================
% Sample-wise regressions, time
for s = 1:size(t_ml,3)-1  % looping through samples
    for t = 1:size(t_ml,2)  % looping through time-points
        m = regstats(t_ml(:,t,s+1),nanzscore([prior_full(:,s) LLR_full(:,s+1) LLR_full(:,s+1).*nanzscore(surprise_full(:,s+1)) LLR_full(:,s+1).*nanzscore(-abs(prior_full(:,s))) sess_r]),'linear',{'tstat','rsquare'});  % signed prior, LLR, LLR*surprise & LLR*uncertainty
        TpriorS_t(t,s) = m.tstat.t(2);
        TllrS_t(t,s) = m.tstat.t(3);
        TllrXsurpriseS_t(t,s) = m.tstat.t(4);
        TllrXuncertS_t(t,s) = m.tstat.t(5);
        Rsq_DV_t(t,s) = m.rsquare;
        
        m = regstats(t_ml(:,t,s+1),[LLR_full(:,s) LLR_full(:,s+1) LLR_full(:,s+1).*nanzscore(surprise_full(:,s+1)) LLR_full(:,s+1).*nanzscore(-abs(prior_full(:,s))) sess_r],'linear',{'rsquare'});  % previous LLR, current LLR and LLR*surprise
        Rsq_evidence_t(t,s) = m.rsquare;
        
        m = regstats(t_ml(:,t,s+1),[psi_full(:,max([2 s]):min([s+2,size(t_ml,3)])) sess_r],'linear',{'rsquare'});  % model including pure belief state after previous, current and next sample
        Rsq_DVpure_t(t,s) = m.rsquare;
        m = regstats(t_ml(:,t,s+1),[LLR_full(:,max([2 s]):min([s+2,size(t_ml,3)])) sess_r],'linear',{'rsquare'});  % model including pure evidence from previous, current and next sample
        Rsq_evidencepure_t(t,s) = m.rsquare;
        
        m = regstats(t_ml(:,t,s+1),[prior_full(:,s) LLR_full(:,s+1) LLR_full(:,s+1).*nanzscore(surprise_full(:,s+1)) LLR_full(:,s+1).*nanzscore(-abs(prior_full(:,s))) LLR_full(:,s+1).*nanzscore(dil_full(:,s+1)) sess_r],'linear',{'tstat'});  % signed prior, LLR and LLR*surprise
        TllrXpupilS_t(t,s) = m.tstat.t(6);
    end
end

fprintf('Done.\n')

% ==================================================================
% PULL PER-SAMPLE, TIME-AVERAGED SINGLE-TRIAL MOTOR SIGNAL, Z-SCORED WITHIN SESSION
% ==================================================================
win2 = [0.4 0.6];  % 2nd window relative to sample onset over which to average
swLI_t = nan(size(t_ml,1),size(t_ml,3));
swLIlpr_t = nan(size(t_ml,1),size(t_ml,3)-1);
swLIllr_t = nan(size(t_ml,1),size(t_ml,3)-1);
for s = 1:length(sessions)
    for smp = 1:size(t_ml,3)
        swLI_t(sess_full==sessions(s),smp) = squeeze(zscore(mean(t_ml(sess_full==sessions(s),smptimes>=win2(1) & smptimes<=win2(2),smp),2)));
    end
    for smp = 2:size(t_ml,3)
        m = regstats(squeeze(zscore(mean(t_ml(sess_full==sessions(s),smptimes>=win2(1) & smptimes<=win2(2),smp),2))),[LLR_full(sess_full==sessions(s),smp) LLR_full(sess_full==sessions(s),smp).*surprise_full(sess_full==sessions(s),smp)],'linear',{'r'});
        swLIlpr_t(sess_full==sessions(s),smp-1) = zscore(m.r);  % storing LI with current LLR & LLR*surprise regressed out (i.e. isolating a measure of pure prior encoding in LI)
        
        m = regstats(squeeze(zscore(mean(t_ml(sess_full==sessions(s),smptimes>=win2(1) & smptimes<=win2(2),smp),2))),[prior_full(sess_full==sessions(s),smp-1)],'linear',{'r'});
        swLIllr_t(sess_full==sessions(s),smp-1) = zscore(m.r);  % storing LI with prior belief regressed out (i.e. a measure of sample-wise *change* in LI)
    end
end

% ==================================================================
% SAVE RESULTS
% ==================================================================
save([loadpath,subject,'_output_appML',strs,'_tonly.mat'],'trltimes','t_avg','t_avg_CP','t_avg_CPnm','t_avg_CProc','t_avg_diffLLR','mEv','stEv','meanEv','n_matches','Tposterior_t','smptimes','smp_t_avg','grad','cfg',...
    'TpriorS_t','TllrS_t','TllrXsurpriseS_t','TllrXuncertS_t','TllrXpupilS_t','Rsq_DV_t','Rsq_evidence_t','Rsq_DVpure_t','Rsq_evidencepure_t',...
    'win2','swLI_t','swLIlpr_t','swLIllr_t','LPR_full','LLR_full','surprise_full','psi_full','prior_full','dil_full','sess_full','choices_full','pswitch_full','fdist_full','distseq_full')


