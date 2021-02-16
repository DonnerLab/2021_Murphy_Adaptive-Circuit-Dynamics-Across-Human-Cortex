function [] = runFreqAnalysis_appMotorLocalizer(n)
  % Applies two baselining operations (pre-trial and peri-sample) to MEG
  % TF data and regresses these data (channel*time) onto model-based
  % variables of interest

allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};
subject = allsubj{n};
  
basewin = [-0.4 -0.2];  % baseline window relative to pre-mask onset (s)
basetype = 'dB_common';  % options: 'dB_st', 'dB_common', 'pc_st' or 'pc_common'

smpwin = [0 1.2];  % window for sample-wise analyses
trlwin = [-0.5 5.8];  % window for full-trial-wise analyses

modeltype = 'fitted_linIU';  % switch b/w 'normative', 'fitted', 'fitted_np', 'fitted_lin', 'fitted_npIU' & 'fitted_linIU'
pupiltype = 'fixed';  % switch between 'fixed' (time-point picked from grand-av surprise encoding) and 'ss' (time-point picked from subject-specific surprise encoding)
priortype = 'psi';  % switch between 'LPR' (prior term will be posterior from previous sample) & 'psi' (prior term will be Glaze-transformed prior for current sample)
surprisetype = 'pCP';  % switch between 'pCP' and 'surprise'

coeftype = 'beta';   % switch b/w 'beta' and 'tscore'

% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun'))
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
if strcmp(modeltype,'normative'), str2 = 'output';
elseif strcmp(modeltype,'fitted'), str2 = 'output_fitted';
elseif strcmp(modeltype,'fitted_np'), str2 = 'output_fitted_np';
elseif strcmp(modeltype,'fitted_lin'), str2 = 'output_fitted_lin';
elseif strcmp(modeltype,'fitted_npIU'), str2 = 'output_fitted_npIU';
elseif strcmp(modeltype,'fitted_linIU'), str2 = 'output_fitted_linIU'; end

megpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/TF/wMotorLoc/';
if ~exist([savepath,basetype,filesep,str2,filesep],'dir'), mkdir([savepath,basetype,filesep,str2,filesep]), end  % making save directory if doesn't exist
savepath = [savepath,basetype,filesep,str2,filesep];

locpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/TF_Motor/',basetype,'/'];  % path to motor localizer weights

modelpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Modelling/'];
if strcmp(modeltype,'fitted'), load([modelpath,'Glaze_npLLR',filesep,'Fits',filesep,subject,'_fixed_fit.mat']);
elseif strcmp(modeltype,'fitted_np'), load([modelpath,'FitPhi_data_npLLR2',filesep,'Fits',filesep,subject,'_fixed_fit.mat']);
elseif strcmp(modeltype,'fitted_lin'), load([modelpath,'Glaze_basic',filesep,'Fits',filesep,subject,'_fixed_fit.mat']);
elseif strcmp(modeltype,'fitted_npIU'), load([modelpath,'Glaze_npLLR_InconUp',filesep,'Fits',filesep,subject,'_fixed_fit.mat']);
elseif strcmp(modeltype,'fitted_linIU'), load([modelpath,'Glaze_basic_InconUp',filesep,'Fits',filesep,subject,'_fixed_fit.mat']); end

subjfiles = dir([megpath,subject,'-*TF.mat']);  % pull all meg filenames for this subject

smp_data=[]; trl_data=[];
LLR_full=[]; LPR_full=[]; surprise_full=[]; psi_full=[];
dil_full=[]; X_full=[]; Y_full=[]; choices_full=[]; sess_full=[];
pswitch_full=[]; fdist_full=[]; distseq_full=[]; stimIn_full=[]; LLRinO_full=[];

% ==================================================================
% LOAD MOTOR LOCALIZER
% ==================================================================
load([locpath,'Motor_loc_weights.mat'])

% ==================================================================
% TRIAL-WISE ANALYSIS
% ==================================================================
fprintf('Beginning trial-wise analysis...\n')
for f = 1:length(subjfiles)

    load([megpath,subjfiles(f).name])  % load meg data
    sess = str2double(subjfiles(f).name(5));
    freq.time = round(freq.time,2);  % rounding time vector to nearest 2nd decimal - otherwise slight inaccuracies can lead to bad timing later
    freq.freq = round(freq.freq);    % same here
    
    % ==================================================================
    % PULL LEFT/RIGHT CHANNEL INDICES
    % ==================================================================
    if f==1
        Rchans=[]; Lchans=[];
        for c = find(strncmp(freq.label,'MR',2))'  % loop through each right channel and find matching left channel if it exists
            if ~isempty(find(strcmp(freq.label,['ML',freq.label{c}(3:end)])))
                Rchans(end+1,1) = c;
                Lchans(end+1,1) = find(strcmp(freq.label,['ML',freq.label{c}(3:end)]));
            end
        end
    end
    
    % ==================================================================
    % PULL DESIRED SEGMENTS OF DATA
    % ==================================================================
    fprintf('Concatenating trial-wise data segments for %s, lo-freq...\n',subjfiles(f).name)
    if f==1
        freqs = freq.freq;
        trltimes = freq.time(freq.time>=trlwin(1) & freq.time<=trlwin(2));
    end
    trl_data = cat(1,trl_data,freq.powspctrm(:,:,freqs>=min(allfreqs) & freqs<=max(allfreqs),freq.time>=trlwin(1) & freq.time<=trlwin(2)));
    
    freqs = freqs(freqs>=min(allfreqs) & freqs<=max(allfreqs));  % trimming freqs to align with motor localizer if necessary (since localizer was fit to pooled freqs)
    if f==1
        w_freq = w_freq(ismember(allfreqs,freqs));  % trimming weight vectors/matrices to align freqs with current dataset (since localizer was fit to pooled freqs)
        w_sens_freq = w_sens_freq(:,ismember(allfreqs,freqs));
    end
    
    % ==================================================================
    % CONCATENATE MODEL-BASED VARIABLES
    % ==================================================================
    % concatenate behaviour/pupil
    stimIn_full = [stimIn_full; freq.mdlvars.stimIn];
    choices_full = [choices_full; freq.mdlvars.choices];
    sess_full = [sess_full; ones(length(freq.mdlvars.choices),1).*sess];
    pswitch_full = [pswitch_full; freq.mdlvars.pswitch];
    fdist_full = [fdist_full; freq.mdlvars.fdist];
    distseq_full = [distseq_full; freq.mdlvars.distseq];
    LLRinO_full = [LLRinO_full; freq.mdlvars.LLR];
    
    if strcmp(modeltype,'normative')
        LLR_full = [LLR_full; freq.mdlvars.LLR];
        LPR_full = [LPR_full; freq.mdlvars.LPR];
        psi_full = [psi_full; freq.mdlvars.psi];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; freq.mdlvars.surpriseW];    % PICK SURPRISE TYPE HERE
        end
    elseif strcmp(modeltype,'fitted')
        LLR_full = [LLR_full; freq.mdlvars.LLR_M];
        LPR_full = [LPR_full; freq.mdlvars.LPR_M];
        psi_full = [psi_full; freq.mdlvars.psi_M];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; freq.mdlvars.surpriseW_M];     % PICK SURPRISE TYPE HERE
        end
    elseif strcmp(modeltype,'fitted_np')
        LLR_full = [LLR_full; freq.mdlvars.LLR_Mnp];
        LPR_full = [LPR_full; freq.mdlvars.LPR_Mnp];
        psi_full = [psi_full; freq.mdlvars.psi_Mnp];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; freq.mdlvars.surpriseW_Mnp];     % PICK SURPRISE TYPE HERE
        end
    elseif strcmp(modeltype,'fitted_lin')
        LLR_full = [LLR_full; freq.mdlvars.LLR_Mlin];
        LPR_full = [LPR_full; freq.mdlvars.LPR_Mlin];
        psi_full = [psi_full; freq.mdlvars.psi_Mlin];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; freq.mdlvars.surpriseW_Mlin];     % PICK SURPRISE TYPE HERE
        end
    end
    
    if strcmp(pupiltype,'ss'), dil_full = [dil_full; freq.pupilvars.dil]; elseif strcmp(pupiltype,'fixed'), dil_full = [dil_full; freq.pupilvars.dilav]; end
    if strcmp(pupiltype,'ss'), X_full = [X_full; freq.pupilvars.X]; elseif strcmp(pupiltype,'fixed'), X_full = [X_full; freq.pupilvars.Xav]; end
    if strcmp(pupiltype,'ss'), Y_full = [Y_full; freq.pupilvars.Y]; elseif strcmp(pupiltype,'fixed'), Y_full = [Y_full; freq.pupilvars.Yav]; end
    
    % ==================================================================
    % STORE FT STRUCTURES FOR LATER PLOTTING
    % ==================================================================
    if f==1,
        cfg = freq.cfg;
    end
end
clear freq

% ================================================
% COMPUTE CHANGE-POINT PROBABILITY IF REQUIRED
% ================================================
if strcmp(surprisetype,'pCP')
    if strcmp(modeltype,'normative')
        pIn = cat(3,normpdf(stimIn_full,17,29),normpdf(stimIn_full,-17,29));
        H = 0.08;
    elseif strcmp(modeltype,'fitted')
        gains = cumsum(pm_fit(3:end)');
        gains = [fliplr(gains) 0 gains];  % creating full vector of gain terms, symmetric around LLR=0
        rLLRfull_norm = [sort(-rLLR) 0 rLLR];  % creating full vector of matching original LLRs
        
        LLRext = log(normpdf(-100:1:100,17,29)./normpdf(-100:1:100,-17,29)); % Fit p(stim) from fitted LLRs under constraint that integral of p(stim) within stimulus bounds = 1
        xy = [interp1(LLRext,-100:1:100,rLLRfull_norm,'linear')' (rLLRfull_norm.*gains)'];  % [x=true stimulus positions, y=fitted LLRs]
        log_pstim = fminsearchbnd(@(params) fit_pstim_from_LLR(params,xy), log(normpdf(xy(:,1),gen.mu(1),gen.sigma(1))), ones(1,size(xy,1)).*-60, ones(1,size(xy,1)).*-1);  % fit p(stim) in log space
        stim=xy(:,1)'; pstim = exp(interp1(stim,log_pstim,[min(stim) -90:1:90 max(stim)],'spline')); pstim = pstim./sum(pstim);  % calculate normalized hi-res p(stim)
        
        for samp = 1:size(LLR_full,2)
            pIn(:,samp,1) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(pstim(subj,:)),stimIn_full(:,samp),'spline'));
            pIn(:,samp,2) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(fliplr(pstim)),stimIn_full(:,samp),'spline'));
        end
        H = pm_fit(1);
        [~,surprise,~] = accGlaze_fast(LLR_full,H,0,'pCP',pIn);
        surprise_full = [surprise_full; log(surprise)];
        
    elseif strcmp(modeltype,'fitted_np')
        gains = cumsum(pm_fit(length(L_nm1)+1:length(L_nm1)+length(rLLR))');
        gains = [fliplr(gains) 0 gains];  % creating full vector of gain terms, symmetric around LLR=0
        rLLRfull = [sort(-rLLR) 0 rLLR];  % creating full vector of matching original LLRs
        
        LLRext = log(normpdf(-100:1:100,17,29)./normpdf(-100:1:100,-17,29)); % Fit p(stim) from fitted LLRs under constraint that integral of p(stim) within stimulus bounds = 1
        xy = [interp1(LLRext,-100:1:100,rLLRfull,'linear')' (rLLRfull.*gains)'];  % [x=true stimulus positions, y=fitted LLRs]
        log_pstim = fminsearchbnd(@(params) fit_pstim_from_LLR(params,xy), log(normpdf(xy(:,1),gen.mu(1),gen.sigma(1))), ones(1,size(xy,1)).*-60, ones(1,size(xy,1)).*-1);  % fit p(stim) in log space
        stim=xy(:,1)'; pstim = exp(interp1(stim,log_pstim,[min(stim) -90:1:90 max(stim)],'spline')); pstim = pstim./sum(pstim);  % calculate normalized hi-res p(stim)
        
        for samp = 1:size(LLR_full,2)
            pIn(:,samp,1) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(pstim),stimIn_full(:,samp),'spline'));
            pIn(:,samp,2) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(fliplr(pstim)),stimIn_full(:,samp),'spline'));
        end
        H = 0.08;
        [~,surprise,~] = accGlaze_fast(LLR_full,H,0,'pCP',pIn);
        surprise_full = [surprise_full; log(surprise)];
        
    elseif strcmp(modeltype,'fitted_lin')  % for linear fits, can infer p(stim) exactly under assumption of symmetric gaussians
        cLLR = LLR_full(1,2); cStim = stimIn_full(1,2);
        sigma = sqrt((-(cStim-17).^2 + (cStim+17).^2)./(2.*cLLR));
        pIn = cat(3,normpdf(stimIn_full,17,sigma),normpdf(stimIn_full,-17,sigma));
        H = pm_fit(1);
        [~,surprise,~] = accGlaze_fast(LLR_full,H,0,'pCP',pIn);
        surprise_full = [surprise_full; log(surprise)];
        
    elseif strcmp(modeltype,'fitted_npIU')
        gains = cumsum(pm_fit(4:end)');
        gains = [fliplr(gains) 0 gains];  % creating full vector of gain terms, symmetric around LLR=0
        rLLRfull_norm = [sort(-rLLR) 0 rLLR];  % creating full vector of matching original LLRs
        
        LLRext = log(normpdf(-100:1:100,17,29)./normpdf(-100:1:100,-17,29)); % Fit p(stim) from fitted LLRs under constraint that integral of p(stim) within stimulus bounds = 1
        xy = [interp1(LLRext,-100:1:100,rLLRfull_norm,'linear')' (rLLRfull_norm.*gains)'];  % [x=true stimulus positions, y=fitted LLRs]
        log_pstim = fminsearchbnd(@(params) fit_pstim_from_LLR(params,xy), log(normpdf(xy(:,1),gen.mu(1),gen.sigma(1))), ones(1,size(xy,1)).*-60, ones(1,size(xy,1)).*-1);  % fit p(stim) in log space
        stim=xy(:,1)'; pstim = exp(interp1(stim,log_pstim,[min(stim) -90:1:90 max(stim)],'spline')); pstim = pstim./sum(pstim);  % calculate normalized hi-res p(stim)
        
        LLR_full = nan(size(LLRinO_full));
        for samp = 1:size(LLR_full,2)
            LLR_full(:,samp) = interp1(rLLRfull_norm,rLLRfull_norm.*gains,LLRinO_full(:,samp),'spline');
            pIn(:,samp,1) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(pstim),stimIn_full(:,samp),'spline'));
            pIn(:,samp,2) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(fliplr(pstim)),stimIn_full(:,samp),'spline'));
        end
        H = pm_fit(1);
        [LPR_full,surprise_full,psi_full] = accGlaze_InconUp_fast(LLR_full,H,pm_fit(3),0,'pCP',pIn);
        surprise_full = log(surprise_full);
        
    elseif strcmp(modeltype,'fitted_linIU')
        LLR_full = LLRinO_full.*pm_fit(2);
        cLLR = LLR_full(1,2); cStim = stimIn_full(1,2);
        sigma = sqrt((-(cStim-17).^2 + (cStim+17).^2)./(2.*cLLR));
        pIn = cat(3,normpdf(stimIn_full,17,sigma),normpdf(stimIn_full,-17,sigma));
        H = pm_fit(1);
        [LPR_full,surprise_full,psi_full] = accGlaze_InconUp_fast(LLR_full,H,pm_fit(4),0,'pCP',pIn);
        surprise_full = log(surprise_full);
    end
end

% ==================================================================
% ENSURE USE OF DESIRED FORM OF PRIOR BELIEF
% ==================================================================
prior_full = LPR_full;
if strcmp(priortype,'psi')
    prior_full(:,1:end-1) = psi_full(:,2:end);  % N.B. this way, final term in LPR_full will always be final, untransformed belief (i.e. the quantity that the model uses to make a decision)
end

% ==================================================================
% APPLY BASELINES
% ==================================================================
fprintf('Applying baselines...\n')
if strcmp(basetype,'dB_st')
    trl_data = 10.*log10(trl_data);  % convert to dB
    bl = nanmean(trl_data(:, :, :, trltimes>=basewin(1) & trltimes<=basewin(2)), 4);  % extract dB-transformed pre-onset baselines per trial, channel & freq
    for t = 1:size(bl,1)   % looping through trials
        for c = 1:size(bl,2)  % looping through channels
            trl_data(t,c,:,:) = squeeze(trl_data(t,c,:,:))-repmat(squeeze(bl(t,c,:)),1,length(trltimes));  % applying baselines to dB-transformed data
        end
    end
    
elseif strcmp(basetype,'dB_common')
    trl_data = 10.*log10(trl_data);  % convert to dB
    bl = squeeze(nanmean(nanmean(trl_data(:, :, :, trltimes>=basewin(1) & trltimes<=basewin(2)),4),1));  % extract dB-transformed pre-onset baselines per channel & freq
    for t = 1:size(trl_data,1)   % looping through trials
        for c = 1:size(bl,1)  % looping through channels
            trl_data(t,c,:,:) = squeeze(trl_data(t,c,:,:))-repmat(bl(c,:)',1,length(trltimes));  % applying baselines to dB-transformed data
        end
    end
    
elseif strcmp(basetype,'pc_st')
    bl = nanmean(trl_data(:, :, :, trltimes>=basewin(1) & trltimes<=basewin(2)), 4);  % extract pre-onset baselines per trial, channel & freq
    for t = 1:size(bl,1)   % looping through trials
        for c = 1:size(bl,2)  % looping through channels
            trl_data(t,c,:,:) = (squeeze(trl_data(t,c,:,:))-repmat(squeeze(bl(t,c,:)),1,length(trltimes)))./repmat(squeeze(bl(t,c,:)),1,length(trltimes)).*100; % applying %-change
        end
    end
    
elseif strcmp(basetype,'pc_common')
    bl = squeeze(nanmean(nanmean(trl_data(:, :, :, trltimes>=basewin(1) & trltimes<=basewin(2)),4),1));  % extract pre-onset baselines per channel & freq
    for t = 1:size(trl_data,1)   % looping through trials
        for c = 1:size(bl,1)  % looping through channels
            trl_data(t,c,:,:) = (squeeze(trl_data(t,c,:,:))-repmat(bl(c,:)',1,length(trltimes)))./repmat(bl(c,:)',1,length(trltimes)).*100; % applying %-change
        end
    end
end

% ==================================================================
% COMPUTE LATERALIZATION INDICES
% ==================================================================
for c = 1:length(Rchans)
    trl_data(:,c,:,:) = trl_data(:,Rchans(c),:,:) - trl_data(:,Lchans(c),:,:);  % store full-trial, half-scalp LI data - trials*chans*freqs*times (overwriting to save memory - works b/c L/R chans are well separated)
end
trl_data(:,length(Rchans)+1:end,:,:) = [];

% ==================================================================
% APPLY MOTOR LOCALIZER WEIGHTS
% ==================================================================
trl_datax = trl_data;  % create duplicate
% Apply only sensor-domain weights (leaving freq*time info)
for c = 1:length(w_sens)
    trl_datax(:,c,:,:) = trl_datax(:,c,:,:).*w_sens(c);
end
tf_ml = squeeze(sum(trl_datax,2));  % leaves trials*freqs*time matrix
clear trl_datax

% Apply joint sensor-/freq-domain weights (leaving only a time series)
t_ml=[];
for tr = 1:size(trl_data,1)
    for t = 1:size(trl_data,4)
        t_ml(tr,t) = sum(sum(squeeze(trl_data(tr,:,:,t)).*w_sens_freq));  % leaves trials*time matrix
    end
end
clear trl_data
    
tf_avg = squeeze(mean(cat(1,tf_ml(choices_full==1,:,:),tf_ml(choices_full==0,:,:).*-1),1));
% t_avg = mean([t_ml(choices_full==1,:); t_ml(choices_full==0,:).*-1],1);

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
TposteriorTF=[]; TposteriorT=[];
fprintf('Running regressions of full-trial power (time*freq) onto final posterior belief...\n Sensor')
for f = 1:size(tf_ml,2)  % looping through freqs
    for t = 1:size(tf_ml,3)  % looping through time-points
        m = regstats(nanzscore(tf_ml(:,f,t)),nanzscore([prior_full(:,end) sess_r]),'linear',{'tstat','beta'});  % signed posterior belief
        if strcmp(coeftype,'beta')
            TposteriorTF(f,t) = m.beta(2);
        elseif strcmp(coeftype,'tscore')
            TposteriorTF(f,t) = m.tstat.t(2);
        end
    end
end
% for t = 1:size(t_ml,2)  % looping through time-points
%     m = regstats(t_ml(:,t),[prior_full(:,end) sess_r],'linear',{'tstat'});  % signed posterior belief
%     TposteriorT(t) = m.tstat.t(2);
% end
fprintf(' Done.\n')

% ==================================================================
% SAVE RESULTS AND CLEAN UP
% ==================================================================
if strcmp(surprisetype,'pCP')
    if strcmp(coeftype,'beta')
        savename = [savepath,subject,'_trialwise_output_appML_pCP_beta.mat'];
    elseif strcmp(coeftype,'tscore')
        savename = [savepath,subject,'_trialwise_output_appML_pCP.mat'];
    end
else
    if strcmp(coeftype,'beta')
       savename = [savepath,subject,'_trialwise_output_appML_beta.mat'];
    elseif strcmp(coeftype,'tscore')
        savename = [savepath,subject,'_trialwise_output_appML.mat'];
    end
end
save(savename,'trltimes','freqs','grad','cfg','tf_avg','TposteriorTF','t_ml')

clear tf_avg t_avg tf_avg t_avg TposteriorTF TposteriorT tf_ml t_ml


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ==================================================================
% SAMPLE-WISE ANALYSIS
% ==================================================================
fprintf('Beginning sample-wise analysis...\n')
for f = 1:length(subjfiles)

    load([megpath,subjfiles(f).name])  % load meg data
    
    % aligning freqs with motor localizer
    freq.freq = round(freq.freq);
    freq.powspctrm = freq.powspctrm(:,:,freq.freq>=min(allfreqs) & freq.freq<=max(allfreqs),:);
    freq.freq = freq.freq(freq.freq>=min(allfreqs) & freq.freq<=max(allfreqs));  % trimming freqs to align with motor localizer if necessary (since localizer was fit to pooled freqs)
    
    % ==================================================================
    % APPLY BASELINES
    % ==================================================================
    fprintf('Applying baseline for %s...\n',subjfiles(f).name)
    
    freq.time = round(freq.time,2);  % rounding time vector to nearest 2nd decimal - otherwise slight inaccuracies can lead to bad timing later
    
    if strcmp(basetype,'dB_st')
        freq.powspctrm = 10.*log10(freq.powspctrm);  % convert to dB
        stbase = nanmean(freq.powspctrm(:, :, :, freq.time>=basewin(1) & freq.time<=basewin(2)), 4);  % extract dB-transformed pre-onset baselines per trial, channel & freq
        for t = 1:size(stbase,1)   % looping through trials
            for c = 1:size(stbase,2)  % looping through channels
                freq.powspctrm(t,c,:,:) = squeeze(freq.powspctrm(t,c,:,:))-repmat(squeeze(stbase(t,c,:)),1,length(freq.time));  % applying baselines to dB-transformed data
            end
        end
        
    elseif strcmp(basetype,'dB_common')   % uses trial-averaged baseline computed earlier
        freq.powspctrm = 10.*log10(freq.powspctrm);  % convert to dB
        for t = 1:size(freq.powspctrm,1)   % looping through trials
            for c = 1:size(bl,1)  % looping through channels
                freq.powspctrm(t,c,:,:) = squeeze(freq.powspctrm(t,c,:,:))-repmat(bl(c,:)',1,length(freq.time));  % applying baselines to dB-transformed data
            end
        end
        
    elseif strcmp(basetype,'pc_st')
        stbase = nanmean(freq.powspctrm(:, :, :, freq.time>=basewin(1) & freq.time<=basewin(2)), 4);  % extract  pre-onset baselines per trial, channel & freq
        for t = 1:size(stbase,1)   % looping through trials
            for c = 1:size(stbase,2)  % looping through channels
                freq.powspctrm(t,c,:,:) = (squeeze(freq.powspctrm(t,c,:,:))-repmat(squeeze(stbase(t,c,:)),1,length(freq.time)))./repmat(squeeze(stbase(t,c,:)),1,length(freq.time)).*100; % applying %-change
            end
        end
        
    elseif strcmp(basetype,'pc_common')   % uses trial-averaged baseline computed earlier
        for t = 1:size(freq.powspctrm,1)   % looping through trials
            for c = 1:size(bl,1)  % looping through channels
                freq.powspctrm(t,c,:,:) = (squeeze(freq.powspctrm(t,c,:,:))-repmat(bl(c,:)',1,length(freq.time)))./repmat(bl(c,:)',1,length(freq.time)).*100; % applying %-change
            end
        end
    end
    
    % ==================================================================
    % PULL DESIRED SEGMENTS OF DATA
    % ==================================================================
    fprintf('Concatenating trial- & sample-wise data segments...\n')
    onsets = 0.4:0.4:0.4*12;  % vector of all sample onset times relative to pre-mask
    if f==1,
        smptimes = freq.time(freq.time>=smpwin(1) & freq.time<=smpwin(2));  % getting vector of sample times relative to dot onset
    end
    
    for s = 1:length(onsets)
        stsmp = find(freq.time>=onsets(s),1,'first');
        if f==1 && s==1, newts = 1:size(freq.powspctrm,1); elseif f>1 && s==1, newts = size(smp_data,1)+1:size(smp_data,1)+size(freq.powspctrm,1); end
        smp_data(newts,:,:,:,s) = freq.powspctrm(:,:,:,stsmp:stsmp+length(smptimes)-1);
    end
end
clear freq

% ==================================================================
% COMPUTE LATERALIZATION INDICES
% ==================================================================
for c = 1:length(Rchans)
    smp_data(:,c,:,:,:) = smp_data(:,Rchans(c),:,:,:) - smp_data(:,Lchans(c),:,:,:);  % store full-trial, half-scalp LI data - trials*chans*freqs*times*samples (overwriting to save memory)
end
smp_data(:,length(Rchans)+1:end,:,:,:) = [];

% ==================================================================
% APPLY MOTOR LOCALIZER WEIGHTS
% ==================================================================
smp_datax = smp_data;  % create duplicate
% Apply only sensor-domain weights (leaving freq*time info)
for c = 1:length(w_sens)
    smp_datax(:,c,:,:,:) = smp_datax(:,c,:,:,:).*w_sens(c);
end
tf_ml = squeeze(sum(smp_datax,2));  % leaves trials*freqs*time*samples matrix
clear smp_datax

% Apply joint sensor-/freq-domain weights (leaving only a time series)
t_ml=[];
for smp = 1:size(smp_data,5)
    for tr = 1:size(smp_data,1)
        for t = 1:size(smp_data,4)
            t_ml(tr,t,smp) = sum(sum(squeeze(smp_data(tr,:,:,t,smp)).*w_sens_freq));  % leaves trials*time*samples matrix
        end
    end
end
clear smp_data

% ==================================================================
% COMPUTE TRIAL-AVERAGED RESPONSES
% ==================================================================
smp_tf_avg = squeeze(mean(cat(1,tf_ml(choices_full==1,:,:,:),tf_ml(choices_full==0,:,:,:).*-1),1));
% smp_t_avg = squeeze(mean(cat(1,t_ml(choices_full==1,:,:),t_ml(choices_full==0,:,:).*-1),1));

% ==================================================================
% RUN SINGLE-TRIAL REGRESSIONS
% ==================================================================
% Sample-wise regressions, time*freq
fprintf('Running regressions of sample-wise power (time*freq*space) onto certainty, surprise and pupil dilation...\n Sample ')
for s = 1:size(tf_ml,4)-1  % looping through samples
    fprintf('%d, ',s)
    for f = 1:size(tf_ml,2)  % looping through freqs
        for t = 1:size(tf_ml,3)  % looping through time-points
            m = regstats(nanzscore(tf_ml(:,f,t,s+1)),nanzscore([prior_full(:,s) LLR_full(:,s+1) LLR_full(:,s+1).*nanzscore(surprise_full(:,s+1)) LLR_full(:,s+1).*nanzscore(-abs(prior_full(:,s))) sess_r]),'linear',{'tstat','beta','rsquare'});  % signed prior, LLR, LLR*surprise & LLR*uncertainty
            if strcmp(coeftype,'beta')
                TpriorS_tf(f,t,s) = m.beta(2);
                TllrS_tf(f,t,s) = m.beta(3);
                TllrXsurpriseS_tf(f,t,s) = m.beta(4);
                TllrXuncertS_tf(f,t,s) = m.beta(5);
            elseif strcmp(coeftype,'tscore')
                TpriorS_tf(f,t,s) = m.tstat.t(2);
                TllrS_tf(f,t,s) = m.tstat.t(3);
                TllrXsurpriseS_tf(f,t,s) = m.tstat.t(4);
                TllrXuncertS_tf(f,t,s) = m.tstat.t(5);
            end
            Rsq_DV_tf(f,t,s) = m.rsquare;
            
            m = regstats(nanzscore(tf_ml(:,f,t,s+1)),nanzscore([LLR_full(:,s) LLR_full(:,s+1) LLR_full(:,s+1).*nanzscore(surprise_full(:,s+1)) LLR_full(:,s+1).*nanzscore(-abs(prior_full(:,s))) sess_r]),'linear',{'rsquare'});  % previous LLR, current LLR and LLR*surprise
            Rsq_evidence_tf(f,t,s) = m.rsquare;
            
            [~,~,stats] = glmfit([prior_full(:,max([2 s]):min([s+2,size(tf_ml,4)])) sess_r],tf_ml(:,f,t,s+1));
            Rsq_DVpure_tf(f,t,s) = 1-(sum(stats.resid.^2)./sum((tf_ml(:,f,t,s+1)-mean(tf_ml(:,f,t,s+1))).^2));
            Rsq_DVpure_tf_BIC(f,t,s) = 2*normlike([0 stats.s],stats.resid) + ((1+size([prior_full(:,max([2 s]):min([s+2,size(tf_ml,4)])) sess_r],2))*log(size(LLR_full,1)));  % BIC = 2L + klog(n), where L is negative log-likelihood
            [~,~,stats] = glmfit([LLR_full(:,max([2 s]):min([s+2,size(tf_ml,4)])) sess_r],tf_ml(:,f,t,s+1));
            Rsq_evidencepure_tf(f,t,s) = 1-(sum(stats.resid.^2)./sum((tf_ml(:,f,t,s+1)-mean(tf_ml(:,f,t,s+1))).^2));
            Rsq_evidencepure_tf_BIC(f,t,s) = 2*normlike([0 stats.s],stats.resid) + ((1+size([LLR_full(:,max([2 s]):min([s+2,size(tf_ml,4)])) sess_r],2))*log(size(LLR_full,1)));  % BIC = 2L + klog(n), where L is negative log-likelihood
            
            [~,~,stats] = glmfit([prior_full(:,s+1) sess_r],tf_ml(:,f,t,s+1));
            Rsq_DVpure_tf_BIC1(f,t,s) = 2*normlike([0 stats.s],stats.resid) + ((1+size([prior_full(:,s+1) sess_r],2))*log(size(LLR_full,1)));  % BIC = 2L + klog(n), where L is negative log-likelihood
            [~,~,stats] = glmfit([LLR_full(:,s+1) sess_r],tf_ml(:,f,t,s+1));
            Rsq_evidencepure_tf_BIC1(f,t,s) = 2*normlike([0 stats.s],stats.resid) + ((1+size([LLR_full(:,s+1) sess_r],2))*log(size(LLR_full,1)));  % BIC = 2L + klog(n), where L is negative log-likelihood
            [~,~,stats] = glmfit([LPR_full(:,s+1) sess_r],tf_ml(:,f,t,s+1));
            Rsq_Lpure_tf_BIC1(f,t,s) = 2*normlike([0 stats.s],stats.resid) + ((1+size([LPR_full(:,s+1) sess_r],2))*log(size(LLR_full,1)));  % BIC = 2L + klog(n), where L is negative log-likelihood
            
            m = regstats(nanzscore(tf_ml(:,f,t,s+1)),nanzscore([prior_full(:,s) LLR_full(:,s+1) LLR_full(:,s+1).*nanzscore(surprise_full(:,s+1)) LLR_full(:,s+1).*nanzscore(-abs(prior_full(:,s))) LLR_full(:,s+1).*nanzscore(dil_full(:,s+1)) sess_r]),'linear',{'tstat','beta'});  % signed prior, LLR and LLR*surprise
            if strcmp(coeftype,'beta')
                TllrXpupilS_tf(f,t,s) = m.beta(6);
            elseif strcmp(coeftype,'tscore')
                TllrXpupilS_tf(f,t,s) = m.tstat.t(6);
            end
        end
    end
end
fprintf('Done.\n')

% ==================================================================
% PULL PER-SAMPLE, TIME-AVERAGED SINGLE-TRIAL MOTOR SIGNAL, Z-SCORED WITHIN SESSION
% ==================================================================
win2 = [0.4 0.6];  % 2nd window relative to sample onset over which to average
swLI = nan(size(tf_ml,1),size(tf_ml,2),size(tf_ml,4));
swLIlpr = nan(size(tf_ml,1),size(tf_ml,2),size(tf_ml,4)-1);
swLIllr = nan(size(tf_ml,1),size(tf_ml,2),size(tf_ml,4)-1);
for s = 1:length(sessions)
    for f = 1:size(tf_ml,2)
        for smp = 1:size(tf_ml,4)
            swLI(sess_full==sessions(s),f,smp) = squeeze(zscore(mean(tf_ml(sess_full==sessions(s),f,smptimes>=win2(1) & smptimes<=win2(2),smp),3)));
        end
        for smp = 2:size(tf_ml,4)
            m = regstats(squeeze(zscore(mean(tf_ml(sess_full==sessions(s),f,smptimes>=win2(1) & smptimes<=win2(2),smp),3))),[LLR_full(sess_full==sessions(s),smp) LLR_full(sess_full==sessions(s),smp).*surprise_full(sess_full==sessions(s),smp)],'linear',{'r'});
            swLIlpr(sess_full==sessions(s),f,smp-1) = zscore(m.r);  % storing LI with current LLR & LLR*surprise regressed out (i.e. isolating a measure of pure prior encoding in LI)
            
            m = regstats(squeeze(zscore(mean(tf_ml(sess_full==sessions(s),f,smptimes>=win2(1) & smptimes<=win2(2),smp),3))),[prior_full(sess_full==sessions(s),smp-1)],'linear',{'r'});
            swLIllr(sess_full==sessions(s),f,smp-1) = zscore(m.r);  % storing LI with prior belief regressed out (i.e. a measure of sample-wise *change* in LI)
        end
    end
end

% ==================================================================
% SAVE RESULTS
% ==================================================================
if strcmp(surprisetype,'pCP')
    if strcmp(coeftype,'beta')
        savename = [savepath,subject,'_samplewise_output_appML_pCP_beta.mat'];
    elseif strcmp(coeftype,'tscore')
        savename = [savepath,subject,'_samplewise_output_appML_pCP.mat'];
    end
else
    if strcmp(coeftype,'beta')
        savename = [savepath,subject,'_samplewise_output_appML_beta.mat'];
    elseif strcmp(coeftype,'tscore')
        savename = [savepath,subject,'_samplewise_output_appML.mat'];
    end
end
save(savename,'smptimes','freqs','smp_tf_avg','grad','cfg',...
    'TpriorS_tf','TllrS_tf','TllrXsurpriseS_tf','TllrXuncertS_tf','TllrXpupilS_tf','Rsq_DV_tf','Rsq_evidence_tf','Rsq_DVpure_tf','Rsq_evidencepure_tf','Rsq_DVpure_tf_BIC','Rsq_evidencepure_tf_BIC',...
    'Rsq_DVpure_tf_BIC1','Rsq_evidencepure_tf_BIC1','Rsq_Lpure_tf_BIC1',...
    'win2','swLI','swLIlpr','swLIllr','LPR_full','LLR_full','surprise_full','psi_full','prior_full','dil_full','sess_full','choices_full','pswitch_full','fdist_full','stimIn_full','distseq_full','t_ml')



