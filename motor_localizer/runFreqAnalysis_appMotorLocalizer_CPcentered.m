function [] = runFreqAnalysis_appMotorLocalizer_CPcentered(n)
  % Applies two baselining operations (pre-trial and peri-sample) to MEG
  % TF data and regresses these data (channel*time) onto model-based
  % variables of interest

allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};
subject = allsubj{n};
  
basewin = [-0.4 -0.2];  % baseline window relative to pre-mask onset (s)
basetype = 'dB_common';  % options: 'dB_st', 'dB_common', 'pc_st' or 'pc_common'

smpwin = [0 1.2];  % window for sample-wise analyses
trlwin = [-0.5 5.8];  % window for full-trial-wise analyses

modeltype = 'fitted_npIU';  % switch b/w 'normative', 'fitted', 'fitted_np', 'fitted_lin', 'fitted_npIU' & 'fitted_linIU'
pupiltype = 'fixed';  % switch between 'fixed' (time-point picked from grand-av surprise encoding) and 'ss' (time-point picked from subject-specific surprise encoding)
priortype = 'psi';  % switch between 'LPR' (prior term will be posterior from previous sample) & 'psi' (prior term will be Glaze-transformed prior for current sample)
surprisetype = 'pCP';  % switch between 'pCP' and 'surprise'

% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun/
addpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts/FMINSEARCHBND')
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

modelpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Modelling/Glaze_npLLR_InconUp/';
load([modelpath,'Fits',filesep,subject,'_fixed_fit.mat']);
leakpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Modelling/Leak_basic/';

if ~exist([savepath,basetype,filesep,str2,filesep],'dir'), mkdir([savepath,basetype,filesep,str2,filesep]), end  % making save directory if doesn't exist
savepath = [savepath,basetype,filesep,str2,filesep];

locpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/TF_Motor/',basetype,'/'];  % path to motor localizer weights

smp_data=[]; trl_data=[];
LLR_full=[]; LPR_full=[]; surprise_full=[]; psi_full=[];
dil_full=[]; X_full=[]; Y_full=[]; choices_full=[]; sess_full=[];
pswitch_full=[]; fdist_full=[]; distseq_full=[]; stimIn_full=[]; LLRinO_full=[];


% --- %%%%%%% --- %
% --- LO-FREQ --- %
% --- %%%%%%% --- %
load([locpath,'Motor_loc_weights.mat'])   % load motor localizer

fprintf('Beginning trial-wise analysis...\n')
subjfiles = dir([megpath,subject,'-*TF.mat']);  % pull all meg filenames for this subject
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
    freqs = freq.freq;
    trltimes = freq.time(freq.time>=trlwin(1) & freq.time<=trlwin(2));
    trl_data = cat(1,trl_data,freq.powspctrm(:,:,freqs>=min(allfreqs) & freqs<=max(allfreqs),freq.time>=trlwin(1) & freq.time<=trlwin(2)));
    
    if f==1
        freqs = freqs(freqs>=min(allfreqs) & freqs<=max(allfreqs));  % trimming freqs to align with motor localizer if necessary (since localizer was fit to pooled freqs)
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
% COMPUTE CHANGE-POINT PROBABILITY (& BELIEF UPDATES FOR SUBSET OF MODELS) IF REQUIRED
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
% PASS FINAL LPR TRHOUGH GLAZE TRANSFER FUNCTION - NB: assumes Glaze-basic fits for simplicity, need to add others if desired
% ==================================================================
load([modelpath,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
psi_full(:,end+1) = LPR_full(:,end)+log(((1-pm_fit(1))/pm_fit(1))+exp(-LPR_full(:,end)))-log(((1-pm_fit(1))/pm_fit(1))+exp(LPR_full(:,end)));

% ==================================================================
% COMPUTE BELIEF TIME COURSES FOR LEAKY & PERFECT ACCUMULATION
% ==================================================================
load([leakpath,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
[LPR_full_leak,~,psi_full_leak] = accLeak_fast(LLRinO_full.*pm_fit(2),pm_fit(1),0,'DY_prior_weighted',pIn,0.08);
psi_full_leak(:,end+1) = LPR_full_leak(:,end).*(1-pm_fit(1));

LPR_full_perf = cumsum(LLRinO_full,2);

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
tf_mlL = squeeze(sum(trl_datax,2));  % leaves trials*freqs*time matrix
clear trl_datax trl_data


% --- %%%%%%% --- %
% --- HI-FREQ --- %
% --- %%%%%%% --- %
trl_data=[];
load([locpath,'Motor_loc_weights.mat'])   % load motor localizer

fprintf('Beginning trial-wise analysis...\n')
subjfiles = dir([megpath,subject,'-*TF_HiFreq.mat']);  % pull all meg filenames for this subject
for f = 1:length(subjfiles)

    load([megpath,subjfiles(f).name])  % load meg data
    sess = str2double(subjfiles(f).name(5));
    freq.time = round(freq.time,2);  % rounding time vector to nearest 2nd decimal - otherwise slight inaccuracies can lead to bad timing later
    freq.freq = round(freq.freq);    % same here
    
    % ==================================================================
    % PULL DESIRED SEGMENTS OF DATA
    % ==================================================================
    fprintf('Concatenating trial-wise data segments for %s, lo-freq...\n',subjfiles(f).name)
    freqs = freq.freq;
    trltimes = freq.time(freq.time>=trlwin(1) & freq.time<=trlwin(2));
    trl_data = cat(1,trl_data,freq.powspctrm(:,:,freqs>=min(allfreqs) & freqs<=max(allfreqs),freq.time>=trlwin(1) & freq.time<=trlwin(2)));
    
    if f==1
        freqs = freqs(freqs>=min(allfreqs) & freqs<=max(allfreqs));  % trimming freqs to align with motor localizer if necessary (since localizer was fit to pooled freqs)
        w_freq = w_freq(ismember(allfreqs,freqs));  % trimming weight vectors/matrices to align freqs with current dataset (since localizer was fit to pooled freqs)
        w_sens_freq = w_sens_freq(:,ismember(allfreqs,freqs));
    end
end
clear freq

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
tf_ml = cat(2,tf_mlL,squeeze(sum(trl_datax,2)));  % leaves trials*freqs*time matrix with ALL FREQUENCIES
clear trl_datax trl_data



% ==================================================================
% MAKE CATEGORICAL SESSION REGRESSORS
% ==================================================================
sessions = unique(sess_full);
sess_r = zeros(length(sess_full),length(sessions)-1);
for s = 1:length(sessions)-1
    sess_r(sess_full==sessions(s+1),s) = ones(length(find(sess_full==sessions(s+1))),1);
end

% ==================================================================
% AVERAGE ACROSS TRIALS PER DIFFERENT CHANGE-POINT POSITIONS
% ==================================================================
% Pull only trials with 0 or 1 change-points
singleCPts = find(sum(pswitch_full,2)<=1);

trl_data = tf_ml(singleCPts,:,:); clear tf_ml
pswitch_full = pswitch_full(singleCPts,:);
fdist_full = fdist_full(singleCPts);
psi_full = psi_full(singleCPts,2:end);  % getting rid of first col here since it's all zeros
LPR_full = LPR_full(singleCPts,:);
surprise_full = surprise_full(singleCPts,:);

psi_full_leak = psi_full_leak(singleCPts,2:end);  % getting rid of first col here since it's all zeros
LPR_full_leak = LPR_full_leak(singleCPts,:);
LPR_full_perf = LPR_full_perf(singleCPts,:);

% Match signs of final generative distributions to enable trial pooling
trl_data(fdist_full==0,:,:) = trl_data(fdist_full==0,:,:).*-1;
psi_full(fdist_full==0,:) = psi_full(fdist_full==0,:).*-1;
LPR_full(fdist_full==0,:) = LPR_full(fdist_full==0,:).*-1;
psi_full_leak(fdist_full==0,:) = psi_full_leak(fdist_full==0,:).*-1;
LPR_full_leak(fdist_full==0,:) = LPR_full_leak(fdist_full==0,:).*-1;
LPR_full_perf(fdist_full==0,:) = LPR_full_perf(fdist_full==0,:).*-1;

% Average trials within each CP position
dataCP=[]; psiCP=[]; lprCP=[]; psiCPleak=[]; lprCPleak=[]; lprCPperf=[]; trlnums=[];
for s = 1:size(pswitch_full,2)
    if s==1  % no-CP trials
        dataCP(s,:,:) = mean(trl_data(sum(pswitch_full,2)==0,:,:),1);
        psiCP(s,:) = mean(psi_full(sum(pswitch_full,2)==0,:),1);
        lprCP(s,:) = mean(LPR_full(sum(pswitch_full,2)==0,:),1);
        psiCPleak(s,:) = mean(psi_full_leak(sum(pswitch_full,2)==0,:),1);
        lprCPleak(s,:) = mean(LPR_full_leak(sum(pswitch_full,2)==0,:),1);
        lprCPperf(s,:) = mean(LPR_full_perf(sum(pswitch_full,2)==0,:),1);
        trlnums(s) = length(find(sum(pswitch_full,2)==0));
    else
        if ~isempty(find(pswitch_full(:,s)==1))
            dataCP(s,:,:) = mean(trl_data(pswitch_full(:,s)==1,:,:),1);
            psiCP(s,:) = mean(psi_full(pswitch_full(:,s)==1,:),1);
            lprCP(s,:) = mean(LPR_full(pswitch_full(:,s)==1,:),1);
            psiCPleak(s,:) = mean(psi_full_leak(pswitch_full(:,s)==1,:),1);
            lprCPleak(s,:) = mean(LPR_full_leak(pswitch_full(:,s)==1,:),1);
            lprCPperf(s,:) = mean(LPR_full_perf(pswitch_full(:,s)==1,:),1);
        else
            dataCP(s,:,:) = nan(1,size(dataCP,2),size(dataCP,3));
            psiCP(s,:) = nan(1,size(psiCP,2));
            lprCP(s,:) = nan(1,size(lprCP,2));
            psiCPleak(s,:) = nan(1,size(psiCPleak,2));
            lprCPleak(s,:) = nan(1,size(lprCPleak,2));
            lprCPperf(s,:) = nan(1,size(lprCPperf,2));
        end
        trlnums(s) = length(find(pswitch_full(:,s)==1));
    end
end


% ==================================================================
% SAVE RESULTS
% ==================================================================
save([savepath,subject,'_CPcentered_output_appML.mat'],'trltimes','dataCP','psiCP','lprCP','psiCPleak','lprCPleak','lprCPperf','trlnums','allfreqs')





