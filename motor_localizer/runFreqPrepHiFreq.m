function [] = runFreqPrepHiFreq(n)
  % Apply planar gradient transformation and calculate TF representations of
  % MEG power for each of two gradiometers per sensor/trial, using a single
  % Hanning taper for freq range 3-35Hz (window length: 400 ms, step size:
  % 50 ms, freq resolution: 2.5 Hz, bin size: 1 Hz). After TFR calculation,
  % power estimates from the two planar gradiometers for each sensor are
  % combined by taking their sum.

  
% Subject/model variant stuff
allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};
sessions = {'1' '2' '3'}; 
recordings = {'01' '02'}; 

% Construct cell arrays for calling jobs
subj_in={}; sess_in={}; rec_in={};
for i = 1:length(allsubj);
    for j = 1:length(sessions);
        for k = 1:length(recordings);
            subj_in{end+1} = allsubj{i};
            sess_in{end+1} = sessions{j};
            rec_in{end+1} = recordings{k};
        end
    end
end

% QNV
subj_in{end+1} = 'QNV'; sess_in{end+1} = '4'; rec_in{end+1} = '01';
subj_in{end+1} = 'QNV'; sess_in{end+1} = '4'; rec_in{end+1} = '02';

%Subjects 3 recordings
% EMB, Session 2
% rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-8
subj_in{end+1} = 'EMB'; sess_in{end+1} = '2'; rec_in{end+1} = '03';

% JTB, Session 2
% rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-7
subj_in{end+1} = 'JTB'; sess_in{end+1} = '2'; rec_in{end+1} = '03';

% PDP, Session 2
% rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-7
subj_in{end+1} = 'PDP'; sess_in{end+1} = '2'; rec_in{end+1} = '03';

% TNB, Session 3
% rec 01: block 1, rec 02: blocks 2-5, rec 03: blocks 6-9
subj_in{end+1} = 'TNB'; sess_in{end+1} = '3'; rec_in{end+1} = '03';

% TSJ, Session 2
% rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-7
subj_in{end+1} = 'TSJ'; sess_in{end+1} = '2'; rec_in{end+1} = '03';

% Pulling 
subject = subj_in{n}; session = sess_in{n}; recording = rec_in{n};

pause(mod(n,30)*15-1)  % implement file-specific pause time to avoid heavy load on cluster
  
% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults 

megpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed/Data/';  % path of preprocessed MEG data
behavpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Data/',subject,filesep];  % path of behavioural data
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed/';

subjfile = dir([megpath,filesep,subject,'-',session,'*',recording,'.mat']);  % pull current meg filename

% Pupil
loadpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Pupil_Analysis/Pupil/3.Interpolated/';
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Pupil_Analysis/Pupil/'));
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Pupil_Analysis/Behaviour
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Pupil_Analysis/Gen_fun

freqs = [0.06 6];  % filter cutoffs [lo hi]
load([loadpath,'Pupil_peak_surprise_encoding.mat'])
tsampav = 0.54;
tsamp = maxts(strcmp(allsubj,subject));  % subject-specific time-point relative to sample onset at which to take pupil measure (derived from pupil surprise-encoding analysis)

% Modelling
modelpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Modelling/';
modeltype = 'Glaze_npLLR';   % GLAZE model variant for calculating fitted belief/surprise metrics
modeltype_np = 'FitPhi_data_npLLR2';   % FULLY NON-PARAMETRIC model variant for calculating fitted belief/surprise metrics
modeltype_linLLR = 'Glaze_basic';   % SIMPLE GLAZE model variant where dot->LLR transfer function is constrained to be linear

% ==================================================================
% LOAD CURRENT MEG/BEHAVIOURAL/PUPIL DATA & DERIVE COMPUTATIONAL VARIABLES
% ==================================================================
fprintf('\nLoading model fit: %s...\n',subject)
load([modelpath,modeltype_np,filesep,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
if strcmp(modeltype_np,'FitPhi_data')
    phi_n = pm_fit(1:end-1);
elseif strcmp(modeltype_np,'FitPhi_data_pwLLR')
    phi_n = pm_fit(1:end-4);
elseif strcmp(modeltype_np,'FitPhi_data_npLLR2')
    phi_n = [fliplr(-pm_fit(1:length(L_nm1))') 0 pm_fit(1:length(L_nm1))'];
    L_nm1 = [sort(-L_nm1) 0 L_nm1];
end

fprintf('\nLoading meg file: %s...\n',subjfile.name)
load([megpath,subjfile.name])  % load meg data

fprintf('\nLoading and processing behavioural data...\n')
stimIn_full=[]; choices_full=[]; pswitch_full=[]; fdist_full=[]; distseq_full=[];
LLR_full=[]; psi_full=[]; LPR_full=[]; surprise_full=[]; surpriseS_full=[]; surpriseW_full=[]; deltaL_full=[];
LLR_M_full=[]; psi_M_full=[]; LPR_M_full=[]; surprise_M_full=[]; surpriseS_M_full=[]; surpriseW_M_full=[]; deltaL_M_full=[];
LLR_Mnp_full=[]; psi_Mnp_full=[]; LPR_Mnp_full=[]; surprise_Mnp_full=[]; surpriseS_Mnp_full=[]; surpriseW_Mnp_full=[]; deltaL_Mnp_full=[];
LLR_Mlin_full=[]; psi_Mlin_full=[]; LPR_Mlin_full=[]; surprise_Mlin_full=[]; surpriseS_Mlin_full=[]; surpriseW_Mlin_full=[]; deltaL_Mlin_full=[];
dil_full=[]; X_full=[]; Y_full=[]; dil_fullav=[]; X_fullav=[]; Y_fullav=[]; dil_raw_full=[];
for b = unique(data.trialinfo(:,1))'  % looping through each block within this meg dataset
    load([behavpath,'S',session,filesep,'Behaviour',filesep,subject,'_',session,'_',num2str(b),'.mat'])
    load([behavpath,'S',session,filesep,'Sample_seqs',filesep,subject,'_',session,'_',num2str(b),'.mat'])
    
    Behav = Behav(unique(data.trialinfo(data.trialinfo(:,1)==b,2)),:);  % dumping any trials not contained in meg dataset
    stimIn = stimIn(unique(data.trialinfo(data.trialinfo(:,1)==b,2)),:);
    pswitch = pswitch(unique(data.trialinfo(data.trialinfo(:,1)==b,2)),:);
    distseqs = distseqs(unique(data.trialinfo(data.trialinfo(:,1)==b,2)),:);
    
    % Converting sample and choice values to appropriate signs for choice regressions
    stimIn = round(stimIn.*-1);
    choices = Behav(:,2)-1;
    
    % Convert stimulus values to LLRs & calculate sample-wise surprise - NORMATIVE GLAZE
    LLRin = log(normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2))./normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
    [LPR,surprise,psi] = accGlaze_fast(LLRin,gen.H,0,'DY',[]);
    [~,surpriseS] = accGlaze_fast(LLRin,gen.H,0,'scaled_prior',[]);
    [~,deltaL] = accGlaze_fast(LLRin,gen.H,0,'absL',[]);
    [~,surpriseW] = accGlaze_fast(LLRin,gen.H,0,'DY_prior_weighted',[]);
    
    % Convert stimulus values to LLRs & calculate sample-wise surprise - FITTED GLAZE
    load([modelpath,modeltype,filesep,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
    
    if isempty(strfind(modeltype,'npLLR')), rLLR = []; end  % empty variable to pass to getLLR if LLR function is not non-parametric
    LLRinM = getLLR(LLRin,pm_fit,modeltype,rLLR);  % calculate model-derived sample-wise LLRs
    
    [LPR_M,surprise_M,psi_M] = accGlaze_fast(LLRinM,pm_fit(1),0,'DY',[]);
    [~,surpriseS_M] = accGlaze_fast(LLRinM,pm_fit(1),0,'scaled_prior',[]);
    [~,deltaL_M] = accGlaze_fast(LLRinM,pm_fit(1),0,'absL',[]);
    [~,surpriseW_M] = accGlaze_fast(LLRinM,pm_fit(1),0,'DY_prior_weighted',[]);
    
    % Convert stimulus values to LLRs & calculate sample-wise surprise - FITTED NON-PARAMETRIC
    load([modelpath,modeltype_np,filesep,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
    if strcmp(modeltype_np,'FitPhi_data')
        phi_n = pm_fit(1:end-1);
    elseif strcmp(modeltype_np,'FitPhi_data_pwLLR')
        phi_n = pm_fit(1:end-4);
    elseif strcmp(modeltype_np,'FitPhi_data_npLLR2')
        phi_n = [fliplr(-pm_fit(1:length(L_nm1))') 0 pm_fit(1:length(L_nm1))'];
        L_nm1 = [sort(-L_nm1) 0 L_nm1];
    end
    
    if isempty(strfind(modeltype_np,'npLLR')), rLLR = []; end  % empty variable to pass to getLLR if LLR function is not non-parametric
    LLRinMnp = getLLR(LLRin,pm_fit,modeltype_np,rLLR);  % calculate model-derived sample-wise LLRs
    
    [LPR_Mnp,surprise_Mnp,psi_Mnp] = accPhi_fast(LLRinMnp,L_nm1,phi_n,'DY');
    [~,surpriseS_Mnp] = accPhi_fast(LLRinMnp,L_nm1,phi_n,'scaled_prior');
    [~,deltaL_Mnp] = accPhi_fast(LLRinMnp,L_nm1,phi_n,'absL');
    [~,surpriseW_Mnp] = accPhi_fast(LLRinMnp,L_nm1,phi_n,'DY_prior_weighted');
    
    % Convert stimulus values to LLRs & calculate sample-wise surprise - FITTED GLAZE
    load([modelpath,modeltype_linLLR,filesep,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
    
    LLRinMlin = LLRin.*pm_fit(2);  % calculate model-derived sample-wise LLRs
    
    [LPR_Mlin,surprise_Mlin,psi_Mlin] = accGlaze_fast(LLRinMlin,pm_fit(1),0,'DY',[]);
    [~,surpriseS_Mlin] = accGlaze_fast(LLRinMlin,pm_fit(1),0,'scaled_prior',[]);
    [~,deltaL_Mlin] = accGlaze_fast(LLRinMlin,pm_fit(1),0,'absL',[]);
    [~,surpriseW_Mlin] = accGlaze_fast(LLRinMlin,pm_fit(1),0,'DY_prior_weighted',[]);
    
    % Collating useable single trials
    stimIn_full = [stimIn_full; stimIn];        % sample sequences
    choices_full = [choices_full; choices];     % trial-by-trial choices
    pswitch_full = [pswitch_full; pswitch];     % change-point positions
    fdist_full = [fdist_full; Behav(:,1)-1];    % generative dist @ trial end
    distseq_full = [distseq_full; distseqs-1];     % sequences of generative distributions
    
    LLR_full = [LLR_full; LLRin];               % sample evidence strength
    psi_full = [psi_full; psi];                 % evolving scaled prior
    LPR_full = [LPR_full; LPR];                 % evolving belief
    surprise_full = [surprise_full; surprise];  % surprise
    surpriseS_full = [surpriseS_full; surpriseS];  % surprise derived from Glaze-scaled prior
    surpriseW_full = [surpriseW_full; surpriseW];  % option-combined, weighted surprise measure
    deltaL_full = [deltaL_full; deltaL];        % change in belief
    
    LLR_M_full = [LLR_M_full; LLRinM];               % sample evidence strength
    psi_M_full = [psi_M_full; psi_M];                 % evolving scaled prior
    LPR_M_full = [LPR_M_full; LPR_M];                 % evolving belief
    surprise_M_full = [surprise_M_full; surprise_M];  % surprise (fitted)
    surpriseS_M_full = [surpriseS_M_full; surpriseS_M];  % surprise derived from Glaze-scaled prior (fitted)
    surpriseW_M_full = [surpriseW_M_full; surpriseW_M];  % option-combined, weighted surprise measure (fitted)
    deltaL_M_full = [deltaL_M_full; deltaL_M];        % change in belief (fitted)
    
    LLR_Mnp_full = [LLR_Mnp_full; LLRinMnp];               % sample evidence strength
    psi_Mnp_full = [psi_Mnp_full; psi_Mnp];                 % evolving scaled prior
    LPR_Mnp_full = [LPR_Mnp_full; LPR_Mnp];                 % evolving belief
    surprise_Mnp_full = [surprise_Mnp_full; surprise_Mnp];  % surprise (fitted)
    surpriseS_Mnp_full = [surpriseS_Mnp_full; surpriseS_Mnp];  % surprise derived from Glaze-scaled prior (fitted)
    surpriseW_Mnp_full = [surpriseW_Mnp_full; surpriseW_Mnp];  % option-combined, weighted surprise measure (fitted)
    deltaL_Mnp_full = [deltaL_Mnp_full; deltaL_Mnp];        % change in belief (fitted)
    
    LLR_Mlin_full = [LLR_Mlin_full; LLRinMlin];               % sample evidence strength
    psi_Mlin_full = [psi_Mlin_full; psi_Mlin];                 % evolving scaled prior
    LPR_Mlin_full = [LPR_Mlin_full; LPR_Mlin];                 % evolving belief
    surprise_Mlin_full = [surprise_Mlin_full; surprise_Mlin];  % surprise (fitted)
    surpriseS_Mlin_full = [surpriseS_Mlin_full; surpriseS_Mlin];  % surprise derived from Glaze-scaled prior (fitted)
    surpriseW_Mlin_full = [surpriseW_Mlin_full; surpriseW_Mlin];  % option-combined, weighted surprise measure (fitted)
    deltaL_Mlin_full = [deltaL_Mlin_full; deltaL_Mlin];        % change in belief (fitted)
    
    
    % Pupil data
    dil_samp=[]; X_samp=[]; Y_samp=[];
    dil_sampav=[]; X_sampav=[]; Y_sampav=[];
    
    % Load eyetracker file
    % load variable data with different name
    pd = 'data';
    pupil_data = load([loadpath,subject,'_',session,'_',num2str(b),'_interp.mat'], pd);
    pupil_data = pupil_data.(pd);
    pupil = zscore(pupil_data.pupil);  % z-scoring within-block
    times = pupil_data.times;
    
    new_events=[]; new_eventsmp=[]; megtrials = unique(data.trialinfo(data.trialinfo(:,1)==b,2));
    for t = 1:length(megtrials)
        new_events = [new_events; pupil_data.event(megtrials(t),:)];
        new_eventsmp = [new_eventsmp; pupil_data.eventsmp(megtrials(t),:)];
    end
    pupil_data.event = new_events; pupil_data.eventsmp = new_eventsmp;
    %pupil_data.event = pupil_data.event(unique(data.trialinfo(data.trialinfo(:,1)==b,2)),:);  % dumping any trials not contained in meg dataset
    
    if length(pupil_data.event(:,1))>length(choices), pupil_data.event = pupil_data.event(1:length(choices),:); end  % trimming EL trials in case they exceed .mat trials (can happen if block was terminated prematurely)

    % Downsampling EL data to 250Hz (speeds processing and aligns all datasets to same Fs - some were not recorded @ desired 1000Hz)
    if pupil_data.fsample~=250
        fsnew = 250;

        pupil = resample(pupil,fsnew,pupil_data.fsample)';
        pupil_data.Xgaze = resample(pupil_data.Xgaze,fsnew,pupil_data.fsample)';    % X-GAZE regressor
        pupil_data.Ygaze = resample(pupil_data.Ygaze,fsnew,pupil_data.fsample)';    % Y-GAZE regressor
        times = (0:(length(pupil)-1))./fsnew; pupil_data.times = times;  % manually creating new times vector

        pupil_data.event(:,2:4) = round(pupil_data.event(:,2:4).*(fsnew/pupil_data.fsample));
        pupil_data.eventsmp = round(pupil_data.eventsmp.*(fsnew/pupil_data.fsample));
        pupil_data.badsmp = unique(round(pupil_data.badsmp.*(fsnew/pupil_data.fsample)));  % log of all bad samples that were previously interpolated
        pupil_data.badsmp(pupil_data.badsmp==0) = [];  % in case a sample index was rounded to zero

        pupil_data.fsample = fsnew;  % replacing stored sampling rate in data structure
    end

    % Ensuring behavioural and pupil datasets contain same trials
    nsampstest=[];
    for t = 1:length(choices)
        nsampstest(t,1) = sum(~isnan(stimIn(t,:)));
    end
    assert(sum(nsampstest-pupil_data.event(:,1))==0,'Sample mismatch for subject %s, session %d, block %d.',subject,session,b);
    if sum(nsampstest-pupil_data.event(:,1))>0, blah, end

%     % Regress out gaze position from pupil
%     % TODO: TRY WITH AND WITHOUT 
%     m = regstats(pupil, [pupil_data.Xgaze pupil_data.Ygaze], 'linear', {'r'});
%     pupil = zscore(m.r);

    % Isolating useable full-sequence trials
%     ts=[];
%     for t = 1:length(choices)
%         if sum(isnan(stimIn(t,:)))==0 && choices(t)<2, ts(end+1) = t; end   %  && Behav(t,6)==0
%     end

    % Filter
    [bfilt,afilt] = butter(3, freqs(1)*2/pupil_data.fsample, 'high');   % hi-pass
    pupil = filtfilt(bfilt,afilt, pupil);

    [bfilt,afilt] = butter(3, freqs(2)*2/pupil_data.fsample, 'low');   % lo-pass
    pupil = filtfilt(bfilt,afilt, pupil);

    pupil = zscore(pupil);

    % taking 1st derivative of pupil signal
    pupilD1 = diff(pupil).*pupil_data.fsample;
    
    raw_av_win = [0.6 1.4];  % window for calculating sample-wise pupil dilation from raw signal
    raw_base_win = [-0.05 0.05];  % window around sample onset to use as baseline
    
    % Loop through trials
    for t = 1:length(choices)
        % Pull individual sample response
        for smp = 1:size(pupil_data.eventsmp,2)
            if ~isnan(pupil_data.eventsmp(t,smp))
                smp1 = find(times>=(times(pupil_data.eventsmp(t,smp))+tsamp),1,'first');  % subject-specific pupil time point
                dil_samp(t,smp) = pupilD1(smp1);
                X_samp(t,smp) = pupil_data.Xgaze(smp1);
                Y_samp(t,smp) = pupil_data.Ygaze(smp1);
                
                smp1av = find(times>=(times(pupil_data.eventsmp(t,smp))+tsampav),1,'first');  % pupil time point picked from grand-average surprise encoding
                dil_sampav(t,smp) = pupilD1(smp1av);
                X_sampav(t,smp) = pupil_data.Xgaze(smp1av);
                Y_sampav(t,smp) = pupil_data.Ygaze(smp1av);
                
                dil_samp_raw(t,smp) = mean(pupil(times>=(times(pupil_data.eventsmp(t,smp))+raw_av_win(1)) & times<=(times(pupil_data.eventsmp(t,smp))+raw_av_win(2)))) - ...
                    mean(pupil(times>=(times(pupil_data.eventsmp(t,smp))+raw_base_win(1)) & times<=(times(pupil_data.eventsmp(t,smp))+raw_base_win(2))));
            else
                dil_samp(t,smp) = NaN;
                X_samp(t,smp) = NaN;
                Y_samp(t,smp) = NaN;
                
                dil_sampav(t,smp) = NaN;
                X_sampav(t,smp) = NaN;
                Y_sampav(t,smp) = NaN;
                
                dil_samp_raw(t,smp) = NaN;
            end
        end
    end
    
    % Concatenate
    dil_full = [dil_full; dil_samp];     % trial-by-trial average pupil dilation
    X_full = [X_full; X_samp];     % trial-by-trial average X gaze position
    Y_full = [Y_full; Y_samp];     % trial-by-trial average X gaze position
    dil_fullav = [dil_fullav; dil_sampav];     % trial-by-trial average pupil dilation
    X_fullav = [X_fullav; X_sampav];     % trial-by-trial average X gaze position
    Y_fullav = [Y_fullav; Y_sampav];     % trial-by-trial average X gaze position
    dil_raw_full = [dil_raw_full; dil_samp_raw];   % trial-by-trial average pupil dilation (raw signal, not derivative)
end

% ==================================================================
% PULL ONLY FULL-LENGTH TRIALS - consider skipping this to maximize trial counts
% ==================================================================
fprintf('Keeping only full-length trials...\n')
assert(length(choices_full)==length(data.trial),'ERROR: Trial counts in MEG/behaviour are unequal')

ts=[];  % useable trials based on behavioural data
nsamps=[];
for t = 1:length(choices_full)
    nsamps(t,1) = length(find(~isnan(stimIn_full(t,:))));
    if sum(isnan(stimIn_full(t,:)))==0 && choices_full(t)<2, ts(end+1) = t; end
end
assert(isempty(find((data.trialinfo(:,end)-nsamps)~=0, 1)),'ERROR: Mismatch in MEG/behaviour number of samples per trial')

nsampsP=[];
for t = 1:length(choices_full)
    nsampsP(t,1) = length(find(~isnan(dil_full(t,:))));
end
assert(isempty(find((data.trialinfo(:,end)-nsampsP)~=0, 1)),'ERROR: Mismatch in MEG/pupil number of samples per trial')

LLR_full = LLR_full(ts,:);
psi_full = psi_full(ts,:);
LPR_full = LPR_full(ts,:);
surprise_full = surprise_full(ts,:);
surpriseS_full = surpriseS_full(ts,:);
surpriseW_full = surpriseW_full(ts,:);
deltaL_full = deltaL_full(ts,:);

LLR_M_full = LLR_M_full(ts,:);
psi_M_full = psi_M_full(ts,:);
LPR_M_full = LPR_M_full(ts,:);
surprise_M_full = surprise_M_full(ts,:);
surpriseS_M_full = surpriseS_M_full(ts,:);
surpriseW_M_full = surpriseW_M_full(ts,:);
deltaL_M_full = deltaL_M_full(ts,:);

LLR_Mnp_full = LLR_Mnp_full(ts,:);
psi_Mnp_full = psi_Mnp_full(ts,:);
LPR_Mnp_full = LPR_Mnp_full(ts,:);
surprise_Mnp_full = surprise_Mnp_full(ts,:);
surpriseS_Mnp_full = surpriseS_Mnp_full(ts,:);
surpriseW_Mnp_full = surpriseW_Mnp_full(ts,:);
deltaL_Mnp_full = deltaL_Mnp_full(ts,:);

LLR_Mlin_full = LLR_Mlin_full(ts,:);
psi_Mlin_full = psi_Mlin_full(ts,:);
LPR_Mlin_full = LPR_Mlin_full(ts,:);
surprise_Mlin_full = surprise_Mlin_full(ts,:);
surpriseS_Mlin_full = surpriseS_Mlin_full(ts,:);
surpriseW_Mlin_full = surpriseW_Mlin_full(ts,:);
deltaL_Mlin_full = deltaL_Mlin_full(ts,:);

stimIn_full = stimIn_full(ts,:);
choices_full = choices_full(ts,:);
pswitch_full = pswitch_full(ts,:);
fdist_full = fdist_full(ts,:);
distseq_full = distseq_full(ts,:);

dil_full = dil_full(ts,:);
X_full = X_full(ts,:);
Y_full = Y_full(ts,:);
dil_fullav = dil_fullav(ts,:);
X_fullav = X_fullav(ts,:);
Y_fullav = Y_fullav(ts,:);
dil_raw_full = dil_raw_full(ts,:);


cfg             = [];
cfg.trials      = ts;
data = ft_selectdata(cfg, data);

% ==================================================================
% PLANAR GRADIENT TRANSFORMATION
% ==================================================================
fprintf('\nRunning planar gradient transformation...\n')

% define neighbours based on CTF template
cfg                 = [];
cfg.method          = 'template';
cfg.layout          = 'CTF275';
neighbours          = ft_prepare_neighbours(cfg);

% compute planar gradiometers for MEG sensors
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.planarmethod    = 'sincos';
cfg.channel         = 'MEG';
cfg.neighbours      = neighbours;
data                = ft_megplanar(cfg, data);

% ==================================================================
% TIME-FREQUENCY DECOMPOSITION
% ==================================================================
fprintf('\nRunning time-frequency decomposition...\n')

cfg                 = [];
cfg.output          = 'pow';
cfg.channel         = 'MEG';
cfg.method          = 'mtmconvol';   % specifying multi-taper method
cfg.taper           = 'dpss';     % with Hanning taper
cfg.keeptrials      = 'yes';
cfg.keeptapers      = 'no';
cfg.precision       = 'single'; % saves disk space
%cfg.feedback        = 'none'; % improve readability of logfiles

% make nice timebins at each 50 ms, will include 0 point of pre-mask onset; last time point will correspond to go cue time for trial with shortest dot-to-go interval
mintime = data.time{1}(1);
minpos = find(data.trialinfo(:,3)==min(data.trialinfo(:,3))); maxtime = data.time{minpos(1)}(end);
toi = floor(mintime) : 0.05 : ceil(maxtime);
toi(toi < mintime) = []; toi(toi > maxtime) = [];

cfg.toi             = toi; % all time within each locking
cfg.pad             = 10; % pad to a fixed number of seconds before TFR

cfg.foi             = 36:4:120;   % frequencies of interest
cfg.t_ftimwin       = ones(1, length(cfg.foi)) .* 0.25;
% cfg.t_ftimwin       = 3*(1./cfg.foi);  % adaptive time window of 3 cycles per estimated frequency
cfg.tapsmofrq       = ones(1, length(cfg.foi)) .* 6;  % specifies half the spectral concentration, which should be between 3 and 11 times the Rayleigh frequency (1/t_ftimwin)

freq                 = ft_freqanalysis(cfg, data);
% assert(isequal(freq.freq, cfg.foi), '! spectral estimation returned different foi, double-check settings !');

% ==================================================================
% COMBINE PLANAR GRADIOMETERS
% ==================================================================
freq = ft_combineplanar([], freq);

% ==================================================================
% SAVE
% ==================================================================
freq.mdlvars.LLR = LLR_full;
freq.mdlvars.psi = psi_full;
freq.mdlvars.LPR = LPR_full;
freq.mdlvars.surprise = surprise_full;
freq.mdlvars.surpriseS = surpriseS_full;
freq.mdlvars.surpriseW = surpriseW_full;
freq.mdlvars.deltaL = deltaL_full;

freq.mdlvars.LLR_M = LLR_M_full;
freq.mdlvars.psi_M = psi_M_full;
freq.mdlvars.LPR_M = LPR_M_full;
freq.mdlvars.surprise_M = surprise_M_full;
freq.mdlvars.surpriseS_M = surpriseS_M_full;
freq.mdlvars.surpriseW_M = surpriseW_M_full;
freq.mdlvars.deltaL_M = deltaL_M_full;

freq.mdlvars.LLR_Mnp = LLR_Mnp_full;
freq.mdlvars.psi_Mnp = psi_Mnp_full;
freq.mdlvars.LPR_Mnp = LPR_Mnp_full;
freq.mdlvars.surprise_Mnp = surprise_Mnp_full;
freq.mdlvars.surpriseS_Mnp = surpriseS_Mnp_full;
freq.mdlvars.surpriseW_Mnp = surpriseW_Mnp_full;
freq.mdlvars.deltaL_Mnp = deltaL_Mnp_full;

freq.mdlvars.LLR_Mlin = LLR_Mlin_full;
freq.mdlvars.psi_Mlin = psi_Mlin_full;
freq.mdlvars.LPR_Mlin = LPR_Mlin_full;
freq.mdlvars.surprise_Mlin = surprise_Mlin_full;
freq.mdlvars.surpriseS_Mlin = surpriseS_Mlin_full;
freq.mdlvars.surpriseW_Mlin = surpriseW_Mlin_full;
freq.mdlvars.deltaL_Mlin = deltaL_Mlin_full;

freq.mdlvars.stimIn = stimIn_full;
freq.mdlvars.choices = choices_full;
freq.mdlvars.pswitch = pswitch_full;
freq.mdlvars.fdist = fdist_full;
freq.mdlvars.distseq = distseq_full;

freq.pupilvars.dil = dil_full;
freq.pupilvars.X = X_full;
freq.pupilvars.Y = Y_full;
freq.pupilvars.dilav = dil_fullav;
freq.pupilvars.Xav = X_fullav;
freq.pupilvars.Yav = Y_fullav;
freq.pupilvars.dil_raw_full = dil_raw_full;

save([savepath,subject,'-',session,'-',recording,'_TF_HiFreq.mat'],'freq')
    
end