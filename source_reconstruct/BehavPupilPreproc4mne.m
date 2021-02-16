function BehavPupilPreproc4mne(subject,session)

% Subject list (for getting corresponding subject #; order is important)
allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};
subjnum = find(strcmp(allsubj,subject));
  
% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts/FMINSEARCHBND
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults 

megpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed4mne/';  % path of preprocessed MEG data
behavpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Data/',subject,filesep];  % path of behavioural data
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed4mne/BehavPupil/';

subjfiles = dir([megpath,filesep,subject,'-',session,'*_preproc4mne.mat']);  % pull current meg filename

% Pupil
loadpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Pupil_Analysis/Pupil/3.Interpolated/';
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Pupil_Analysis/Pupil/'));
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Pupil_Analysis/Behaviour
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Pupil_Analysis/Gen_fun

freqs = [0.06 6];  % filter cutoffs [lo hi]
load([loadpath,'Pupil_peak_surprise_encoding.mat'])
tsampav = 0.54;
tsamp = maxts(strcmp(allsubj,subject));  % subject-specific time-point relative to sample onset at which to take pupil measure (derived from pupil surprise-encoding analysis)
raw_av_win = [0.6 1.4];  % window for calculating sample-wise pupil dilation from raw signal
raw_base_win = [-0.05 0.05];  % window around sample onset to use as baseline

% Modelling
modelpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Modelling/';
modeltype = 'Glaze_npLLR';   % GLAZE model variant for calculating fitted belief/surprise metrics
modeltype_np = 'FitPhi_data_npLLR2';   % FULLY NON-PARAMETRIC model variant for calculating fitted belief/surprise metrics
modeltype_linLLR = 'Glaze_basic';   % SIMPLE GLAZE model variant where dot->LLR transfer function is constrained to be linear
modeltype_npIU = 'Glaze_npLLR_InconUp';   % Non-parametric LLR function W/ inconsistent sample upweighting model variant for calculating fitted belief/surprise metrics
modeltype_linIU = 'Glaze_basic_InconUp';

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

stimIn_full=[]; choices_full=[]; pswitch_full=[]; fdist_full=[]; distseq_full=[]; trialinfo_full=[];
LLR_full=[]; psi_full=[]; LPR_full=[]; surprise_full=[]; surpriseS_full=[]; surpriseW_full=[]; deltaL_full=[]; pCP_full=[];
LLR_M_full=[]; psi_M_full=[]; LPR_M_full=[]; surprise_M_full=[]; surpriseS_M_full=[]; surpriseW_M_full=[]; deltaL_M_full=[]; pCP_M_full=[];
LLR_Mnp_full=[]; psi_Mnp_full=[]; LPR_Mnp_full=[]; surprise_Mnp_full=[]; surpriseS_Mnp_full=[]; surpriseW_Mnp_full=[]; deltaL_Mnp_full=[]; pCP_Mnp_full=[];
LLR_Mlin_full=[]; psi_Mlin_full=[]; LPR_Mlin_full=[]; surprise_Mlin_full=[]; surpriseS_Mlin_full=[]; surpriseW_Mlin_full=[]; deltaL_Mlin_full=[]; pCP_Mlin_full=[];
LLR_MnpIU_full=[]; psi_MnpIU_full=[]; LPR_MnpIU_full=[]; surprise_MnpIU_full=[]; surpriseS_MnpIU_full=[]; surpriseW_MnpIU_full=[]; deltaL_MnpIU_full=[]; pCP_MnpIU_full=[];
LLR_MlinIU_full=[]; psi_MlinIU_full=[]; LPR_MlinIU_full=[]; surprise_MlinIU_full=[]; surpriseS_MlinIU_full=[]; surpriseW_MlinIU_full=[]; deltaL_MlinIU_full=[]; pCP_MlinIU_full=[];

dil_full=[]; X_full=[]; Y_full=[]; dil_fullav=[]; X_fullav=[]; Y_fullav=[]; dil_raw_full=[]; base_raw_full=[]; dil_pret_full=[]; base_pret_full=[];
for f = 1:length(subjfiles)
    fprintf('\nLoading meg file: %s...\n',subjfiles(f).name)
    load([megpath,subjfiles(f).name]);  % load meg data
    trialinfo = data.trialinfo; clear data  % keep only trial info matrix
    trialinfo_full = [trialinfo_full; trialinfo]; % concatenate trial infor matrices across recs for sanity checks
    
    fprintf('\nLoading and processing behavioural data...\n')
    for b = unique(trialinfo(:,1))'  % looping through each block within this meg dataset
        load([behavpath,'S',session,filesep,'Behaviour',filesep,subject,'_',session,'_',num2str(b),'.mat'])
        load([behavpath,'S',session,filesep,'Sample_seqs',filesep,subject,'_',session,'_',num2str(b),'.mat'])
        
        Behav = Behav(unique(trialinfo(trialinfo(:,1)==b,2)),:);  % dumping any trials not contained in meg dataset
        stimIn = stimIn(unique(trialinfo(trialinfo(:,1)==b,2)),:);
        pswitch = pswitch(unique(trialinfo(trialinfo(:,1)==b,2)),:);
        distseqs = distseqs(unique(trialinfo(trialinfo(:,1)==b,2)),:);
        
        % Converting sample and choice values to appropriate signs for choice regressions
        stimIn = round(stimIn.*-1);
        choices = Behav(:,2)-1;
        
        % Convert stimulus values to LLRs & calculate sample-wise surprise - NORMATIVE GLAZE
        LLRin = log(normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2))./normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
        pIn = cat(3,normpdf(stimIn,17,29),normpdf(stimIn,-17,29));
        H = 0.08;
        
        [LPR,surprise,psi] = accGlaze_fast(LLRin,gen.H,0,'DY',[]);
        [~,surpriseS] = accGlaze_fast(LLRin,gen.H,0,'scaled_prior',[]);
        [~,deltaL] = accGlaze_fast(LLRin,gen.H,0,'absL',[]);
        [~,surpriseW] = accGlaze_fast(LLRin,gen.H,0,'DY_prior_weighted',[]);
        [~,pCP,~] = accGlaze_fast(LLRin,H,0,'pCP',pIn);
        
        % Convert stimulus values to LLRs & calculate sample-wise surprise - FITTED GLAZE
        load([modelpath,modeltype,filesep,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
        
        if isempty(strfind(modeltype,'npLLR')), rLLR = []; end  % empty variable to pass to getLLR if LLR function is not non-parametric
        LLRinM = getLLR(LLRin,pm_fit,modeltype,rLLR);  % calculate model-derived sample-wise LLRs
        
        gains = cumsum(pm_fit(3:end)');
        gains = [fliplr(gains) 0 gains];  % creating full vector of gain terms, symmetric around LLR=0
        rLLRfull_norm = [sort(-rLLR) 0 rLLR];  % creating full vector of matching original LLRs
        
        LLRext = log(normpdf(-100:1:100,17,29)./normpdf(-100:1:100,-17,29)); % Fit p(stim) from fitted LLRs under constraint that integral of p(stim) within stimulus bounds = 1
        xy = [interp1(LLRext,-100:1:100,rLLRfull_norm,'linear')' (rLLRfull_norm.*gains)'];  % [x=true stimulus positions, y=fitted LLRs]
        log_pstim = fminsearchbnd(@(params) fit_pstim_from_LLR(params,xy), log(normpdf(xy(:,1),gen.mu(1),gen.sigma(1))), ones(1,size(xy,1)).*-60, ones(1,size(xy,1)).*-1);  % fit p(stim) in log space
        stim=xy(:,1)'; pstim = exp(interp1(stim,log_pstim,[min(stim) -90:1:90 max(stim)],'spline')); pstim = pstim./sum(pstim);  % calculate normalized hi-res p(stim)
        
        pIn=[];
        for samp = 1:size(LLRinM,2)
            pIn(:,samp,1) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(pstim),stimIn(:,samp),'spline'));
            pIn(:,samp,2) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(fliplr(pstim)),stimIn(:,samp),'spline'));
        end
        H = pm_fit(1);
        
        [LPR_M,surprise_M,psi_M] = accGlaze_fast(LLRinM,pm_fit(1),0,'DY',[]);
        [~,surpriseS_M] = accGlaze_fast(LLRinM,pm_fit(1),0,'scaled_prior',[]);
        [~,deltaL_M] = accGlaze_fast(LLRinM,pm_fit(1),0,'absL',[]);
        [~,surpriseW_M] = accGlaze_fast(LLRinM,pm_fit(1),0,'DY_prior_weighted',[]);
        [~,pCP_M,~] = accGlaze_fast(LLRinM,H,0,'pCP',pIn);
        
        % Convert stimulus values to LLRs & calculate sample-wise surprise - FITTED NON-PARAMETRIC
        load([modelpath,modeltype_np,filesep,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
        
        gains = cumsum(pm_fit(length(L_nm1)+1:length(L_nm1)+length(rLLR))');
        gains = [fliplr(gains) 0 gains];  % creating full vector of gain terms, symmetric around LLR=0
        rLLRfull = [sort(-rLLR) 0 rLLR];  % creating full vector of matching original LLRs
        
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
        
        LLRext = log(normpdf(-100:1:100,17,29)./normpdf(-100:1:100,-17,29)); % Fit p(stim) from fitted LLRs under constraint that integral of p(stim) within stimulus bounds = 1
        xy = [interp1(LLRext,-100:1:100,rLLRfull,'linear')' (rLLRfull.*gains)'];  % [x=true stimulus positions, y=fitted LLRs]
        log_pstim = fminsearchbnd(@(params) fit_pstim_from_LLR(params,xy), log(normpdf(xy(:,1),gen.mu(1),gen.sigma(1))), ones(1,size(xy,1)).*-60, ones(1,size(xy,1)).*-1);  % fit p(stim) in log space
        stim=xy(:,1)'; pstim = exp(interp1(stim,log_pstim,[min(stim) -90:1:90 max(stim)],'spline')); pstim = pstim./sum(pstim);  % calculate normalized hi-res p(stim)
        
        pIn=[];
        for samp = 1:size(LLRinMnp,2)
            pIn(:,samp,1) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(pstim),stimIn(:,samp),'spline'));
            pIn(:,samp,2) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(fliplr(pstim)),stimIn(:,samp),'spline'));
        end
        H = 0.08;
        
        [LPR_Mnp,surprise_Mnp,psi_Mnp] = accPhi_fast(LLRinMnp,L_nm1,phi_n,'DY');
        [~,surpriseS_Mnp] = accPhi_fast(LLRinMnp,L_nm1,phi_n,'scaled_prior');
        [~,deltaL_Mnp] = accPhi_fast(LLRinMnp,L_nm1,phi_n,'absL');
        [~,surpriseW_Mnp] = accPhi_fast(LLRinMnp,L_nm1,phi_n,'DY_prior_weighted');
        [~,pCP_Mnp,~] = accPhi_fast(LLRinMnp,L_nm1,phi_n,'pCP',pIn,H);
        
        % Convert stimulus values to LLRs & calculate sample-wise surprise - FITTED GLAZE
        load([modelpath,modeltype_linLLR,filesep,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
        
        LLRinMlin = LLRin.*pm_fit(2);  % calculate model-derived sample-wise LLRs
        
        c = 2; while stimIn(1,c)==0, c=c+1; end  % get sample ~= zero
        cLLR = LLRinMlin(1,c); cStim = stimIn(1,c);
        sigma = sqrt((-(cStim-17).^2 + (cStim+17).^2)./(2.*cLLR));
        pIn = cat(3,normpdf(stimIn,17,sigma),normpdf(stimIn,-17,sigma));
        H = pm_fit(1);
        
        [LPR_Mlin,surprise_Mlin,psi_Mlin] = accGlaze_fast(LLRinMlin,pm_fit(1),0,'DY',[]);
        [~,surpriseS_Mlin] = accGlaze_fast(LLRinMlin,pm_fit(1),0,'scaled_prior',[]);
        [~,deltaL_Mlin] = accGlaze_fast(LLRinMlin,pm_fit(1),0,'absL',[]);
        [~,surpriseW_Mlin] = accGlaze_fast(LLRinMlin,pm_fit(1),0,'DY_prior_weighted',[]);
        [~,pCP_Mlin,~] = accGlaze_fast(LLRinMlin,H,0,'pCP',pIn);
        
        
        % Convert stimulus values to LLRs & calculate sample-wise surprise - FITTED GLAZE W/ INCONSISTENCY GAIN FACTOR
        load([modelpath,modeltype_npIU,filesep,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
        
        if isempty(strfind(modeltype_npIU,'npLLR')), rLLR = []; end  % empty variable to pass to getLLR if LLR function is not non-parametric
        gains = cumsum(pm_fit(4:end)');
        gains = [fliplr(gains) 0 gains];  % creating full vector of gain terms, symmetric around LLR=0
        rLLRfull_norm = [sort(-rLLR) 0 rLLR];  % creating full vector of matching original LLRs
        
        LLRext = log(normpdf(-100:1:100,17,29)./normpdf(-100:1:100,-17,29)); % Fit p(stim) from fitted LLRs under constraint that integral of p(stim) within stimulus bounds = 1
        xy = [interp1(LLRext,-100:1:100,rLLRfull_norm,'linear')' (rLLRfull_norm.*gains)'];  % [x=true stimulus positions, y=fitted LLRs]
        log_pstim = fminsearchbnd(@(params) fit_pstim_from_LLR(params,xy), log(normpdf(xy(:,1),gen.mu(1),gen.sigma(1))), ones(1,size(xy,1)).*-60, ones(1,size(xy,1)).*-1);  % fit p(stim) in log space
        stim=xy(:,1)'; pstim = exp(interp1(stim,log_pstim,[min(stim) -90:1:90 max(stim)],'spline')); pstim = pstim./sum(pstim);  % calculate normalized hi-res p(stim)
        
        LLRinMnpIU = nan(size(LLRin)); pIn=[];
        for samp = 1:size(LLRinMnpIU,2)
            LLRinMnpIU(:,samp) = interp1(rLLRfull_norm,rLLRfull_norm.*gains,LLRin(:,samp),'spline');
            pIn(:,samp,1) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(pstim),stimIn(:,samp),'spline'));
            pIn(:,samp,2) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(fliplr(pstim)),stimIn(:,samp),'spline'));
        end
        H = pm_fit(1);
        
        [LPR_MnpIU,surprise_MnpIU,psi_MnpIU] = accGlaze_InconUp_fast(LLRinMnpIU,H,pm_fit(3),0,'DY',[]);
        [~,surpriseS_MnpIU] = accGlaze_InconUp_fast(LLRinMnpIU,H,pm_fit(3),0,'scaled_prior',[]);
        [~,deltaL_MnpIU] = accGlaze_InconUp_fast(LLRinMnpIU,H,pm_fit(3),0,'absL',[]);
        [~,surpriseW_MnpIU] = accGlaze_InconUp_fast(LLRinMnpIU,H,pm_fit(3),0,'DY_prior_weighted',[]);
        [~,pCP_MnpIU,~] = accGlaze_InconUp_fast(LLRinMnpIU,H,pm_fit(3),0,'pCP',pIn);
        
        
        % Convert stimulus values to LLRs & calculate sample-wise surprise - FITTED GLAZE W/ INCONSISTENCY GAIN FACTOR
        load([modelpath,modeltype_linIU,filesep,'Fits',filesep,subject,'_fixed_fit.mat'])  % load meg data
        
        LLRinMlinIU = LLRin.*pm_fit(2);  % calculate model-derived sample-wise LLRs
        
        c = 2; while stimIn(1,c)==0, c=c+1; end  % get sample ~= zero
        cLLR = LLRinMlinIU(1,c); cStim = stimIn(1,c);
        sigma = sqrt((-(cStim-17).^2 + (cStim+17).^2)./(2.*cLLR));
        pIn = cat(3,normpdf(stimIn,17,sigma),normpdf(stimIn,-17,sigma));
        H = pm_fit(1);
        
        [LPR_MlinIU,surprise_MlinIU,psi_MlinIU] = accGlaze_InconUp_fast(LLRinMlinIU,H,pm_fit(4),0,'DY',[]);
        [~,surpriseS_MlinIU] = accGlaze_InconUp_fast(LLRinMlinIU,H,pm_fit(4),0,'scaled_prior',[]);
        [~,deltaL_MlinIU] = accGlaze_InconUp_fast(LLRinMlinIU,H,pm_fit(4),0,'absL',[]);
        [~,surpriseW_MlinIU] = accGlaze_InconUp_fast(LLRinMlinIU,H,pm_fit(4),0,'DY_prior_weighted',[]);
        [~,pCP_MlinIU,~] = accGlaze_InconUp_fast(LLRinMlinIU,H,pm_fit(4),0,'pCP',pIn);
        
        
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
        pCP_full = [pCP_full; log(pCP)];        % change-point probability
        
        LLR_M_full = [LLR_M_full; LLRinM];               % sample evidence strength
        psi_M_full = [psi_M_full; psi_M];                 % evolving scaled prior
        LPR_M_full = [LPR_M_full; LPR_M];                 % evolving belief
        surprise_M_full = [surprise_M_full; surprise_M];  % surprise (fitted)
        surpriseS_M_full = [surpriseS_M_full; surpriseS_M];  % surprise derived from Glaze-scaled prior (fitted)
        surpriseW_M_full = [surpriseW_M_full; surpriseW_M];  % option-combined, weighted surprise measure (fitted)
        deltaL_M_full = [deltaL_M_full; deltaL_M];        % change in belief (fitted)
        pCP_M_full = [pCP_M_full; log(pCP_M)];        % change-point probability
        
        LLR_Mnp_full = [LLR_Mnp_full; LLRinMnp];               % sample evidence strength
        psi_Mnp_full = [psi_Mnp_full; psi_Mnp];                 % evolving scaled prior
        LPR_Mnp_full = [LPR_Mnp_full; LPR_Mnp];                 % evolving belief
        surprise_Mnp_full = [surprise_Mnp_full; surprise_Mnp];  % surprise (fitted)
        surpriseS_Mnp_full = [surpriseS_Mnp_full; surpriseS_Mnp];  % surprise derived from Glaze-scaled prior (fitted)
        surpriseW_Mnp_full = [surpriseW_Mnp_full; surpriseW_Mnp];  % option-combined, weighted surprise measure (fitted)
        deltaL_Mnp_full = [deltaL_Mnp_full; deltaL_Mnp];        % change in belief (fitted)
        pCP_Mnp_full = [pCP_Mnp_full; log(pCP_Mnp)];        % change-point probability
        
        LLR_Mlin_full = [LLR_Mlin_full; LLRinMlin];               % sample evidence strength
        psi_Mlin_full = [psi_Mlin_full; psi_Mlin];                 % evolving scaled prior
        LPR_Mlin_full = [LPR_Mlin_full; LPR_Mlin];                 % evolving belief
        surprise_Mlin_full = [surprise_Mlin_full; surprise_Mlin];  % surprise (fitted)
        surpriseS_Mlin_full = [surpriseS_Mlin_full; surpriseS_Mlin];  % surprise derived from Glaze-scaled prior (fitted)
        surpriseW_Mlin_full = [surpriseW_Mlin_full; surpriseW_Mlin];  % option-combined, weighted surprise measure (fitted)
        deltaL_Mlin_full = [deltaL_Mlin_full; deltaL_Mlin];        % change in belief (fitted)
        pCP_Mlin_full = [pCP_Mlin_full; log(pCP_Mlin)];        % change-point probability
        
        LLR_MnpIU_full = [LLR_MnpIU_full; LLRinMnpIU];               % sample evidence strength
        psi_MnpIU_full = [psi_MnpIU_full; psi_MnpIU];                 % evolving scaled prior
        LPR_MnpIU_full = [LPR_MnpIU_full; LPR_MnpIU];                 % evolving belief
        surprise_MnpIU_full = [surprise_MnpIU_full; surprise_MnpIU];  % surprise (fitted)
        surpriseS_MnpIU_full = [surpriseS_MnpIU_full; surpriseS_MnpIU];  % surprise derived from Glaze-scaled prior (fitted)
        surpriseW_MnpIU_full = [surpriseW_MnpIU_full; surpriseW_MnpIU];  % option-combined, weighted surprise measure (fitted)
        deltaL_MnpIU_full = [deltaL_MnpIU_full; deltaL_MnpIU];        % change in belief (fitted)
        pCP_MnpIU_full = [pCP_MnpIU_full; log(pCP_MnpIU)];        % change-point probability
        
        LLR_MlinIU_full = [LLR_MlinIU_full; LLRinMlinIU];               % sample evidence strength
        psi_MlinIU_full = [psi_MlinIU_full; psi_MlinIU];                 % evolving scaled prior
        LPR_MlinIU_full = [LPR_MlinIU_full; LPR_MlinIU];                 % evolving belief
        surprise_MlinIU_full = [surprise_MlinIU_full; surprise_MlinIU];  % surprise (fitted)
        surpriseS_MlinIU_full = [surpriseS_MlinIU_full; surpriseS_MlinIU];  % surprise derived from Glaze-scaled prior (fitted)
        surpriseW_MlinIU_full = [surpriseW_MlinIU_full; surpriseW_MlinIU];  % option-combined, weighted surprise measure (fitted)
        deltaL_MlinIU_full = [deltaL_MlinIU_full; deltaL_MlinIU];        % change in belief (fitted)
        pCP_MlinIU_full = [pCP_MlinIU_full; log(pCP_MlinIU)];        % change-point probability
        
        % Pupil data
        dil_samp=[]; X_samp=[]; Y_samp=[];
        dil_sampav=[]; X_sampav=[]; Y_sampav=[];
        dil_samp_raw=[]; base_samp_raw=[]; pret_raw=[]; pret_D1=[];
        
        % Load eyetracker file
        % load variable data with different name
        pd = 'data';
        pupil_data = load([loadpath,subject,'_',session,'_',num2str(b),'_interp.mat'], pd);
        pupil_data = pupil_data.(pd);
        pupil = zscore(pupil_data.pupil);  % z-scoring within-block
        times = pupil_data.times;
        
        new_events=[]; new_eventsmp=[]; megtrials = unique(trialinfo(trialinfo(:,1)==b,2));
        for t = 1:length(megtrials)
            new_events = [new_events; pupil_data.event(megtrials(t),:)];
            new_eventsmp = [new_eventsmp; pupil_data.eventsmp(megtrials(t),:)];
        end
        pupil_data.event = new_events; pupil_data.eventsmp = new_eventsmp;
        %pupil_data.event = pupil_data.event(unique(trialinfo(trialinfo(:,1)==b,2)),:);  % dumping any trials not contained in meg dataset
        
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
            pupil_data.badsmp(pupil_data.badsmp==0) = [];  % in case a sample index was rounded to zerofdist_full = [fdist_full; Behav(:,1)-1];    % generative dist @ trial end
            
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
        
        % Loop through trials
        for t = 1:length(choices)
            % calculate pre-trial baseline window
                %          window rel. to first sample onset          minus      inter-sample interval
                pret_win = [-0.1 0]   +   times(pupil_data.eventsmp(t,1)-(pupil_data.eventsmp(t,2)-pupil_data.eventsmp(t,1)));
                pret_raw(t,1) = mean(pupil(times>=pret_win(1) & times<=pret_win(2)));
                pret_D1(t,1) = pupilD1(find(times<=pret_win(2),1,'last'));
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
                    
                    base_samp_raw(t,smp) = mean(pupil(times>=(times(pupil_data.eventsmp(t,smp))+raw_base_win(1)) & times<=(times(pupil_data.eventsmp(t,smp))+raw_base_win(2))));
                else
                    dil_samp(t,smp) = NaN;
                    X_samp(t,smp) = NaN;
                    Y_samp(t,smp) = NaN;
                    
                    dil_sampav(t,smp) = NaN;
                    X_sampav(t,smp) = NaN;
                    Y_sampav(t,smp) = NaN;
                    
                    dil_samp_raw(t,smp) = NaN;
                    
                    base_samp_raw(t,smp) = NaN;
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
        dil_raw_full = [dil_raw_full; dil_samp_raw];   % trial-by-trial average pupil dilation per sample (raw signal, not derivative)
        base_raw_full = [base_raw_full; base_samp_raw];   % trial-by-trial average baseline pupil per sample (raw signal, not derivative)
        dil_pret_full = [dil_pret_full; pret_D1];   % trial-by-trial avg pre-trial pupil dilation (derivative)
        base_pret_full = [base_pret_full; pret_raw];   % trial-by-trial avg pre-trial baseline pupil (raw signal, not derivative)
    end
    
end

% ==================================================================
% SANITY CHECKS
% ==================================================================
fprintf('Keeping only full-length trials...\n')
assert(length(choices_full)==size(trialinfo_full,1),'ERROR: Trial counts in MEG/behaviour are unequal')

ts=[];  % useable trials based on behavioural data
nsamps=[];
for t = 1:length(choices_full)
    nsamps(t,1) = length(find(~isnan(stimIn_full(t,:))));
    if sum(isnan(stimIn_full(t,:)))==0 && choices_full(t)<2, ts(end+1) = t; end
end
assert(isempty(find((trialinfo_full(:,end-1)-nsamps)~=0, 1)),'ERROR: Mismatch in MEG/behaviour number of samples per trial')

nsampsP=[];
for t = 1:length(choices_full)
    nsampsP(t,1) = length(find(~isnan(dil_full(t,:))));
end
assert(isempty(find((trialinfo_full(:,end-1)-nsampsP)~=0, 1)),'ERROR: Mismatch in MEG/pupil number of samples per trial')

% ==================================================================
% SAVE
% ==================================================================
trialID = trialinfo_full(:,end);  % saving unique trial identifier in mne-converted data (not there's a bug in this, but can still be used within-subject for appropriate trial matching)

mdlvars.LLR = LLR_full;
mdlvars.psi = psi_full;
mdlvars.LPR = LPR_full;
mdlvars.surprise = surprise_full;
mdlvars.surpriseS = surpriseS_full;
mdlvars.surpriseW = surpriseW_full;
mdlvars.deltaL = deltaL_full;
mdlvars.pCP = pCP_full;

mdlvars.LLR_M = LLR_M_full;
mdlvars.psi_M = psi_M_full;
mdlvars.LPR_M = LPR_M_full;
mdlvars.surprise_M = surprise_M_full;
mdlvars.surpriseS_M = surpriseS_M_full;
mdlvars.surpriseW_M = surpriseW_M_full;
mdlvars.deltaL_M = deltaL_M_full;
mdlvars.pCP_M = pCP_M_full;

mdlvars.LLR_Mnp = LLR_Mnp_full;
mdlvars.psi_Mnp = psi_Mnp_full;
mdlvars.LPR_Mnp = LPR_Mnp_full;
mdlvars.surprise_Mnp = surprise_Mnp_full;
mdlvars.surpriseS_Mnp = surpriseS_Mnp_full;
mdlvars.surpriseW_Mnp = surpriseW_Mnp_full;
mdlvars.deltaL_Mnp = deltaL_Mnp_full;
mdlvars.pCP_Mnp = pCP_Mnp_full;

mdlvars.LLR_Mlin = LLR_Mlin_full;
mdlvars.psi_Mlin = psi_Mlin_full;
mdlvars.LPR_Mlin = LPR_Mlin_full;
mdlvars.surprise_Mlin = surprise_Mlin_full;
mdlvars.surpriseS_Mlin = surpriseS_Mlin_full;
mdlvars.surpriseW_Mlin = surpriseW_Mlin_full;
mdlvars.deltaL_Mlin = deltaL_Mlin_full;
mdlvars.pCP_Mlin = pCP_Mlin_full;

mdlvars.LLR_MnpIU = LLR_MnpIU_full;
mdlvars.psi_MnpIU = psi_MnpIU_full;
mdlvars.LPR_MnpIU = LPR_MnpIU_full;
mdlvars.surprise_MnpIU = surprise_MnpIU_full;
mdlvars.surpriseS_MnpIU = surpriseS_MnpIU_full;
mdlvars.surpriseW_MnpIU = surpriseW_MnpIU_full;
mdlvars.deltaL_MnpIU = deltaL_MnpIU_full;
mdlvars.pCP_MnpIU = pCP_MnpIU_full;

mdlvars.LLR_MlinIU = LLR_MlinIU_full;
mdlvars.psi_MlinIU = psi_MlinIU_full;
mdlvars.LPR_MlinIU = LPR_MlinIU_full;
mdlvars.surprise_MlinIU = surprise_MlinIU_full;
mdlvars.surpriseS_MlinIU = surpriseS_MlinIU_full;
mdlvars.surpriseW_MlinIU = surpriseW_MlinIU_full;
mdlvars.deltaL_MlinIU = deltaL_MlinIU_full;
mdlvars.pCP_MlinIU = pCP_MlinIU_full;

mdlvars.stimIn = stimIn_full;
mdlvars.choices = choices_full;
mdlvars.pswitch = pswitch_full;
mdlvars.fdist = fdist_full;
mdlvars.distseq = distseq_full;

pupilvars.dil = dil_full;
pupilvars.X = X_full;
pupilvars.Y = Y_full;
pupilvars.dilav = dil_fullav;
pupilvars.Xav = X_fullav;
pupilvars.Yav = Y_fullav;
pupilvars.dilraw = dil_raw_full;
pupilvars.baseraw = base_raw_full;
pupilvars.dilpret = dil_pret_full;
pupilvars.basepret = base_pret_full;

save([savepath,subject,'-',session,'_vars.mat'],'mdlvars','pupilvars','trialID')
    
