function runSRanalysis_LIhalf_PPI(nsubj,nclust)

allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};  % ORDER IS IMPORTANT FOR TRIAL-IDs!!!!
nsess   = [   3     3     3     3     3     3     3     3     3     3     3     3     3     4     3     3     3];
subj = allsubj{nsubj};

proc = 'F';
type = 'Lateralized';  % possibilities: Averaged, Lateralized, Pair

clusters = {'vfcPrimary';'vfcEarly';'vfcVO';'vfcPHC';'vfcTO';'vfcLO';'vfcV3ab';'vfcIPS01';'vfcIPS23';'vfcFEF';...   % Wang
    'JWG_aIPS';'JWG_IPS_PCeS';'JWG_M1';...                                                                        % JW
    'HCPMMP1_cingulate_pos';'HCPMMP1_paracentral_midcingulate';'HCPMMP1_insular_front_opercular';'HCPMMP1_premotor';'HCPMMP1_dlpfc';'HCPMMP1_frontal_inferior';'HCPMMP1_frontal_orbital_polar';... % Glasser
    'post_medial_frontal';'ant_medial_frontal';'vent_medial_frontal'};  % PM-defined medial PFC clusters

smpwin = [0 1.0];  % window for sample-wise analyses
trlwin = [-0.5 5.801];  % window for full-trial-wise analyses

induced = 1;
modeltype = 'fitted_linIU';  % switch b/w 'normative', 'fitted', 'fitted_np' & 'fitted_lin'
pupiltype = 'fixed';  % switch between 'fixed' (time-point picked from grand-av surprise encoding) and 'ss' (time-point picked from subject-specific surprise encoding)
priortype = 'psi';  % switch between 'LPR' (prior term will be posterior from previous sample) & 'psi' (prior term will be Glaze-transformed prior for current sample)
surprisetype = 'pCP';  % switch between 'pCP' and 'surprise'

coeftype = 'beta';   % switch b/w 'beta' and 'tscore'

% Paths dervied from processing options
if induced, str1 = 'Conv2mne_induced'; else str1 = 'Conv2mne'; end
if strcmp(modeltype,'normative'), str2 = 'output';
elseif strcmp(modeltype,'fitted'), str2 = 'output_fitted';
elseif strcmp(modeltype,'fitted_np'), str2 = 'output_fitted_np';
elseif strcmp(modeltype,'fitted_lin'), str2 = 'output_fitted_lin';
elseif strcmp(modeltype,'fitted_npIU'), str2 = 'output_fitted_npIU';
elseif strcmp(modeltype,'fitted_linIU'), str2 = 'output_fitted_linIU'; end
megpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/',str1,'/agg/'];
varpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed4mne/BehavPupil/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/';
if ~exist([savepath,str1,filesep,'agg',filesep,str2,filesep],'dir'), mkdir([savepath,str1,filesep,'agg',filesep,str2,filesep]), end  % making save directory if doesn't exist
savepath = [savepath,str1,filesep,'agg',filesep,str2,filesep];

% ==================================================================
% LOAD DATA
% ==================================================================
% loop through sessions and load data
trl_data=[];
LLR_full=[]; psi_full=[]; LPR_full=[]; surprise_full=[];
dil_full=[]; X_full=[]; Y_full=[]; choices_full=[]; sess_full=[];
for sess = 1:nsess(strcmp(allsubj,subj))
    
    fname = [megpath,'S',subj,'_SESS',num2str(sess),'_',proc,'_agg.hdf'];
    fprintf('Loading %s...\n',fname)
        
    % pull times & frequencies vectors, trial IDs
    if sess==1
        info = h5info(fname);  % pull information about file structure
        freqs = sort(cellfun(@str2num,[{info.Groups(1).Groups(1).Datasets(:).Name}]));
        times = h5readatt(fname,['/',type,'/',clusters{nclust},'/',num2str(freqs(1)),'/'],'cols_time')';
        times = round(times,2);   % can be very minor error in timings which can throw off indexing - rounding to remove this
    end
    
    % construct data matrix for this cluster (trials*freqs*times)
    sess_data = [];
    for f = 1:length(freqs)
        data_in = h5read(fname,['/',type,'/',clusters{nclust},'/',num2str(freqs(f)),'/'])';
        
        mne_ids = double(h5readatt(fname,['/',type,'/',clusters{nclust},'/',num2str(freqs(f)),'/'],'rows_trial'))'- nsubj - sess; % pull trial IDs (transformed to deal with earlier bug)
        
        if ~issorted(mne_ids)  % of trials aren't in ascending order, sort
            [mne_ids,sortI] = sort(mne_ids);
            data_in = data_in(sortI,:);
        end
        sess_data(:,f,:) = data_in;
    end
    
    % load and concatenate behaviour/pupil
    load([varpath,subj,'-',num2str(sess),'_vars.mat'])
    
    trialID = trialID - nsubj - sess;  % transform trial IDs to deal with earlier bug
    assert(isempty(find((mne_ids'-trialID)~=0, 1)),'ERROR: Mismatch in MEG/behaviour trial IDs')
    assert(length(mdlvars.choices)==size(sess_data,1),'ERROR: Trial counts in MEG/behaviour are unequal')
    
    if strcmp(modeltype,'normative')
        LLR_full = [LLR_full; mdlvars.LLR];
        LPR_full = [LPR_full; mdlvars.LPR];
        psi_full = [psi_full; mdlvars.psi];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; mdlvars.surpriseW];    % PICK SURPRISE TYPE HERE
        else surprise_full = [surprise_full; mdlvars.pCP];
        end
    elseif strcmp(modeltype,'fitted')
        LLR_full = [LLR_full; mdlvars.LLR_M];
        LPR_full = [LPR_full; mdlvars.LPR_M];
        psi_full = [psi_full; mdlvars.psi_M];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; mdlvars.surpriseW_M];     % PICK SURPRISE TYPE HERE
        else surprise_full = [surprise_full; mdlvars.pCP_M];
        end
    elseif strcmp(modeltype,'fitted_np')
        LLR_full = [LLR_full; mdlvars.LLR_Mnp];
        LPR_full = [LPR_full; mdlvars.LPR_Mnp];
        psi_full = [psi_full; mdlvars.psi_Mnp];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; mdlvars.surpriseW_Mnp];     % PICK SURPRISE TYPE HERE
        else surprise_full = [surprise_full; mdlvars.pCP_Mnp];
        end
    elseif strcmp(modeltype,'fitted_lin')
        LLR_full = [LLR_full; mdlvars.LLR_Mlin];
        LPR_full = [LPR_full; mdlvars.LPR_Mlin];
        psi_full = [psi_full; mdlvars.psi_Mlin];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; mdlvars.surpriseW_Mlin];     % PICK SURPRISE TYPE HERE
        else surprise_full = [surprise_full; mdlvars.pCP_Mlin];
        end
    elseif strcmp(modeltype,'fitted_npIU')
        LLR_full = [LLR_full; mdlvars.LLR_MnpIU];
        LPR_full = [LPR_full; mdlvars.LPR_MnpIU];
        psi_full = [psi_full; mdlvars.psi_MnpIU];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; mdlvars.surpriseW_MnpIU];     % PICK SURPRISE TYPE HERE
        else surprise_full = [surprise_full; mdlvars.pCP_MnpIU];
        end
    elseif strcmp(modeltype,'fitted_linIU')
        LLR_full = [LLR_full; mdlvars.LLR_MlinIU];
        LPR_full = [LPR_full; mdlvars.LPR_MlinIU];
        psi_full = [psi_full; mdlvars.psi_MlinIU];
        if strcmp(surprisetype,'surprise')
            surprise_full = [surprise_full; mdlvars.surpriseW_MlinIU];     % PICK SURPRISE TYPE HERE
        else surprise_full = [surprise_full; mdlvars.pCP_MlinIU];
        end
    end
    
    choices_full = [choices_full; mdlvars.choices];
    sess_full = [sess_full; ones(length(mdlvars.choices),1).*sess];
    if strcmp(pupiltype,'ss'), dil_full = [dil_full; pupilvars.dil]; elseif strcmp(pupiltype,'fixed'), dil_full = [dil_full; pupilvars.dilav]; end
    if strcmp(pupiltype,'ss'), X_full = [X_full; pupilvars.X]; elseif strcmp(pupiltype,'fixed'), X_full = [X_full; pupilvars.Xav]; end
    if strcmp(pupiltype,'ss'), Y_full = [Y_full; pupilvars.Y]; elseif strcmp(pupiltype,'fixed'), Y_full = [Y_full; pupilvars.Yav]; end
    
    % trim data to desired epoch length and concatenate across sessions
    sess_data = sess_data(:,:,times>=trlwin(1) & times<=trlwin(2));
    trl_data = cat(1,trl_data,sess_data); clear sess_data
end

% ==================================================================
% MAKE CATEGORICAL SESSION REGRESSORS
% ==================================================================
sessions = unique(sess_full);
sess_r = zeros(length(sess_full),length(sessions)-1);
for s = 1:length(sessions)-1
    sess_r(sess_full==sessions(s+1),s) = ones(length(find(sess_full==sessions(s+1))),1);
end

% ==================================================================
% ENSURE USE OF DESIRED FORM OF PRIOR BELIEF
% ==================================================================
prior_full = LPR_full;
if strcmp(priortype,'psi')
    prior_full(:,1:end-1) = psi_full(:,2:end);  % N.B. this way, final term in LPR_full will always be final, untransformed belief (i.e. the quantity that the model uses to make a decision)
end

% ==================================================================
% SEGMENT DATA AROUND INDIVIDUAL SAMPLES
% ==================================================================
fprintf('Concatenating trial- & sample-wise data segments...\n')
times = times(times>=trlwin(1) & times<=trlwin(2));
onsets = 0.4:0.4:0.4*12; onsets=round(onsets,1);  % vector of all sample onset times relative to pre-mask
smptimes = times(times>=smpwin(1) & times<=smpwin(2));  % getting vector of sample times relative to dot onset

for s = 1:length(onsets)
    stsmp = find(times>=onsets(s)+smpwin(1),1,'first');
    smp_data(:,:,:,s) = trl_data(:,:,stsmp:stsmp+length(smptimes)-1);
end
clear trl_data

% ==================================================================
% THROW OUT ANY BAD RESPONSE TRIALS
% ==================================================================
LLR_full = LLR_full(choices_full<=1,:);
prior_full = prior_full(choices_full<=1,:);
surprise_full = surprise_full(choices_full<=1,:);
sess_r = sess_r(choices_full<=1,:);
smp_data = smp_data(choices_full<=1,:,:,:);
choices_full = choices_full(choices_full<=1,:);

% ==================================================================
% RUN SINGLE-TRIAL REGRESSIONS
% ==================================================================
% Sample-wise regressions, time*freq*space
TfullMod=[];

fprintf('Running PPIs...')
for f = 1:size(smp_data,2)  % looping through freqs
    for t = 1:size(smp_data,3)  % looping through time-points

        full_r=[];
        for s = 1:size(smp_data,4)-1  % looping through samples
            m = regstats(smp_data(:,f,t,s+1),[prior_full(:,s) LLR_full(:,s+1) LLR_full(:,s+1).*nanzscore(surprise_full(:,s+1)) LLR_full(:,s+1).*nanzscore(-abs(prior_full(:,s))) sess_r],'linear',{'r'});
            full_r(:,s) = nanzscore(m.r);
        end
        
        [b,~,m] = glmfit([LLR_full LLR_full(:,2:end).*nanzscore(surprise_full(:,2:end)) LLR_full(:,2:end).*nanzscore(-abs(prior_full(:,1:end-1))) full_r_wC],[choices_full ones(length(choices_full),1)],'binomial');
        if strcmp(coeftype,'beta')
            TfullMod(f,t,1:size(LLR_full,2)-1) = b((end-size(LLR_full,2)+2):end);
        elseif strcmp(coeftype,'tscore')
            TfullMod(f,t,1:size(LLR_full,2)-1) = m.t((end-size(LLR_full,2)+2):end);
        end
    end
end
fprintf('Done.\n')

% ==================================================================
% SAVE RESULTS
% ==================================================================
if strcmp(coeftype,'beta'), bstr='_beta'; else bstr=''; end
if strcmp(surprisetype,'pCP')
    if strcmp(priortype,'psi')
        savename = [savepath,subj,'_',proc,'_',clusters{nclust},'_psi_LIhalf_pCP_PPI',bstr,'.mat'];
    else
        savename = [savepath,subj,'_',proc,'_',clusters{nclust},'_LIhalf_pCP_PPI',bstr,'.mat'];
    end
else
    if strcmp(priortype,'psi')
        savename = [savepath,subj,'_',proc,'_',clusters{nclust},'_psi_LIhalf_PPI',bstr,'.mat'];
    else
        savename = [savepath,subj,'_',proc,'_',clusters{nclust},'_LIhalf_PPI',bstr,'.mat'];
    end
end

save(savename,'smptimes','freqs','TfullMod')

