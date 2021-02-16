function [] = runFreqAnalysis_appMotorLocalizer_CPcentered_tonly(n)
  % Applies two baselining operations (pre-trial and peri-sample) to MEG
  % TF data and regresses these data (channel*time) onto model-based
  % variables of interest

allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};
subject = allsubj{n};

basetype = 'dB_common';  % options: 'dB_st', 'dB_common', 'pc_st' or 'pc_common'

modeltype = 'fitted_linIU';  % switch b/w 'normative', 'fitted', 'fitted_np', 'fitted_lin', 'fitted_npIU' & 'fitted_linIU'
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

loadpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/TF/wMotorLoc/',basetype,filesep,str2,filesep];

if strcmp(surprisetype,'pCP'), strs='_pCP'; else strs=''; end

% ==================================================================
% TRIAL-WISE ANALYSIS
% ==================================================================
fprintf('Beginning trial-wise analysis...\n')
load([loadpath,subject,'_samplewise_output_appML',strs,'.mat'],'LPR_full','LLR_full','surprise_full','psi_full','prior_full','dil_full','sess_full',...
    'choices_full','pswitch_full','fdist_full','distseq_full','win2');   % load sample-wise computational variables

load([loadpath,subject,'_CPcentered_output_appML.mat'],'psiCP','lprCP','psiCPleak','lprCPleak','lprCPperf')  % load already-computed belief time-courses binned by CP position (info not fullt available in existing t_only .mats)

load([loadpath,subject,'_trialwise_output_appML',strs,'.mat'],'trltimes','freqs','grad','cfg','t_ml');   % load trial-wise data - lo freq
t_mlL = t_ml;
load([loadpath,subject,'_trialwise_output_appML_HiFreq.mat'],'t_ml');   % load trial-wise data - hi freq
t_ml = squeeze(sum(cat(3,t_mlL,t_ml),3));  % sum across freqs (already weighted in previous processing step)

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

trl_data = t_ml(singleCPts,:);
pswitch_full = pswitch_full(singleCPts,:);
fdist_full = fdist_full(singleCPts);
surprise_full = surprise_full(singleCPts,:);

% Match signs of final generative distributions to enable trial pooling
trl_data(fdist_full==0,:) = trl_data(fdist_full==0,:).*-1;

% Average trials within each CP position
dataCP=[]; beliefCP=[]; beliefCPleak=[]; beliefCPperf=[]; trlnums=[];
for s = 1:size(pswitch_full,2)
    if s==1  % no-CP trials
        dataCP(s,:) = mean(trl_data(sum(pswitch_full,2)==0,:),1);
        trlnums(s) = length(find(sum(pswitch_full,2)==0));
    else
        if ~isempty(find(pswitch_full(:,s)==1))
            dataCP(s,:) = mean(trl_data(pswitch_full(:,s)==1,:),1);
        else
            dataCP(s,:) = nan(1,size(dataCP,2));
        end
        trlnums(s) = length(find(pswitch_full(:,s)==1));
    end
end


% ==================================================================
% SAVE RESULTS
% ==================================================================
save([loadpath,subject,'_CPcentered_output_appML_tonly.mat'],'trltimes','dataCP','psiCP','lprCP','psiCPleak','lprCPleak','lprCPperf','trlnums')


