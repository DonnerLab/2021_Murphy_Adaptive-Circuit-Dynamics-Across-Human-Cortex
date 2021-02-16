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

if strcmp(surprisetype,'pCP'), surpstr = '_pCP'; else surpstr = ''; end
if strcmp(coeftype,'beta'), bstr = '_beta'; else bstr = ''; end


%%%%%%%%%%%%%%%%%%%
%%% Sample-wise %%%
%%%%%%%%%%%%%%%%%%%
sess_r = [2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2];  % # session regressors per subject
pred_r = [2 3 3 3 3 3 3 3 3 3 2];  % # other predictors per sample postion

fprintf('Loading sample-wise data...\nSubject ')

samps = 2:12;  % samples to include for pupil*LLR analyses
DVpure_nLL_sum=[]; EVpure_nLL_sum=[]; DVpure1_nLL_sum=[]; EVpure1_nLL_sum=[]; Lpure1_nLL_sum=[];
for s = 1:length(allsubj)
    fprintf('%d, ',s)
    % Load data
    load([loadpath,allsubj{s},'_samplewise_output_appML',surpstr,bstr,'.mat'])
    freqsL = freqs;
        
    for smp = 1:size(Rsq_DVpure_tf_BIC,3)
        DVpure_nLL_L(:,:,smp) = (Rsq_DVpure_tf_BIC(:,:,smp)-((1+pred_r(smp)+sess_r(s))*log(size(LLR_full,1))))./2;  % retrieve -LL for selected regression models, used later for superBIC
        EVpure_nLL_L(:,:,smp) = (Rsq_evidencepure_tf_BIC(:,:,smp)-((1+pred_r(smp)+sess_r(s))*log(size(LLR_full,1))))./2;
        
        DVpure1_nLL_L(:,:,smp) = (Rsq_DVpure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;  % retrieve -LL for selected regression models, used later for superBIC
        EVpure1_nLL_L(:,:,smp) = (Rsq_evidencepure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;
        Lpure1_nLL_L(:,:,smp) = (Rsq_Lpure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;
    end
        
    load([loadpath,allsubj{s},'_samplewise_output_appML_HiFreq',surpstr,bstr,'.mat'])
    allfreqs = [freqsL freqs];
    
    for smp = 1:size(Rsq_DVpure_tf_BIC,3)
        DVpure_nLL(:,:,smp) = (Rsq_DVpure_tf_BIC(:,:,smp)-((1+pred_r(smp)+sess_r(s))*log(size(LLR_full,1))))./2;  % retrieve -LL for selected regression models, used later for superBIC
        EVpure_nLL(:,:,smp) = (Rsq_evidencepure_tf_BIC(:,:,smp)-((1+pred_r(smp)+sess_r(s))*log(size(LLR_full,1))))./2;
        
        DVpure1_nLL(:,:,smp) = (Rsq_DVpure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;  % retrieve -LL for selected regression models, used later for superBIC
        EVpure1_nLL(:,:,smp) = (Rsq_evidencepure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;
        Lpure1_nLL(:,:,smp) = (Rsq_Lpure_tf_BIC1(:,:,smp)-((1+1+sess_r(s))*log(size(LLR_full,1))))./2;
        
        npreds(smp) = 1+pred_r(smp)+sess_r(s);
        npreds1(smp) = 1+1+sess_r(s);
    end
    
    nObs(s) = size(LLR_full,1);  % recording # observations for use in later super BIC
    nPreds(s) = mean(npreds(samps-1)); nPreds1(s) = mean(npreds1(samps-1));  % recording average # predictors for use in later super BIC
    
    m = [squeeze(mean(DVpure_nLL_L(:,:,samps-1),3)); squeeze(mean(DVpure_nLL(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    DVpure_nLL_sum = sum(cat(3,DVpure_nLL_sum,m),3); % sum averaged nLL across subjects
    m = [squeeze(mean(EVpure_nLL_L(:,:,samps-1),3)); squeeze(mean(EVpure_nLL(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    EVpure_nLL_sum = sum(cat(3,EVpure_nLL_sum,m),3);
    
    m = [squeeze(mean(DVpure1_nLL_L(:,:,samps-1),3)); squeeze(mean(DVpure1_nLL(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    DVpure1_nLL_sum = sum(cat(3,DVpure1_nLL_sum,m),3); % sum averaged nLL across subjects
    m = [squeeze(mean(EVpure1_nLL_L(:,:,samps-1),3)); squeeze(mean(EVpure1_nLL(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    EVpure1_nLL_sum = sum(cat(3,EVpure1_nLL_sum,m),3);
    m = [squeeze(mean(Lpure1_nLL_L(:,:,samps-1),3)); squeeze(mean(Lpure1_nLL(:,:,samps-1),3))];  % average nLL across samples and concatenate frequencies
    Lpure1_nLL_sum = sum(cat(3,Lpure1_nLL_sum,m),3);
end

% Computing super BIC difference
ga_superBIC_pure_diff = ((2*DVpure_nLL_sum)+(sum(nPreds).*log(sum(nObs))))-((2*EVpure_nLL_sum)+(sum(nPreds).*log(sum(nObs))));
ga_superBIC_pure_diff1 = ((2*DVpure1_nLL_sum)+(sum(nPreds1).*log(sum(nObs))))-((2*EVpure1_nLL_sum)+(sum(nPreds1).*log(sum(nObs))));
ga_superBIC_purePsiL_diff1 = ((2*DVpure1_nLL_sum)+(sum(nPreds1).*log(sum(nObs))))-((2*Lpure1_nLL_sum)+(sum(nPreds1).*log(sum(nObs))));

% --- %%%%%%%%%% --- %
% --- %% SAVE %% --- %
% --- %%%%%%%%%% --- %
save([savepath,'SuperBIC_TF',surpstr,bstr,'.mat'],'ga*','allfreqs')

