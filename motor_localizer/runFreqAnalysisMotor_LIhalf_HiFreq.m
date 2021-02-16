function [] = runFreqAnalysisMotor_LIhalf_HiFreq(subject)
  % Applies two baselining operations (pre-trial and peri-sample) to MEG
  % TF data and regresses these data (channel*time) onto model-based
  % variables of interest

% allsubj = {'DCB' 'DHB' 'EMB' 'EXF' 'HBC' 'JTB' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB'};
% subject = allsubj{n};
  
basewin = [-0.3 -0.2];  % baseline window relative to pre-mask onset (s)
basetype = 'dB_common';  % options: 'dB_st', 'dB_common', 'pc_st' or 'pc_common'

trlwin = [-0.3 1.3];  % window for full-trial-wise analyses

RTcutoff = 1.3;  % discard trials with RTs quicker than this value (s) - value is relative to direction cue, so go cue will appear always at 1.3s

% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
megpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/PreprocessedMotor/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/TF_Motor/';

if ~exist([savepath,basetype,filesep],'dir'), mkdir([savepath,basetype,filesep]), end  % making save directory if doesn't exist
savepath = [savepath,basetype,filesep];

subjfiles = dir([megpath,subject,'-*TF_HiFreq.mat']);  % pull all meg filenames for this subject

trl_data=[];
choices_full=[]; RT_full=[]; sess_full=[];

% ==================================================================
% TRIAL-WISE ANALYSIS
% ==================================================================
fprintf('Beginning analysis...\n')
for f = 1:length(subjfiles)

    load([megpath,subjfiles(f).name])  % load meg data
    sess = str2double(subjfiles(f).name(5));
    freq.time = round(freq.time,2);  % rounding time vector to nearest 2nd decimal - otherwise slight inaccuracies can lead to bad timing later
    
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
    fprintf('Concatenating trial-wise data segments for %s...\n',subjfiles(f).name)
    if f==1,
        freqs = freq.freq;
        trltimes = freq.time(freq.time>=trlwin(1) & freq.time<=trlwin(2));
    end
    trl_data = cat(1,trl_data,freq.powspctrm(:,:,:,freq.time>=trlwin(1) & freq.time<=trlwin(2)));
    
    % ==================================================================
    % CONCATENATE BEHAVIOURAL DATA
    % ==================================================================    
    choices_full = [choices_full; freq.resps-56];  % 0=left response, 1=right response
    RT_full = [RT_full; freq.RT];
    sess_full = [sess_full; ones(length(freq.resps),1).*sess];
    
    % ==================================================================
    % STORE FT STRUCTURES FOR LATER PLOTTING
    % ==================================================================
    if f==1,
        grad_full = freq.grad;
        
        grad = freq.grad;
        grad.chanpos = grad.chanpos(Rchans,:);
        grad.chanori = grad.chanori(Rchans,:);
        grad.chantype = grad.chantype(Rchans);
        grad.chanunit = grad.chanunit(Rchans,:);
        grad.label = grad.label(Rchans);
        
        cfg = freq.cfg;
    end
end
clear freq

% ==================================================================
% THROW OUT PREMATURE RESPONSE TRIALS
% ==================================================================
trl_data = trl_data(RT_full>=RTcutoff,:,:,:);
choices_full = choices_full(RT_full>=RTcutoff);
sess_full = sess_full(RT_full>=RTcutoff);
RT_full = RT_full(RT_full>=RTcutoff);

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
trl_avg_full = squeeze(mean(trl_data(choices_full==1,:,:,:),1)) - squeeze(mean(trl_data(choices_full==0,:,:,:),1));

for c = 1:length(Rchans)
    trl_data(:,c,:,:) = trl_data(:,Rchans(c),:,:) - trl_data(:,Lchans(c),:,:);  % store full-trial, half-scalp LI data - trials*chans*freqs*times (overwriting to save memory - works b/c L/R chans are well separated)
end
trl_data(:,length(Rchans)+1:end,:,:) = [];

trl_avg = squeeze(mean(cat(1,trl_data(choices_full==1,:,:,:),trl_data(choices_full==0,:,:,:).*-1),1));

% ==================================================================
% MAKE CATEGORICAL SESSION REGRESSORS
% ==================================================================
sessions = unique(sess_full);
sess_r = zeros(length(sess_full),length(sessions)-1);
for s = 1:length(sessions)-1
    sess_r(sess_full==sessions(s+1),s) = ones(length(find(sess_full==sessions(s+1))),1);
end

% ==================================================================
% RUN SINGLE-TRIAL REGRESSION
% ==================================================================
% Full-trial regression, time*freq*space
TrespF=[];
fprintf('Running regressions of full-trial power (time*freq*space) onto response direction...\n Sensor')
for c = 1:size(trl_data,2)  % looping through channels
    if mod(c,10)==0, fprintf('%d,',c), end
    for f = 1:size(trl_data,3)  % looping through freqs
        for t = 1:size(trl_data,4)  % looping through time-points
            m = regstats(trl_data(:,c,f,t),[choices_full sess_r],'linear',{'tstat'});  % signed posterior belief
            TrespF(c,f,t) = m.tstat.t(2);
        end
    end
end
fprintf(' Done.\n')

% ==================================================================
% SAVE RESULTS AND CLEAN UP
% ==================================================================
save([savepath,subject,'_output_LIhalf_HiFreq.mat'],'trltimes','freqs','grad','grad_full','cfg','trl_avg','trl_avg_full','TrespF')





