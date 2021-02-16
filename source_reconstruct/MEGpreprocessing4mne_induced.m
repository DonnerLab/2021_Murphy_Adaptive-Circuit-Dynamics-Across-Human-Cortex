function [] = MEGpreprocessing4mne_induced(n)

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
subj_in{end+1} = 'QNV';
sess_in{end+1} = '4';
rec_in{end+1} = '01';

subj_in{end+1} = 'QNV';
sess_in{end+1} = '4';
rec_in{end+1} = '02';

%Subjects 3 reccordings
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

% Pulling file ID and pausing
subject = subj_in{n}; session = sess_in{n}; recording = rec_in{n};

% pause(mod(n,30)*15-1)  % implement file-specific pause time to avoid heavy load on cluster

addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
cd /mnt/homes/home024/pmurphy/meg_data/surprise/  % specify path where dataset is located
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts/PreprocScripts/

subjfile = dir([cd,'/',subject,'-',session,'_Surprise_','*',recording,'.ds']);

art_times = [-0.5 0];  % times within which to check for artifacts; (1)=secs relative to pre-mask; (2)=secs relative to go cue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD PREVIOUSLY-PROCESSED FILE AND STORE IDs OF CLEAN TRIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed4mne/';
namefile = strcat(subject,'-',session,'_',recording,'_preproc4mne.mat');
load([loadpath,namefile]);
trialIDs = data.trialinfo(:,end);
times = data.time;
sampleinfo = data.sampleinfo;
clear data

% ==================================================================
% 1. READ EVENTS FROM ENTIRE RECORDING & DEFINE TRIALS
% ================================================================== 
cfgtrl = [];
cfgtrl.dataset = subjfile(1).name;
cfgtrl.subj = subject; cfgtrl.sess = session; cfgtrl.rec = recording;

cfgtrl.trialfun                = 'ft_trialfun_surprise'; % our own function
cfgtrl.trialdef.eventtype      = 'UPPT001';
cfgtrl.trialdef.eventvalue     = [11]; % onset of pre-sequence mask
cfgtrl.trialdef.prestim        = 1.0; % in seconds; for 3 cycles per freq, min(freq)=3Hz, basewin = [-500 -300], this value should be 1.0
cfgtrl.trialdef.postgo         = 0.5; % in seconds; specifies extra time after go cue (marking end of epoch of interest) to include in epoch

cfgtrl.event = ft_read_event(cfgtrl.dataset);
cfgtrl = ft_definetrial(cfgtrl);

prestim = cfgtrl.trialdef.prestim;

% ==================================================================
% 2. PER-BLOCK HIGH-PASS & LINE NOISE FILTERING
% ==================================================================
% Loop through blocks, load data, filter, epoch, and concatenate
hdr = ft_read_header(cfgtrl.dataset);
blocks = unique(cfgtrl.trl(:,4));
filtpad = 10;  % padding around block bounds for filtering (in seconds)
for b = 1:length(blocks)
    fprintf('\n\nPROCESSING BLOCK %d of %d...\n\n',b,length(blocks))
    ctrl = cfgtrl.trl(cfgtrl.trl(:,4)==blocks(b),:); % pull trial info for only this block
    
    cfg = [];
    cfg.trl = [max([1 ctrl(1,1)+ctrl(1,3)-(filtpad*hdr.Fs)]) min([ctrl(end,2)+(filtpad*hdr.Fs) cfgtrl.event(end).sample]) 0];  % block bounds
    cfg.dataset = subjfile(1).name;
    cfg.event = cfgtrl.event;
    cfg.channel = {'M*','HLC*','UADC*','UPPT*','EEG057'};   % keeping only useful channels
    cfg.continuous = 'yes';
    dataB = ft_preprocessing(cfg);
    
    % Store unfiltered EyeLink and headloc channels
    ELchans = ft_channelselection({'UADC*'}, dataB.label);
    ELdata = dataB.trial{1}(ismember(dataB.label,ELchans),:);
    
    HLchans = ft_channelselection({'HLC*'}, dataB.label);
    HLdata = dataB.trial{1}(ismember(dataB.label,HLchans),:);
    
    % Rectifying EOG for QNV-3
    dataB.trial{1}(ismember(dataB.label,'EEG057'),:) = dataB.trial{1}(ismember(dataB.label,'EEG057'),:).*-1;
    
    % Filter
    cfg=[];
    cfg.demean = 'yes';
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 0.5;
    cfg.hpfilttype = 'firws';
    % cfg.hpfiltord = 8;
    
    cfg.bsfilter    = 'yes';%bandstop filter
    cfg.bsfreq      = [49 51; 99 101; 149 151];%bandstop frequeny range,specified as [low high] in hz
    
    dataB = ft_preprocessing(cfg,dataB);
    
    % Add EL and headloc data back in
    dataB.trial{1}(ismember(dataB.label,ELchans),:) = ELdata;
    dataB.trial{1}(ismember(dataB.label,HLchans),:) = HLdata;
    clear ELdata HLdata
    
    % Extract trial epochs
    cfg=[];
    cfg.trl = ctrl;
    dataB = ft_redefinetrial(cfg,dataB);
    
    % Append blocks
    if b == 1
        data = dataB;
    else
        data = ft_appenddata([],data,dataB);
    end
    clear dataB
end

cnt = 1;
remaining_tr = length(data.trial);
data.hdr = hdr;

% ==================================================================
% 3. RESAMPLE TO 400Hz
% ==================================================================
cfg3.resample = 'yes';
cfg3.fsample = 1200;
cfg3.resamplefs = 400;
cfg3.detrend = 'yes';  % detrending apparently good for removing edge artifacts during resampling - SHOULDN'T DO THIS IF WANT TO LOOK @ ERFs

data = ft_resampledata(cfg3,data);

% ==================================================================
% 4. ADD EXTRA COLUMN INTO TRIALINFO MATRIX WITH UNIQUE TRIAL ID
% ==================================================================
subjnum = find(strcmp(allsubj,subject));
sessnum = str2double(session);

data.trialinfo(:,end+1) = [ones(size(data.trialinfo,1),1).*subjnum + ones(size(data.trialinfo,1),1).*sessnum + data.trialinfo(:,1).*100 + data.trialinfo(:,2)];

% ==================================================================
% 5. KEEP ONLY CLEAN TRIALS FROM ORIGINAL PREPROCESSING
% ==================================================================
assert(length(unique(data.trialinfo(:,end)))==length(data.trialinfo(:,end)),'Replicates of trial ID!!!')  % testing whether all trial IDs are unique

cfg = [];
cfg.trials = ismember(data.trialinfo(:,end),trialIDs);
fprintf('Keeping only full-length trials...\n');
data = ft_selectdata(cfg, data);

assert(length(find(data.trialinfo(:,12)<12))==0,'Trials included that are not full length!!!')  % testing whether all selected trials are full-length

% ==================================================================
% 6. RESTRUCTURE DATA STRUCTURE FOR MNE
% ==================================================================
min_tr = find(data.trialinfo(:,3)==min(data.trialinfo(:,3)),1,'first');

assert(length(times)==length(data.time{min_tr}),'Data of different length than original preprocessed!!!')

data_mne = struct;
data_mne.time = data.time{min_tr};
data_mne.label = data.label;
data_mne.dimord = 'rpt_chan_time';
data_mne.trialinfo = data.trialinfo;
data_mne.sampleinfo = sampleinfo;
data_mne.grad = data.grad;

trial_mne = [];
for t = 1:length(data.trial)
    trial_mne(t,:,:) = data.trial{t}(:,1:length(data_mne.time));
end
clear data

data_mne.trial = trial_mne;

data = data_mne;
clear data_mne

% ==================================================================
% 7. BASELINING SINGLE-TRIAL DATA
% ==================================================================
basewin = [-0.25 0];  % matches that used by pymeg
MEGchans = ft_channelselection({'M*'}, data.label);
for c = find(ismember(data.label,MEGchans))'
    data.trial(:,c,:) = squeeze(data.trial(:,c,:))-repmat(squeeze(mean(data.trial(:,c,data.time>=basewin(1) & data.time<=basewin(2)),3)),1,size(data.trial,3));
end

% ==================================================================
% 8. REGRESS OUT CHANNEL-WISE TRIAL-AVERAGED RESPONSE
% ==================================================================
model_deets=[];
for c = find(ismember(data.label,MEGchans))'
    st = reshape(squeeze(data.trial(:,c,:))',[size(data.trial,1)*size(data.trial,3) 1]);  % concatenated single-trial responses
    av = repmat(squeeze(mean(data.trial(:,c,:),1)),size(data.trial,1),1);  % vector of trial-averaged response replicated for each trial
    m = regstats(st,av,'linear',{'tstat','rsquare','r'});
    model_deets(end+1,1:2) = [m.rsquare m.tstat.beta(2)];
    data.trial(:,c,:) = reshape(m.r,[size(data.trial,3),size(data.trial,1)])';
end

% ==================================================================
% 9. SAVE DATA
% ==================================================================
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed4mne/';
namefile = strcat(subject,'-',session,'_',recording,'_preproc4mne_induced.mat');
save([savepath,namefile],'data','-v7.3');
save([savepath,strcat(subject,'-',session,'_',recording,'_induced_model_deets.mat')],'model_deets','-v7.3')

end


