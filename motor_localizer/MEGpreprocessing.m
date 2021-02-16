function [] = MEGpreprocessing(n)

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
subj_in{end+1} = 'TSJ'; sess_in{end+1} = '4'; rec_in{end+1} = '03';

% Pulling 
subject = subj_in{n}; session = sess_in{n}; recording = rec_in{n};

pause(mod(n,30)*15-1)  % implement file-specific pause time to avoid heavy load on cluster

addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
cd /mnt/homes/home024/pmurphy/meg_data/surprise/  % specify path where dataset is located
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts/PreprocScripts/
cfg = [];

subjfile = dir([cd,'/',subject,'-',session,'_Surprise_','*',recording,'.ds']);
cfg.dataset = subjfile(1).name;
%cfg.dataset = 'TFD-1_Surprise_20170315_01.ds';

cfg.trialfun                = 'ft_trialfun_surprise'; % our own function
cfg.trialdef.eventtype      = 'UPPT001';
cfg.trialdef.eventvalue     = [11]; % onset of pre-sequence mask
cfg.trialdef.prestim        = 1.0; % in seconds; for 3 cycles per freq, min(freq)=3Hz, basewin = [-500 -300], this value should be 1.0
cfg.trialdef.postgo         = 0.5; % in seconds; specifies extra time after go cue (marking end of epoch of interest) to include in epoch

art_times = [-0.5 0];  % times within which to check for artifacts; (1)=secs relative to pre-mask; (2)=secs relative to go cue
prestim = cfg.trialdef.prestim;

cfg.subj = subject;
cfg.sess = session;
cfg.rec = recording;

cfg.event = ft_read_event(cfg.dataset);
cfg = ft_definetrial(cfg);
cfg.continuous = 'yes';
data = ft_preprocessing(cfg);

cnt = 1;
remaining_tr = length(data.trial);

% ==================================================================
% 1. REMOVE TRIALS WITH EXCESSIVE HEAD MOTION
% see
% http://www.fieldtriptoolbox.org/example/how_to_incorporate_head_movements_in_meg_analysis
% Head motion in reference to the start of your measurement can be computed
% 
% based on the location of the 3 fiducial coils. Trials with excessive head
% 
% movement can be detected using an outlier test (Grubb?s test). This will
% not 
% detect gradual changes in head position (e.g. sliding down in the
% helmet).
% ==================================================================
% 
% compute head rotation wrt first trial
art_times_t = [ones(length(data.trial),1).*find(data.time{1}>=art_times(1),1,'first') (data.trialinfo(:,5)-prestim+art_times(2))*data.hdr.Fs];  % matrix of per-trial times to apply artifact rejection
cc_rel = computeHeadRotation(data,art_times_t);
% plot the rotation of the head
subplot(2,4,cnt); cnt = cnt + 1;
plot(cc_rel); ylabel('HeadM');
axis tight; box off;

% find outliers
[~, idx] = deleteoutliers(cc_rel);
[t,~]    = ind2sub(size(cc_rel),idx);

% only take those where the deviation is more than 6 mm
t = t(any(abs(cc_rel(t, :)) > 6, 2));

% show those on the plot
hold on;
for thist = 1:length(t),
        plot([t(thist) t(thist)], [max(get(gca, 'ylim')) max(get(gca, 'ylim'))], 'k.');
end

% remove those trials
cfg                     = [];
cfg.trials              = true(1, length(data.trial));
cfg.trials(unique(t))   = false; % remove these trials
data                    = ft_selectdata(cfg, data);
fprintf('removing %d excessive head motion trials \n', length(find(cfg.trials == 0)));

subplot(2,4,cnt); cnt = cnt + 1;
if isempty(t),
    title('No motion'); axis off;
else
    % show head motion without those removed
    cc_rel = computeHeadRotation(data,art_times_t);

    % plo1t the rotation of the head
    plot(cc_rel); ylabel('Motion resid');
    axis tight; box off;
end

% plot a quick power spectrum
% save those cfgs for later plotting
cfgfreq             = [];
cfgfreq.method      = 'mtmfft';
cfgfreq.output      = 'pow';
cfgfreq.taper       = 'hanning';
cfgfreq.channel     = 'MEG';
cfgfreq.foi         = 1:130;
cfgfreq.keeptrials  = 'no';
freq                = ft_freqanalysis(cfgfreq, data);

% plot those data and save for visual inspection
subplot(2,4,cnt); cnt = cnt + 1;
loglog(freq.freq, freq.powspctrm, 'linewidth', 0.1); hold on;
loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
axis tight; axis square; box off;
set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'xticklabel', []);

remaining_tr = [remaining_tr, length(data.trial)];

% ==================================================================
% 2. REMOVE TRIALS WITH EYEBLINKS (only during beginning of trial)
% Bandpass filter the vertical EOG channel between 1-15 Hz and z-transform 
% this filtered time course. Select complete trials that exceed a threshold
% of  
% z =4 (alternatively you can set the z-threshold per data file or per
% subject 
% with the ?interactive? mode in ft_artifact_zvalue function). Reject
% trials 
% that contain blink artifacts before going on to the next step. For
% monitoring 
% purposes, plot the time courses of your trials before and after blink
% rejection.
% ================================================================== 

%ft artifact rejection
cfg                              = [];
cfg.continuous                   = 'no'; % data has been epoched

% channel selection, cutoff and padding
%cfg.artfctdef.zvalue.channel     = {'EOGV'};
%our EOG chanel is just the vertical, and corresponds to EEG057
cfg.artfctdef.zvalue.channel     = {'EEG057'};

% 001, 006, 0012 and 0018 are the vertical and horizontal eog chans
cfg.artfctdef.zvalue.trlpadding  = 0; % padding doesnt work for data thats already on disk
cfg.artfctdef.zvalue.fltpadding  = 0; % 0.2; this crashes the artifact func!
cfg.artfctdef.zvalue.artpadding  = 0.05; % go a bit to the sides of blinks

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';

% set cutoff
cfg.artfctdef.zvalue.cutoff     = 2;
% cfg.artfctdef.zvalue.interactive = 'yes';
[~, artifact_eog]               = ft_artifact_zvalue(cfg, data);

cfg                             = [];
%cfg.artfctdef.reject            = 'complete';  %ch when complete, just the artifacts included in crittoilim will be rejected
cfg.artfctdef.eog.artifact      = artifact_eog;

% ch we , since we need to consider the whole trial
% % reject blinks when they occur between ref and the offset of the stim
% crittoilim = [data.trialinfo(:,2) - data.trialinfo(:,1) - 0.4*data.fsample ...
%         data.trialinfo(:,5) - data.trialinfo(:,1) + 0.8*data.fsample] / data.fsample;
% cfg.artfctdef.crittoilim        = crittoilim;

%TMP skip
%data                            = ft_rejectartifact(cfg, data);

remaining_tr = [remaining_tr, length(data.trial)];  

% ==================================================================
% 3. REMOVE TRIALS WITH SACCADES (only during beginning of trial)
% Remove trials with (horizontal) saccades (EOGH). Use the same settings as
% 
% for the EOGV-based blinks detection. The z-threshold can be set a bit
% higher 
% (z = [4 6]). Reject all trials that contain saccades before going
% further.
% ==================================================================

art_times_t = [ones(length(data.trial),1).*find(data.time{1}>=art_times(1),1,'first') (data.trialinfo(:,5)-prestim+art_times(2))*data.hdr.Fs];  % matrix of per-trial times to apply artifact rejection

% Alternative saccades detection based on eyelink channels
ranges = [5 -5];
screen_x = [0 1920];
screen_y= [0 1080]; 
ch_mapping= ['UADC002' 'UADC003' 'UADC004'];
ppd = estimate_pixels_per_degree();
xcenter = screen_x(2)/2;
ycenter = screen_y(2)/2;
tr_sacc = [];

% Alternative saccades detection with velocity acceleration approach
hz = 1200; % sampling frequency, hdr = ft_read_header(cfg.dataset);
threshold = 30; % taken from Niklas's script
acc_thresh = 2000; % taken from Niklas's script
amp_thresh = 1.5; % our own saccadic amplitude threshold - meant to exclude all but most extreme microsaccades

% For plotting
max_tr = find(data.trialinfo(:,3)==max(data.trialinfo(:,3)));
ypos = nan(length(data.trialinfo),length(data.time{max_tr(1)}));
ypos_times = data.time{max_tr(1)};

for i=1:length(data.trial),
    sacc = 0;
    [x, y, z] = eye_voltage2gaze(data.trial{i}(:,art_times_t(i,1):art_times_t(i,2)), ranges, screen_x, screen_y, ch_mapping);
    ypos(i,1:length(y)) = y;
    sacc = check_saccade(x, y, xcenter, ycenter, ppd);
    try
        if sacc == 0
            % Detect saccades with velocity acceleration approach
            sacc = check_saccade_vel_acc(x, y, hz, threshold, acc_thresh, amp_thresh, ppd);
        end
    catch
        savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed/Data/Saccades/';
        namefile = strcat('ErrorTrial_',subject,'-',session,'_sacc_',recording,'.mat');
        save([savepath,namefile],'i','x', 'y', 'hz', 'threshold', 'acc_thresh', 'amp_thresh', 'ppd','-v7.3'); 
    end
    if sacc > 0
        tr_sacc = [tr_sacc, i];
    end
end
if ~isempty(tr_sacc)
    subplot(2,4,8), hold on;
    plot(ypos_times,ypos(tr_sacc,:)','Color',[1 0 0]);
    plot(ypos_times,ypos(~ismember([1:length(data.trial)],tr_sacc),:)','Color',[0.6 0.6 0.6]);
    xlim([art_times(1) max(art_times_t(:,2))./data.hdr.Fs]);

    fprintf('Saccades detected based on Eyelink channels: %d \n',length(tr_sacc));
    % remove those trials
    cfg                 = [];
    cfg.trials          = true(1, length(data.trial));
    cfg.trials(tr_sacc) = false; % remove these trials
    data                = ft_selectdata(cfg, data);
end

remaining_tr = [remaining_tr, length(data.trial)]; 


% ==================================================================
% 3b. RESAMPLE TO 500Hz
% ==================================================================
cfg3.resample = 'yes';
cfg3.fsample = 1200;
cfg3.resamplefs = 500;
cfg3.detrend = 'yes';  % detrending apparently good for removing edge artifacts during resampling - SHOULDN'T DO THIS IS WANT TO LOOK @ ERFs

data = ft_resampledata(cfg3,data);


% ==================================================================
% 4. REMOVE TRIALS WITH JUMPS
% Compute the power spectrum of all trials and a linear line on the loglog-
% transformed power spectrum. Jumps cause broad range increase in the power
% 
% spectrum so trials containing jumps can be selected by detecting outliers
% 
% in the intercepts of the fitted lines (using Grubb?s test for outliers).
% ==================================================================
%detrend and demean
cfg             = [];
cfg.detrend     = 'yes';
cfg.demean      = 'yes';
cfg.channel     = 'MEG';%ch start to process just the MEG channels

data            = ft_preprocessing(cfg,data);

% get the fourier spectrum per trial and sensor
cfgfreq.keeptrials  = 'yes';
freq                = ft_freqanalysis(cfgfreq, data);

% compute the intercept of the loglog fourier spectrum on each trial
disp('searching for trials with squid jumps...');
intercept       = nan(size(freq.powspctrm, 1), size(freq.powspctrm, 2));
x = [ones(size(freq.freq))' log(freq.freq)'];

for t = 1:size(freq.powspctrm, 1),
    for c = 1:size(freq.powspctrm, 2),
        b = x\log(squeeze(freq.powspctrm(t,c,:)));
        intercept(t,c) = b(1);
    end
end

% detect jumps as outliers
[~, idx] = deleteoutliers(intercept(:));
subplot(2,4,cnt); cnt = cnt + 1;
if isempty(idx),
    fprintf('no squid jump trials found \n');
    title('No jumps'); axis off;
else
    fprintf('removing %d squid jump trials \n', length(unique(t)));
    [t,~] = ind2sub(size(intercept),idx);

    % remove those trials
    cfg                 = [];
    cfg.trials          = true(1, length(data.trial));
    cfg.trials(unique(t)) = false; % remove these trials
    data                = ft_selectdata(cfg, data);

    % plot the spectrum again
    cfgfreq.keeptrials = 'no';
    freq            = ft_freqanalysis(cfgfreq, data);
    loglog(freq.freq, freq.powspctrm, 'linewidth', 0.1); hold on;
    loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
    axis tight; axis square; box off;
    set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'xticklabel', []);
    title(sprintf('%d jumps removed', length(unique(t))));
end
remaining_tr = [remaining_tr, length(data.trial)];

% ==================================================================
% 5. REMOVE LINE NOISE
% ==================================================================
cfg             = [];
cfg.padding     = 9;
cfg.bsfilter    = 'yes';%bandstop filter
cfg.bsfreq      = [49 51; 99 101; 149 151];%bandstop frequeny range,specified as [low high] in hz
data            = ft_preprocessing(cfg, data);

% plot power spectrum
cfgfreq.keeptrials = 'no';
freq            = ft_freqanalysis(cfgfreq, data);
subplot(2,4,cnt); cnt = cnt + 1;
loglog(freq.freq, freq.powspctrm, 'linewidth', 0.5); hold on;
loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
axis tight; %ylim(ylims); %ch todo define ylims
axis square; box off;
title('After bandstop');
set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'xticklabel', []);
remaining_tr = [remaining_tr, length(data.trial)];

% ==================================================================
% 6. REMOVE CARS BASED ON THRESHOLD
% Cars moving past the MEG lab cause big slow signal changes. Trials 
% containing these artifacts can be selected and removed by computing 
% the maximum range of the data for every trial. Trials with a larger 
% range than a threshold (standard = 0.75e-11) can be rejected (the
% standard 
% threshold might be low if you have long trials).
% ==================================================================

disp('Looking for CAR artifacts...');
cfg = [];
cfg.trials = true(1, length(data.trial));
worstChanRange = nan(1, length(data.trial));
for t = 1:length(data.trial),
    % compute the range as the maximum of the peak-to-peak values
    % within each channel
    ptpval = max(data.trial{t}, [], 2) - min(data.trial{t}, [], 2);
    % determine range and index of 'worst' channel
    worstChanRange(t) = max(ptpval);
end

% default range for peak-to-peak
artfctdef.range           = 0.75e-11;

% decide whether to reject this trial
cfg.trials = (worstChanRange < artfctdef.range);
fprintf('removing %d CAR trials \n', length(find(cfg.trials == 0)));
data = ft_selectdata(cfg, data);
remaining_tr = [remaining_tr, length(data.trial)];

% ==================================================================
% 7. REMOVE TRIALS WITH MUSCLE BURSTS BEFORE RESPONSE
% Remove muscle using the same z-value-based approach as for the eye 
% channels. Filter the data between 110-140 Hz and use a z-value threshold
% of 10.
% ==================================================================
cfg                              = [];
cfg.continuous                   = 'no'; % data has been epoched

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = {'MEG'}; % make sure there are no NaNs
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.fltpadding  = 0; % 0.2; - this crashes ft_artifact_zvalue!
cfg.artfctdef.zvalue.artpadding  = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter    = 'yes';
cfg.artfctdef.zvalue.bpfreq      = [110 140];
cfg.artfctdef.zvalue.bpfiltord   = 9;
cfg.artfctdef.zvalue.bpfilttype  = 'but';
cfg.artfctdef.zvalue.hilbert     = 'yes';
cfg.artfctdef.zvalue.boxcar      = 0.2;

% set cutoff
cfg.artfctdef.zvalue.cutoff      = 20;
[~, artifact_muscle]             = ft_artifact_zvalue(cfg, data);

cfg                              = [];
cfg.artfctdef.reject             = 'complete';
cfg.artfctdef.muscle.artifact    = artifact_muscle;

% only remove muscle bursts before the response
crittoilim = [ones(length(data.trial),1).*art_times(1) data.trialinfo(:,5)-prestim+art_times(2)];  % matrix of per-trial times to apply artifact rejection
cfg.artfctdef.crittoilim        = crittoilim;
data                            = ft_rejectartifact(cfg, data);

% plot final power spectrum
freq            = ft_freqanalysis(cfgfreq, data);
subplot(2,4,cnt); 
loglog(freq.freq, freq.powspctrm, 'linewidth', 0.5); hold on;
loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
axis tight; axis square; box off; %ylim(ylims);%ch todo define ylims
set(gca, 'xtick', [10 50 100], 'tickdir', 'out');
remaining_tr = [remaining_tr, length(data.trial)];

data.remaining_tr = remaining_tr;
trl = data.cfg.trl;
%save data file
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed/Data/';
namefile = strcat(subject,'-',session,'_Surprise_Preprocessed_',recording,'.mat');
save([savepath,namefile],'data','remaining_tr','trl','-v7.3');

%save plots
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Preprocessed/Data/Plots/';
namefile = strcat(subject,'-',session,'_Surprise_Preprocessed_Plot_',recording);
print([savepath,namefile],'-dpng');

% % % this is to verify the sacc new approach (to save only info trials
% % % rejected)
% % %save data file
% % savepath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';
% % namefile = strcat(subject,'-',session,'_Surprise_Preprocessed_sacc_',recording,'.mat');
% % save([savepath,namefile],'remaining_tr','trl','-v7.3');
% % 
% % %save plots
% % savepath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Plots/';
% % namefile = strcat(subject,'-',session,'_Surprise_Preprocessed_Plot_sacc_',recording);
% % print([savepath,namefile],'-dpng');
end

function cc_rel = computeHeadRotation(data,ts)

% take only head position channels
cfg         = [];
cfg.channel = {'HLC0011','HLC0012','HLC0013', ...
        'HLC0021','HLC0022','HLC0023', ...
            'HLC0031','HLC0032','HLC0033'};
hpos        = ft_selectdata(cfg, data);

% calculate the mean coil position per trial
coil1 = nan(3, length(hpos.trial));
coil2 = nan(3, length(hpos.trial));
coil3 = nan(3, length(hpos.trial));

for t = 1:length(hpos.trial),
    coil1(:,t) = [mean(hpos.trial{1,t}(1,ts(1):ts(2))); mean(hpos.trial{1,t}(2,ts(1):ts(2))); mean(hpos.trial{1,t}(3,ts(1):ts(2)))];
    coil2(:,t) = [mean(hpos.trial{1,t}(4,ts(1):ts(2))); mean(hpos.trial{1,t}(5,ts(1):ts(2))); mean(hpos.trial{1,t}(6,ts(1):ts(2)))];
    coil3(:,t) = [mean(hpos.trial{1,t}(7,ts(1):ts(2))); mean(hpos.trial{1,t}(8,ts(1):ts(2))); mean(hpos.trial{1,t}(9,ts(1):ts(2)))];
end

% calculate the headposition and orientation per trial (function at
% the
% bottom of this script)
cc     = circumcenter(coil1, coil2, coil3);

% compute relative to the first trial
cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))]';
cc_rel = 1000*cc_rel(:, 1:3); % translation in mm
        
end

function [cc] = circumcenter(coil1,coil2,coil3)

% CIRCUMCENTER determines the position and orientation of the circumcenter
% of the three fiducial markers (MEG headposition coils).
%
% Input: X,y,z-coordinates of the 3 coils [3 X N],[3 X N],[3 X N] where N
% is timesamples/trials.
%
% Output: X,y,z-coordinates of the circumcenter [1-3 X N], and the
% orientations to the x,y,z-axes [4-6 X N].
%
% A. Stolk, 2012

% number of timesamples/trials
N = size(coil1,2);

% x-, y-, and z-coordinates of the circumcenter
% use coordinates relative to point `a' of the triangle
xba = coil2(1,:) - coil1(1,:);
yba = coil2(2,:) - coil1(2,:);
zba = coil2(3,:) - coil1(3,:);
xca = coil3(1,:) - coil1(1,:);
yca = coil3(2,:) - coil1(2,:);
zca = coil3(3,:) - coil1(3,:);

% squares of lengths of the edges incident to `a'
balength = xba .* xba + yba .* yba + zba .* zba;
calength = xca .* xca + yca .* yca + zca .* zca;

% cross product of these edges
xcrossbc = yba .* zca - yca .* zba;
ycrossbc = zba .* xca - zca .* xba;
zcrossbc = xba .* yca - xca .* yba;

% calculate the denominator of the formulae
denominator = 0.5 ./ (xcrossbc .* xcrossbc + ycrossbc .* ycrossbc + zcrossbc .* zcrossbc);

% calculate offset (from `a') of circumcenter
xcirca = ((balength .* yca - calength .* yba) .* zcrossbc - (balength .* zca - calength .* zba) .* ycrossbc) .* denominator;
ycirca = ((balength .* zca - calength .* zba) .* xcrossbc - (balength .* xca - calength .* xba) .* zcrossbc) .* denominator;
zcirca = ((balength .* xca - calength .* xba) .* ycrossbc - (balength .* yca - calength .* yba) .* xcrossbc) .* denominator;

cc(1,:) = xcirca + coil1(1,:);
cc(2,:) = ycirca + coil1(2,:);
cc(3,:) = zcirca + coil1(3,:);

% orientation of the circumcenter with respect to the x-, y-, and z-axis
% coordinates
v = [cc(1,:)', cc(2,:)', cc(3,:)'];
vx = [zeros(1,N)', cc(2,:)', cc(3,:)']; % on the x-axis
vy = [cc(1,:)', zeros(1,N)', cc(3,:)']; % on the y-axis
vz = [cc(1,:)', cc(2,:)', zeros(1,N)']; % on the z-axis

for j = 1:N
    % find the angles of two vectors opposing the axes
    thetax(j) = acos(dot(v(j,:),vx(j,:))/(norm(v(j,:))*norm(vx(j,:))));
    thetay(j) = acos(dot(v(j,:),vy(j,:))/(norm(v(j,:))*norm(vy(j,:))));
    thetaz(j) = acos(dot(v(j,:),vz(j,:))/(norm(v(j,:))*norm(vz(j,:))));

    % convert to degrees
    cc(4,j) = (thetax(j) * (180/pi));
    cc(5,j) = (thetay(j) * (180/pi));
    cc(6,j) = (thetaz(j) * (180/pi));
end

end

function sacc = check_saccade(xgaze, ygaze, xcenter, ycenter, ppd)
    sacc = 0;
    x = (xgaze-xcenter)/ppd;
    y = (ygaze-ycenter)/ppd;
    d = (x.^2 + y.^2).^.5;
    a=d(2:length(d));
    if any(a>4)
            sacc = 1;
    end
end

function [sacc, sacc_times] = check_saccade_vel_acc(xgaze, ygaze, Hz, threshold, acc_thresh, amp_thresh, ppd)
    
    sacc = 0;  % initializing binary output flag that indicates whether this trial contains a saccade
    sacc_times = [];
    
    % get x and  y in degrees
    x = xgaze/ppd;
    y = ygaze/ppd;
    [velocity, acceleration] = get_velocity (x, y, Hz);
    saccades = double(velocity > threshold);  % getting vector of samples that violate velocity threshold
    
    borders = diff(saccades);   % start points of candidate saccades will be 1, end points will be -1, all others 0
    if saccades(1)>threshold, borders = [1 borders]; else borders = [0 borders]; end  % in case first sample violates threshold
    if saccades(end)>threshold, borders(end+1) = -1; end  % in case last sample violates threshold
    
    window_size = 3;
    win = ones(1,window_size)/double(window_size);
    x = conv(x, win, 'same'); y = conv(y, win, 'same');   % lightly smoothing gaze time series before applying amplitude threshold
    
    starts = find(borders==1); ends = find(borders==-1)-1;  % getting all start and end points of candidate saccades
    if length(starts)>length(ends), starts(end) = []; end   % if last saccade occurs right at trial end (and so doesn't have corresponding ends), discard
    for i = 1:length(starts)  % looping through all candidate saccades and only accepting if they also violate acceleration/amplitude thresholds
        if ~isempty(find(acceleration(starts(i):ends(i))>acc_thresh))  % applying acceleration threshold
            % applying amplitude threshold based on average gaze position over 3 samples before vs after saccade
            if i>1 && starts(i)-3>0  % a few exceptions to try to make sure preceding 3 samples are useable
                t1 = max([starts(i)-3,ends(i-1)]):starts(i)-1;
            elseif starts(i)<=2
                t1 = 1;
            else t1 = starts(i)-1;
            end
            if i<length(starts) && ends(i)+3<=length(x)  % a few exceptions to try to make sure following 3 samples are useable
                t2 = ends(i)+1:min([ends(i)+3,starts(i+1)]);
            elseif ends(i)>=length(x)-1
                t2 = length(x);
            else t2 = ends(i)+1;
            end
            p1 = [mean(x(t1)) mean(y(t1))];  % x,y coords of pre-saccade gaze position (in d.v.a., always positive)
            p2 = [mean(x(t2)) mean(y(t2))];  % x,y coords of post-saccade gaze position (in d.v.a., always positive)
            if ((p1(1)-p2(1)).^2 + (p1(2)-p2(2)).^2).^.5 > amp_thresh;  % applying amplitude threshold
                sacc = 1;
                sacc_times(end+1,1:2) = [starts(i) ends(i)];
            end
        end
    end
        
end

function [velocity,acceleration] = get_velocity (x, y, Hz)

    % Compute velocity and acceleration of eye movements
    % Based on Niklas py script "The function asumes that the values in x,y are
    % sampled continuously at a rate specified by Hz"
    velocity_window_size = 3;
    Hz = double(Hz);
    distance = (diff(x).^2 + diff(y).^2).^.5;
    distance = [distance(1) distance];
    win = ones(1,velocity_window_size)/double(velocity_window_size);
    velocity = conv(distance, win, 'same');
    velocity = velocity / (velocity_window_size/Hz);
    acceleration = diff(velocity) / (1/Hz);
    acceleration = abs([acceleration(1) acceleration]);

end
