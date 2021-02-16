function [] = MEGpreprocessingMotor(subject, session, recording)
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
% %addpath '/Users/caro/Documents/UKE/Matlab/Toolboxes/fieldtrip-20170322'
ft_defaults
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts/PreprocScripts/
cd /mnt/homes/home024/pmurphy/meg_data/surprise/  % specify path where dataset is located
dir_behav = '/mnt/homes/home024/pmurphy/Surprise_accumulation/DataSensLoc/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/PreprocessedMotor/Data/';

subjfile = dir([cd,'/',subject,'-',session,'_Surprise_','*',recording,'.ds']);

% ==================================================================
% 1. READ EVENTS FROM ENTIRE RECORDING & DEFINE TRIALS
% ==================================================================
cfgtrl = [];
cfgtrl.dataset = subjfile(1).name;
cfgtrl.subj = subject; cfgtrl.sess = session; cfgtrl.rec = recording;

cfgtrl.trialfun                = 'ft_trialfun_surprise_motor'; % our own function
cfgtrl.trialdef.eventtype      = 'UPPT001';
cfgtrl.trialdef.eventvalue     = 36;  % onset of response direction cue   %%% events: 5 = block start, 16 = rest on, 26 = fixation on, 46 = go cue on, 56/57/58 = resp L/R/bad, 6 = block end
cfgtrl.trialdef.prestim        = 0.5; % in seconds
cfgtrl.trialdef.poststim       = 1.6; % in seconds

art_times = [-0.4 1.3];  % times within which to check for artifacts (secs relative to response direction cue)
art_times_blink = [-0.3 1.3];

cfgtrl.event = ft_read_event(subjfile(1).name);
cfgtrl = ft_definetrial(cfgtrl);

%%%%% need to at some point incorporate exception that will appropriately deal with QNV-3 (where triggers need to be recovered from EL file)
if sum(ismember([cfgtrl.event.value],36))>0 % only proceeding if localizer block for this subject/session combination is contained in this dataset
        
    % ==================================================================
    % 2. HIGH-PASS & LINE NOISE FILTERING, & JUMP ID
    % ==================================================================
    % Loop through blocks, load data, filter, epoch, and concatenate
    hdr = ft_read_header(cfgtrl.dataset);
    
    jumpts = [];
    filtpad = 5;  % padding around block bounds for filtering (in seconds)
    
    cfg = [];
    cfg.trl = [max([1 cfgtrl.trl(1,1)+cfgtrl.trl(1,3)-(filtpad*hdr.Fs)]) min([cfgtrl.trl(end,2)+(filtpad*hdr.Fs) cfgtrl.event(end).sample]) 0];  % block bounds
    cfg.dataset = subjfile(1).name;
    cfg.event = cfgtrl.event;
    cfg.channel = {'M*','HLC*','UADC*','UPPT*','EEG057'};   % keeping only useful channels
    cfg.continuous = 'yes';
    data = ft_preprocessing(cfg);
    
    % Pull only MEG channels for filtering
    cfg=[];
    cfg.channel = {'M*'};
    datatemp = ft_selectdata(cfg,data);
    
    %%% ID and store time stamps with putative JUMP artifacts before filtering
    disp('Identifying squid jump artifacts...')
    intercept = [];
    winsize = 2;  % size, in s, of sliding window within which to test for jumps
    slidesize = 1;  % gap, in s, to jump sliding window along timeseries
    rejextend = 2;  % extent, in s, to which jump time window should be extended (forwards and backwards) for artifact rejection
    
    winstarts = (1:data.fsample*slidesize:(length(data.time{1})-data.fsample*winsize))';
    for w = 1:length(winstarts)
        pow = abs(fft(datatemp.trial{1}(:,winstarts(w):(winstarts(w)+data.fsample*winsize))'))';  % calculate power spectra
        freqs = [0:size(pow,2)-1]*data.fsample/size(pow,2);
        pow = pow(:,freqs>=1 & freqs<=100); freqs = freqs(:,freqs>=1 & freqs<=100);  % pull only power b/w 1-100Hz
        x = [ones(length(freqs),1) log(freqs)'];
        for c = 1:size(datatemp.trial{1},1)
            beta = x\log(pow(c,:)');  % fit line to log-log power spectrum for this trial/sensor
            intercept(w,c) = beta(1);
        end
    end
    [~, idx] = deleteoutliers(intercept(:));
    if isempty(idx),
        fprintf('no squid jump trials found \n');
    else
        [t,~] = ind2sub(size(intercept),idx);
        jumpts = [winstarts(t)-rejextend*data.fsample winstarts(t)+(winsize+rejextend)*data.fsample] + data.sampleinfo(1);
        jumpts(jumpts(:,1)<1,1) = 1; jumpts(jumpts(:,1)>cfgtrl.event(end).sample,1) = cfgtrl.event(end).sample;  % making sure ID'd periods + rejextend lie within data bounds
    end
    
    % Filter
    cfg=[];
    cfg.demean = 'yes';
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 0.5;
    cfg.hpfilttype = 'firws';
    % cfg.hpfiltord = 8;
    
    cfg.bsfilter    = 'yes';%bandstop filter
    cfg.bsfreq      = [49 51; 99 101; 149 151];%bandstop frequeny range,specified as [low high] in hz
    
    datatemp = ft_preprocessing(cfg,datatemp);
    
    % Add filtered MEG data back in
    data.trial{1}(ismember(data.label,ft_channelselection({'M*'}, data.label)),:) = datatemp.trial{1};
    clear datatemp
    
    % Extract trial epochs
    cfg=[];
    cfg.trl = cfgtrl.trl;
    data = ft_redefinetrial(cfg,data);
    
    % Set up for artifact identification routines
    cnt = 1;
    remaining_tr = length(data.trial);
    data.hdr = hdr;    
    
    % =================================================================
    % 3. REMOVE TRIALS WITH PRE-ID'D JUMP ARTIFACTS
    % =================================================================
    % colesce overlapping jump periods
    if size(jumpts,1)>1
        cmsmp = jumpts(1,:);
        for w = 1:size(jumpts,1)-1
            if jumpts(w+1,1) - cmsmp(end,2) < 1
                cmsmp(end,2) = jumpts(w+1,2);
            else
                cmsmp(end+1,:) = jumpts(w+1,:);
            end
        end
        jumpts = cmsmp; clear cmsmp
    end
    
    % loop through trials and throw out if they overlap with car/jump periods
    tr_bad = [];
    for t = 1:length(data.trial)
        j = 1;
        while j<=size(jumpts,1)
            if sum(ismember(jumpts(j,1):jumpts(j,2),data.sampleinfo(t,1):data.sampleinfo(t,2)))>0 || sum(ismember(data.sampleinfo(t,1):data.sampleinfo(t,2),jumpts(j,1):jumpts(j,2)))>0
                tr_bad(end+1) = t;
                break
            else j = j+1;
            end
        end
    end
    
    % Plot & remove bad trials
    subplot(2,3,cnt); cnt = cnt + 1;
    if isempty(tr_bad),
        fprintf('No jump trials found \n');
        title('No jumps'); axis off;
    else
        fprintf('Removing %d squid jump trials \n', length(tr_bad));
        
        % plot trace of most outlying channel (max-min) for each artifactual trial
        hold on, maxvals=[];
        for t = 1:length(tr_bad)
            cdata = data.trial{tr_bad(t)}(ismember(data.label,ft_channelselection({'M*'}, data.label)),:);
            maxvals = max(cdata,[],2)-min(cdata,[],2);
            plot(data.time{tr_bad(t)},cdata(find(maxvals==max(maxvals),1,'first'),:),'Color',rand(1,3))
        end
        axis tight; axis square; box off;
        set(gca, 'tickdir', 'out');
        title(sprintf('%d jumps removed', length(tr_bad)));
        
        % remove
        cfg                 = [];
        cfg.trials          = true(1, length(data.trial));
        cfg.trials(tr_bad)  = false; % remove these trials
        data                = ft_selectdata(cfg, data);
    end
    
    remaining_tr = [remaining_tr, length(data.trial)];
    
    
    % ==================================================================
    % 4. REMOVE TRIALS WITH EXCESSIVE HEAD MOTION
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
    art_times_t = [ones(length(data.trial),1).*find(data.time{1}>=art_times(1),1,'first') ones(length(data.trial),1).*find(data.time{1}<=art_times(2),1,'last')];  % matrix of per-trial times to apply artifact rejection
    cc_rel = computeHeadRotation(data,art_times_t);
    % plot the rotation of the head
    subplot(2,3,cnt); cnt = cnt + 1;
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
    
    subplot(2,3,cnt); cnt = cnt + 1;
    if isempty(t),
        title('No motion'); axis off;
    else
        % show head motion without those removed
        cc_rel = computeHeadRotation(data,art_times_t);
        
        % plot the rotation of the head
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
    subplot(2,3,cnt); cnt = cnt + 1;
    loglog(freq.freq, freq.powspctrm, 'linewidth', 0.1); hold on;
    loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
    axis tight; axis square; box off;
    set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'xticklabel', []);
    
    remaining_tr = [remaining_tr, length(data.trial)];
    
    if strcmp(subject,'QNV') && strcmp(session,'3')   % use EOG for blink detection here, since eyelink is missing
    % ==================================================================
    % 5. REMOVE TRIALS WITH EYEBLINKS (only during beginning of trial)
    % Bandpass filter the vertical EOG channel between 1-15 Hz and z-transform
    % this filtered time course. Select complete trials that exceed a threshold
    % of z =4 (alternatively you can set the z-threshold per data file or per subject
    % with the ?interactive? mode in ft_artifact_zvalue function). Reject trials
    % that contain blink artifacts before going on to the next step. For monitoring
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
    cfg.artfctdef.zvalue.cutoff     = 2.5;
    % cfg.artfctdef.zvalue.interactive = 'yes';
    [~, artifact_eog]               = ft_artifact_zvalue(cfg, data);
    
    cfg                             = [];
    %cfg.artfctdef.reject            = 'complete';  %ch when complete, just the artifacts included in crittoilim will be rejected
    cfg.artfctdef.eog.artifact      = artifact_eog;
    
    crittoilim = [ones(length(data.trial),1).*art_times_blink(1) ones(length(data.trial),1).*art_times_blink(2)];  % matrix of per-trial times to apply artifact rejection
    cfg.artfctdef.crittoilim        = crittoilim;
    
    data                            = ft_rejectartifact(cfg, data);
    
    remaining_tr = [remaining_tr, length(data.trial)];
    
    else
    % ==================================================================
    % 6. REMOVE TRIALS WITH SACCADES DETECTED VIA EYELINK DATA
    % ==================================================================    
    % Alternative saccades detection based on eyelink channels
    art_times_t = [ones(length(data.trial),1).*find(data.time{1}>=art_times_blink(1),1,'first') ones(length(data.trial),1).*find(data.time{1}<=art_times_blink(2),1,'last')];  % matrix of per-trial times to apply artifact rejection
    
    % Alternative saccades detection based on eyelink channels
    ranges = [5 -5];
    screen_x = [0 1920];
    screen_y= [0 1080];
    ch_mapping= [find(ismember(data.label,{'UADC002'})) find(ismember(data.label,{'UADC004'})) find(ismember(data.label,{'UADC003'}))];   % [x-gaze y-gaze pupil]
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
    max_tr = 1;
    ypos = nan(size(data.trialinfo,1),length(data.time{max_tr(1)}));
    ypos_times = data.time{max_tr(1)};
    
    for i=1:length(data.trial)
        sacc = 0;
        [x, y, ~] = eye_voltage2gaze(data.trial{i}(:,art_times_t(i,1):art_times_t(i,2)), ranges, screen_x, screen_y, ch_mapping);
        ypos(i,1:length(y)) = y;
        sacc = check_saccade(x, y, xcenter, ycenter, ppd);
        try
            if sacc == 0
                % Detect saccades with velocity acceleration approach
                sacc = check_saccade_vel_acc(x, y, hz, threshold, acc_thresh, amp_thresh, ppd);
            end
        catch
            savepath2 = [savepath,'Saccades/'];
            namefile = strcat('ErrorTrial_',subject,'-',session,'_sacc_',recording,'.mat');
            save([savepath2,namefile],'i','x', 'y', 'hz', 'threshold', 'acc_thresh', 'amp_thresh', 'ppd','-v7.3');
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
    end
    
    % ==================================================================
    % 7. REMOVE CARS BASED ON THRESHOLD
    % Cars moving past the MEG lab cause big slow signal changes. Trials
    % containing these artifacts can be selected and removed by computing
    % the maximum range of the data for every trial. Trials with a larger
    % range than a threshold (standard = 0.75e-11) can be rejected (the
    % standard threshold might be low if you have long trials).
    % ==================================================================
    % demean and keep only MEG channels
    cfg             = [];
    cfg.detrend     = 'no';    %%%%%% PM: set this to NO for all ERF analyses
    cfg.demean      = 'yes';
    cfg.channel     = 'MEG';%ch start to process just the MEG channels
    
    data            = ft_preprocessing(cfg,data);
    
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
    % 8. REMOVE TRIALS WITH MUSCLE BURSTS
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
    crittoilim = [ones(length(data.trial),1).*art_times(1) ones(length(data.trial),1).*art_times(2)];  % matrix of per-trial times to apply artifact rejection
    cfg.artfctdef.crittoilim        = crittoilim;
    data                            = ft_rejectartifact(cfg, data);
    
    % ==================================================================
    % 9. RESAMPLE TO 500Hz
    % ==================================================================
    trl = data.cfg.trl;
    
    cfg3=[];
    cfg3.resample = 'yes';
    cfg3.fsample = 1200;
    cfg3.resamplefs = 500;
    cfg3.detrend = 'yes';  % detrending apparently good for removing edge artifacts during resampling - SHOULDN'T DO THIS IF WANT TO LOOK @ ERFs
    
    data = ft_resampledata(cfg3,data);
    
    % ==================================================================
    % 10. PLOT AND SAVE
    % ==================================================================
    freq            = ft_freqanalysis(cfgfreq, data);
    subplot(2,3,cnt);
    loglog(freq.freq, freq.powspctrm, 'linewidth', 0.5); hold on;
    loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
    axis tight; axis square; box off; %ylim(ylims);%ch todo define ylims
    set(gca, 'xtick', [10 50 100], 'tickdir', 'out');
    remaining_tr = [remaining_tr, length(data.trial)];
        
    data.remaining_tr = remaining_tr;
    %save data file
    namefile = strcat(subject,'-',session,'_Surprise_PreprocessedMotor.mat');
    save([savepath,namefile],'data','remaining_tr','trl','-v7.3');
    
    %save plots
    savepath2 = [savepath,'Plots/'];
    namefile = strcat(subject,'-',session,'_Surprise_PreprocessedMotor_Plot');
    print([savepath2,namefile],'-dpng');
    
end
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
