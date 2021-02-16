function [newpupil, newXgaze, newYgaze, newblinksmp, nanIdx] = blink_interpolatePM(dat, inFile, plotme)
% interpolates blinks and missing data
% Anne Urai, 2016
% Peter Murphy, 2016 - added adaptive interpolation windows 

% interpolation parameters
padding1         = [-0.150 0.150]; % padding before/after EL-defined blinks for initial rough interpolation
wfilt            = 44;             % width (in samples) of Hanning window used for lo-pass filtering

diffthresh       = 1.0;            % z-scored PUPIL derivative threshold that samples before/after a bad period must consistently satisfy   (default = 0.75 for quick blinkers, 0.55 for slow)
gazethresh       = 6.0;            % z-scored GAZE derivative threshold that samples before/after a bad period must consistently satisfy   (default = 6.0 - gaze tends to be noisier/higher freq than pupil)
t_pre            = 0.05;           % window (s) before a bad period in which all samples must satisfy diffthresh   (default = 0.05)
t_post           = 0.06;           % window (s) after a bad period in which all samples must satisfy diffthresh   (default = 0.06)
tmax             = [0.15 0.25];     % max time (s) [start end] points of interpolation window can go from identified bad samples   (default = 0.25 for quick blinkers, 0.5 for slow)
diffthresh2      = 3.5;            % z-scored derivative threshold for identification of bad periods not covered by Eyelink   (default = 3.5)

coalesce1        = 0.250;          % merge 2 blinks into 1 if they are below this distance (in s) apart (default = 0.250)
coalesce2        = 0.500;          % merge 2 bad periods (IDd from step 2) into 1 if they are below this distance (in s) apart (default = 0.500)

% make copies of some stuff
blinksmp = dat.blinksmp;
pupilcopy = dat.pupil;
Xgazecopy = dat.Xgaze;
Ygazecopy = dat.Ygaze;
nanIdx = [];

% plot if unspecified
if ~exist('plotme', 'var'); plotme = true; end

% ====================================================== %
% STEP 1: INTERPOLATE EL-DEFINED BLINKS USING WIDE FIXED WINDOW
% ====================================================== %
% plot raw time-series
if plotme,
    figure('Position', [80, 100, 1500, 840]);
    sp1 = subplot(511); plot(dat.times,dat.pupil, 'color', [0.5 0.5 0.5]);
    axis tight; box off; ylabel('Raw'); title(inFile)
    set(gca, 'xtick', []);
end

% merge consecutive blinks into 1 if they are X ms together
win          = hanning(wfilt);

if ~isempty(blinksmp)
cblinksmp = blinksmp(1,:);
for b = 1:size(blinksmp,1)-1,
    if blinksmp(b+1,1) - cblinksmp(end,2) < coalesce1 * dat.fsample,
        cblinksmp(end,2) = blinksmp(b+1,2);
    else
        cblinksmp(end+1,:) = blinksmp(b+1,:);
    end
end
blinksmp = cblinksmp; clear cblinksmp

% pad the blinks
padblinksmp(:,1) = round(blinksmp(:,1) + padding1(1) * dat.fsample);
padblinksmp(:,2) = round(blinksmp(:,2) + padding1(2) * dat.fsample);

% avoid idx outside range
if any(padblinksmp(:) < 1), padblinksmp(padblinksmp < 1) = 1; end
if any(padblinksmp(:) > length(dat.pupil)), padblinksmp(padblinksmp > length(dat.pupil)) = length(dat.pupil); end

% interpolate
[pupilcopy,~] = interp_nans(pupilcopy,padblinksmp);
[Xgazecopy2,~] = interp_nans(Xgazecopy,padblinksmp);
[Ygazecopy2,~] = interp_nans(Ygazecopy,padblinksmp);

% check to make sure all nans have been dealt with
assert(~any(isnan(pupilcopy)));

% low-pass filter so later-derived threshold is on same scale as signal it's being applied to
pupilcopy    = filter2(win.',pupilcopy,'same');


% ====================================================== %
% STEP 2: USE DERIVATIVE OF STEP 1 OUTPUT TO DEFINE INTERPOLATION WINDOWS
% ====================================================== %
% low-pass filter original time-series (otherwise particularly noisy measurements will yield very bad results)
pupilsmooth  = filter2(win.',dat.pupil,'same');

% calculate raw-unit derivative threshold from roughly interpolated time-series
raw_dthresh = diffthresh*std(diff(pupilcopy(wfilt:end-wfilt+1)));
raw_Xthresh = gazethresh*std(diff(Xgazecopy2(wfilt:end-wfilt+1)));
raw_Ythresh = gazethresh*std(diff(Ygazecopy2(wfilt:end-wfilt+1)));

% use threshold to find window start/end points
padblinksmp = [];
for b = 1:size(blinksmp,1),
    s1 = blinksmp(b,1);   % starting sample
    s2 = blinksmp(b,2);   % ending sample
    if s1-round(t_pre*dat.fsample)<wfilt  % if this starting sample is very close to start of timeseries
        s1 = 1;  % just take first sample as starting point
    else  % otherwise, searching for appropriate starting sample by applying threhsold to pupil derivative
        %while max(abs(diff(pupilsmooth((s1-round(t_pre*dat.fsample)):s1-1))))>raw_dthresh && s1-round(t_pre*dat.fsample)>wfilt && (s1-blinksmp(b,1))>-tmax*dat.fsample
        while (max(abs(diff(pupilsmooth((s1-round(t_pre*dat.fsample)):s1-1))))>raw_dthresh ||...
                 max(abs(diff(dat.Xgaze((s1-round(t_pre*dat.fsample)):s1-1))))>raw_Xthresh ||...
                 max(abs(diff(dat.Ygaze((s1-round(t_pre*dat.fsample)):s1-1))))>raw_Ythresh)...
                 && s1-round(t_pre*dat.fsample)>wfilt && (s1-blinksmp(b,1))>-tmax(1)*dat.fsample
            s1 = s1-1;
        end
        if s1-(t_pre*dat.fsample) == wfilt
            s1 = 1;
        end
    end
    if s2+(t_post*dat.fsample)>length(pupilsmooth)-wfilt+1  % if this ending sample is very close to end of timeseries
        s2 = length(pupilsmooth);  % just take last sample as ending point
    else  % otherwise, searching for appropriate starting sample by applying threhsold to pupil derivative
        %while max(abs(diff(pupilsmooth((s2+2:s2+round(t_post*dat.fsample))))))>raw_dthresh && s2+round(t_post*dat.fsample)<length(pupilsmooth)-wfilt+1 && (s2-blinksmp(b,2))<tmax*dat.fsample
        while (max(abs(diff(pupilsmooth((s2+2:s2+round(t_post*dat.fsample))))))>raw_dthresh ||...
                 max(abs(diff(dat.Xgaze((s2+2:s2+round(t_post*dat.fsample))))))>raw_Xthresh ||...
                 max(abs(diff(dat.Ygaze((s2+2:s2+round(t_post*dat.fsample))))))>raw_Ythresh)...
                 && s2+round(t_post*dat.fsample)<length(pupilsmooth)-wfilt+1 && (s2-blinksmp(b,2))<tmax(2)*dat.fsample
            s2 = s2+1;
        end
        if s2+round(t_post*dat.fsample) == length(pupilsmooth)-wfilt+1
            s2 = length(pupilsmooth);
        end
    end
    padblinksmp(b,1:2) = [s1 s2];
end

% interpolate
[dat.pupil,nanIdx] = interp_nans(dat.pupil,padblinksmp);
[dat.Xgaze,~] = interp_nans(dat.Xgaze,padblinksmp);
[dat.Ygaze,~] = interp_nans(dat.Ygaze,padblinksmp);

% check to make sure all nans have been dealt with
assert(~any(isnan(dat.pupil)));

% plot initial interpolation pass
if plotme, sp2 = subplot(511); hold on;
    plot(dat.times, dat.pupil, 'b');
    axis tight; box off; ylabel('Interp');
    set(gca, 'xtick', []);
end
end

% ====================================================== %
% STEP 3: USE DERIVATIVE OF STEP 2 OUTPUT TO IDENTIFY REMAINING BAD SAMPLES
% ====================================================== %
% low-pass filtering 
pupilsmooth  = filter2(win.',dat.pupil,'same');

% identify periods with derivative that exceeds harsh threshold
pupildiff = zscore(diff(pupilsmooth));   % calculate derivative of filtered pupil signal
Xgazediff = zscore(diff(dat.Xgaze));   % calculate derivative of X-gaze position
Ygazediff = zscore(diff(dat.Ygaze));   % calculate derivative of Y-gaze position
badsmp = find(abs(pupildiff) > diffthresh2*std(pupildiff));  % get positions of all outlying data points
badsmp = badsmp(badsmp>=wfilt & badsmp<=length(pupilsmooth)-wfilt+1);  % chopping off first and last n samples, which are contaminated by filtering

if ~isempty(badsmp)
    
    badpts = [badsmp(1) badsmp(find(abs(diff(badsmp))>1)+1)]';  % get samples where periods of outlying data points begin
    badpts(:,2) = [badsmp(abs(diff(badsmp))>1) badsmp(end)]';  % get samples where periods of outlying data points end
    
    % merge 2 windows into 1 if they are X ms together (since there will usually be a gap between down- and up-turns in derivative)
    cbadpts = badpts(1,:);
    for b = 1:size(badpts,1)-1,
        if badpts(b+1,1) - cbadpts(end,2) < coalesce2 * dat.fsample,
            cbadpts(end,2) = badpts(b+1,2);
        else
            cbadpts(end+1,:) = badpts(b+1,:);
        end
    end
    badpts = cbadpts; clear cbadpts
    
%     for b = 1:size(badpts, 1)-1,
%         if badpts(b+1, 1) - badpts(b, 2) < coalesce2 * dat.fsample,
%             badpts(b, 2) = badpts(b+1, 2);
%             badpts(b+1, :) = nan;
%         end
%     end
%     badpts(isnan(nanmean(badpts, 2)), :) = []; % remove those duplicates
    
    % calculate derivative thresholds for this round
    raw_dthresh = diffthresh*std(abs(diff(pupilsmooth(wfilt:end-wfilt+1))));
    raw_Xthresh = gazethresh*std(abs(diff(dat.Xgaze(wfilt:end-wfilt+1))));
    raw_Ythresh = gazethresh*std(abs(diff(dat.Ygaze(wfilt:end-wfilt+1))));
    
    % use threshold to find window start/end points
    padblinksmp = [];
    for b = 1:size(badpts,1),
        s1 = badpts(b,1);   % starting sample
        s2 = badpts(b,2);   % ending sample
        if s1-round(t_pre*dat.fsample)<wfilt  % if this starting sample is very close to start of timeseries
            s1 = 1;  % just take first sample as starting point
        else  % otherwise, searching for appropriate starting sample by applying threhsold to pupil derivative
            while (max(abs(diff(pupilsmooth((s1-round(t_pre*dat.fsample)):s1-1))))>raw_dthresh ||...
                    max(abs(diff(dat.Xgaze((s1-round(t_pre*dat.fsample)):s1-1))))>raw_Xthresh ||...
                    max(abs(diff(dat.Ygaze((s1-round(t_pre*dat.fsample)):s1-1))))>raw_Ythresh)...
                    && s1-round(t_pre*dat.fsample)>wfilt && (s1-badpts(b,1))>-tmax(1)*dat.fsample
                s1 = s1-1;
            end
            if s1-(t_pre*dat.fsample) == wfilt
                s1 = 1;
            end
        end
        if s2+(t_post*dat.fsample)>length(pupilsmooth)-wfilt+1  % if this ending sample is very close to end of timeseries
            s2 = length(pupilsmooth);  % just take last sample as ending point
        else  % otherwise, searching for appropriate starting sample by applying threhsold to pupil derivative
            while (max(abs(diff(pupilsmooth((s2+2:s2+round(t_post*dat.fsample))))))>raw_dthresh ||...
                    max(abs(diff(dat.Xgaze((s2+2:s2+round(t_post*dat.fsample))))))>raw_Xthresh ||...
                    max(abs(diff(dat.Ygaze((s2+2:s2+round(t_post*dat.fsample))))))>raw_Ythresh)...
                    && s2+round(t_post*dat.fsample)<length(pupilsmooth)-wfilt+1 && (s2-badpts(b,2))<tmax(2)*dat.fsample
                s2 = s2+1;
            end
            if s2+round(t_post*dat.fsample) == length(pupilsmooth)-wfilt+1
                s2 = length(pupilsmooth);
            end
        end
        padblinksmp(b,1:2) = [s1 s2];
    end
    
    % interpolate
    old_dat = dat.pupil;
    [dat.pupil,new_nans] = interp_nans(dat.pupil,padblinksmp);
    [dat.Xgaze,~] = interp_nans(dat.Xgaze,padblinksmp);
    [dat.Ygaze,~] = interp_nans(dat.Ygaze,padblinksmp);
    nanIdx(end+1:end+length(new_nans)) = new_nans;
    
    % check to make sure all nans have been dealt with
    assert(~any(isnan(dat.pupil)));
    
    % Plotting derivative and final interpolated timeseries
    if plotme, sp2 = subplot(512); hold on;
        plot(dat.times(11:end-10), Xgazediff(10:end-10), 'color', [0 1 0]);
        plot(dat.times(11:end-10), Ygazediff(10:end-10), 'color', [1 0 0]);
        plot(dat.times(11:end-10), pupildiff(10:end-10));  % trimming few edge samples because these can be very extreme
        plot(dat.times(new_nans), zeros(1,length(new_nans)), '.');
        box off; ylabel('Derivative');
        set(gca, 'xtick', []); ylim([-max(abs(pupildiff(11:end-10)))*1.02 max(abs(pupildiff(11:end-10)))*1.02]);
        
        sp3 = subplot(513); hold on;
        plot(dat.times, old_dat, 'color', [0.5 0.5 0.5]);
        plot(dat.times, dat.pupil, 'b');
        ylim([min(dat.pupil)*0.9 max(dat.pupil)*1.08]);
        box off; ylabel('Clean');
        
        sp4 = subplot(514); hold on;
        plot(dat.times, Xgazecopy, 'color', [0.5 0.5 0.5]);
        plot(dat.times, dat.Xgaze, 'color', [0 1 0]);
        ylim([min(dat.Xgaze)*0.8 max(dat.Xgaze)*1.2]);
        box off; ylabel('X-gaze');
        
        sp5 = subplot(515); hold on;
        plot(dat.times, Ygazecopy, 'color', [0.5 0.5 0.5]);
        plot(dat.times, dat.Ygaze, 'color', [1 0 0]);
        ylim([min(dat.Ygaze)*0.8 max(dat.Ygaze)*1.2]);
        box off; ylabel('Y-gaze');
        
        try
            linkaxes([sp1 sp2 sp3 sp4 sp5], 'x');
            set([sp1 sp2 sp3 sp4 sp5], 'tickdir', 'out');
        end
        xlim([dat.times(1) dat.times(end)]);
    end
    
    newblinksmp = padblinksmp;
else newblinksmp = [];
end

% specify output
newpupil = dat.pupil;
newXgaze = dat.Xgaze;
newYgaze = dat.Ygaze;




