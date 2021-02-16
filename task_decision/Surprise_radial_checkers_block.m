% TASK: deccide which of two spatial distributions a sequence of dots was
% generated from. Distributions are Gaussian with shared SD but different
% means, and the generative distribution at a given point in time switches
% at a fixed hazard rate. Subjects task is to estimate which is the
% generative distribution at the *end* of a sequence. Sequence lengths are
% variable.


sca, clear
commandwindow;  % stops task response keys being written into task script

testing = 0;  % set to 1 for debugging version: brief stimuli, few trials


%% Path details
cd('C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Task')
datadir = 'C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Data\';


%% Basic PTB setup
options = setup;
if strcmp(options.do_trigger,'yes'), addpath matlabtrigger\, else addpath faketrigger\, end
trigger_enc = setup_trigger;


%% Subject details
if ~testing
    subj = input('Subject initials? ', 's');
    sess = input('Session number? ', 's');
else
    subj = 'test';
    sess = '1';
    options.do_trigger = 'no';
end

behavpath = [datadir,subj,filesep,'S',sess,filesep,'Behaviour',filesep];
stimpath = [datadir,subj,filesep,'S',sess,filesep,'Sample_seqs',filesep];
ELpath = [datadir,subj,filesep,'S',sess,filesep,'Eyetracking',filesep];

if ~isdir([datadir,subj,filesep,'S',sess])  % making new directories for this session if they don't exist, and setting block number to 1
      mkdir(behavpath)
      mkdir(stimpath)
      mkdir(ELpath)
      b = '1';  fprintf('\nNew subject/session combination...\n')
      WaitSecs(2);
else
    b = count_blocks([behavpath,'*.mat']);  % otherwise, counting number of existing files to determine current block number
    if isempty(b), b='1'; end   % catch just in case there are no saved behavioural files (if script errored before b1 end and now rerunning)
end
fprintf('\nCurrent block number is %s...\n',b)
WaitSecs(2);

% Suppress keyboard input to command line
ListenChar(2)


%% Fixed task settings
% Generative settings
gen.mu       =       [17 -17];     % means of generative distributions (polar angles relative to downward vertical midline of zero; + is left of midline)
gen.sigma    =        [29 29];     % standard deviations
gen.range    =   [-89.99 89.99];   % range for truncation
gen.H        =           0.08;     % hazard rate of distribution switches
maxsamps     =             12;     % maximum number of samples per trial

% Trial counts
ntrials   =    64;    % number of trials
pshort    =  0.25;    % proportion of trials with less than max number of samples (uniformly distributed between [2,maxsamps-1] samples)

% Timing (all in seconds)
timing.s         =     0.3;    % individual sample presentation time
timing.sgap      =     0.1;    % blank gap between samples
timing.s_refresh =     0.1;    % time between switch in checkerboard polarity within each sample (i.e. checker flicker rate)
timing.post      =       1;    % minimum time between final sample and response cue
timing.postdist  =     0.5;    % maximum increment to time between final sample and response cue (determines upper bound on uniform distribution)
timing.prefb     =     0.1;    % time between response and associated auditory feedback
timing.fb        =   0.125;    % duration of auditory feedback (for each of two consecutive tones)
timing.ITI       =     2.5;    % minimum time between feedback and next-trial onset
timing.rest      =     2.0;    % amount of fixed ITI time where blinks/rest is allowed
timing.ITIdist   =     1.5;    % maximum increment to time between feedback and next-trial onset (determines upper bound on uniform distribution)

% Stimulus appearance/positioning
stim.s_r             =    round(8.1*options.ppd);   % radial offset of regular samples (in d.v.a converted to pixels)
stim.cbar_r          =                       8.8;   % radial offset of LLR colour bar (in d.v.a)
stim.cbar_rext       =                       0.4;   % radial extent of LLR colour bar
stim.yscaling        =                         0;   % multiplicative scaling factor by which to decrease y-offset of stimuli (to mimic x/y asymmetry in human vision)
stim.fix_in_r        =   round(0.18*options.ppd);   % radius of inner fixation point
stim.fix_out_r       =   round(0.36*options.ppd);   % radius of outer fixation point
stim.fix_in_c        =             [  0   0   0];   % default colour of inner fixation point
stim.fix_out_c       =             [168 142 160];   % colour of outer fixation point for active trial period (REDISH)
stim.fix_out_c_resp  =             [134 171 149];   % colour of outer fixation point for response cueing (GREENISH)
stim.fix_out_rest    =             [128 153 169];   % colour of outer fixation point for rest period (BLUEISH)
stim.fb_freqs        =             [350 950 nan];   % frequencies of auditory feedback tones (nan = white noise)

% Checkerboard patch settings
cb.size     =     0.8;   % size of one side of square within which checkerboard will be drawn (d.v.a.)
cb.freq     =       2;   % spatial frequency (cycles per d.v.a.)

% Make LLR colour bar
load('C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Task\Teufel_rgb.mat')  % loading pre-computed Tuefel colourbar
ncols = size(rgb,1);   % number of discrete colours in LLR colourbar
rgb = rgb(size(rgb,1):-1:1,:);  % re-arranging rgbs such that green = left, red = right

[cbar_img,cbar_alpha] = make_gaussian_LLR_cbar_tex(gen.mu,gen.sigma,gen.range,options.ppd,options.resolution,stim.cbar_r,stim.cbar_rext,stim.yscaling);  % creating LLR colourbar image
cbar_img = round(cbar_img.*(ncols-1))+1;  % rescaling colorbar images for appropriate colouration
cols = rgb';
tmp = nan(size(cbar_img)); tmp(~isnan(cbar_img)) = cols(1,cbar_img(~isnan(cbar_img))); cbar_rgbs = tmp;  % translating LLRs into desired colourmap
tmp = nan(size(cbar_img)); tmp(~isnan(cbar_img)) = cols(2,cbar_img(~isnan(cbar_img))); cbar_rgbs(:,:,2) = tmp;
tmp = nan(size(cbar_img)); tmp(~isnan(cbar_img)) = cols(3,cbar_img(~isnan(cbar_img))); cbar_rgbs(:,:,3) = tmp;
cbar_rgbs(isnan(cbar_rgbs)) = 127;
cbar_rgbs(:,:,4) = cbar_alpha;

% Construct start/end points for vertical midline
stim.mid_xy      = [0; stim.fix_out_r + 0.35*options.ppd];                % start1
stim.mid_xy(:,2) = [0; stim.s_r - ceil(cb.size*(1+stim.yscaling+0.35)*options.ppd/2)];           % end1
stim.mid_xy(:,3) = [0; round((stim.cbar_r)*options.ppd)];           % start2
stim.mid_xy(:,4) = [0; round((stim.cbar_r+stim.cbar_rext)*options.ppd)];  % end2
stim.mid_xy = stim.mid_xy.*(1-stim.yscaling);   % rescaling y-offset
stim.mid_w = 3;   % linewidth

% If only testing, overwrite some of the above
if testing == 1
    ntrials = 6;
    timing.s = 0.3;
    timing.post = 0.25;
    timing.ITIdist = 0.1;
end


%% Initialize Psychtoolbox and create textures
setup_ptb;
if strcmp(options.et,'yes'), PsychPortAudio('Close'); end
stim.audio = setup_audio(stim.fb_freqs,timing.fb);
timing.ifi = options.frameDur;     % inter-frame interval (1/sampling rate)

[gabor.tex,gabor.rect] = createCircularChecker(window, cb.size, cb.freq, options.ppd, 0);
gabor.rect = gabor.rect-(gabor.rect(4)/2);    % making sure subsequent coordinates are always wrt screen center

stim.cbartex = Screen('MakeTexture', window, cbar_rgbs);  % create colourbar texture


try
    %% Present intro screen
    block_instructions(b,window,windowRect);
    
    
    %% Generate sample sequences for all trials
    % Seed random number generator
    rng('default')
    rng('shuffle')
    
    % Generate sequences of distribution switch positions
    pswitch = [ones(ntrials,1) rand(ntrials,maxsamps-1)];  % first draw uniformly distributed random probabilities b/w 0 and 1 that will determine switch positions (first sample is never a switch)
    pswitch(pswitch>gen.H) = 0; pswitch(pswitch~=0) = 1;  % binarize matrix to mark only switch positions
    
    % Generate sequences of which distributions will be drawn from at which times (1=left, 2=right)
    distseqs = zeros(ntrials,maxsamps); dists = [1 2];
    for t = 1:ntrials
        if t<=ntrials/2, cdist = 1; else cdist = 2; end  % making sure each distribution is starting distribution an equal number of times
        s = 0;
        while s<maxsamps
            s = s+1;
            if ~pswitch(t,s), distseqs(t,s) = cdist;  % if switch in distribution has not occured
            else cdist = dists(dists~=cdist);  % if switch in distribution has occured
                distseqs(t,s) = cdist;
            end
        end
    end
    
    % Generate actual evidence samples from distribution sequences
    stimIn = distseqs;
    stimIn(distseqs==1) = round(gen.mu(1)+(randn(size(stimIn(distseqs==1))).*gen.sigma(1)));
    stimIn(distseqs==2) = round(gen.mu(2)+(randn(size(stimIn(distseqs==2))).*gen.sigma(2)));
    stimIn(stimIn<gen.range(1)) = gen.range(1); stimIn(stimIn>gen.range(2)) = gen.range(2);  % in case drawn values exceed range limits
    
    % Shuffle trial order
    shuforder = randperm(ntrials);
    pswitch = pswitch(shuforder,:);
    stimIn = stimIn(shuforder,:);
    distseqs = distseqs(shuforder,:);
    
    % Store generative distributions @ sequence end (i.e. correct choices; 1=left, 2=right)
    fdists = distseqs(:,end);
    
    % Truncate sequence length of random selection of trials
    tlengths = randsample(2:(maxsamps-1),round(ntrials*pshort),true);  % randomly draw lengths of truncated sequences
    ttrials = sort(randsample(ntrials,round(ntrials*pshort),false))';  % randomly draw trials to be truncated
    for t = 1:length(ttrials)
        pswitch(ttrials(t),tlengths(t)+1:end) = nan;
        stimIn(ttrials(t),tlengths(t)+1:end) = nan;
        distseqs(ttrials(t),tlengths(t)+1:end) = nan;
        fdists(ttrials(t)) = distseqs(ttrials(t),tlengths(t));
    end
    
    % Save trial-by-trial sample sequences, switch positions, etc.
    save([stimpath,subj,'_',sess,'_',b,'.mat'],'pswitch','stimIn','distseqs','gen','pshort','timing','stim')
    
    
    %% Transform dot positions from polar to cartesian (x,y) coordinates
    [X,Y] = RectCenter(windowRect);
    [XstimIn,YstimIn] = pol2cart(deg2rad(stimIn+90),stim.s_r);  % get XY coords relative to origin of zero
    YstimIn = YstimIn.*(1-stim.yscaling);  % decrease y-offset of dots by multiplicative scaling factor
    XstimIn = XstimIn+ X; YstimIn = YstimIn+Y;  % reference to screen center
    
    
    %% Start Eyelink recording
    if strcmp(options.et,'yes');
        Eyelink('StartRecording');
        WaitSecs(0.1);
        Eyelink('message', 'Start recording Eyelink');
    end
    
    
    %% Countdown and first fixation point
    str1 = 'Ready to go!';
    cd_sec = 3;
    
    while cd_sec>0
        str2 = num2str(cd_sec);
        DrawFormattedText(window,str1,'center',Y-25,[255 255 255]);
        DrawFormattedText(window,str2,'center',Y+25,[255 255 255]);
        Screen('Flip', window);
        cd_sec=cd_sec-1;
        WaitSecs(0.96);
    end
    
    pos_in = [X-stim.fix_in_r; Y-stim.fix_in_r; X+stim.fix_in_r; Y+stim.fix_in_r];
    pos_out = [X-stim.fix_out_r; Y-stim.fix_out_r; X+stim.fix_out_r; Y+stim.fix_out_r];
    
    Screen('FillOval', window, stim.fix_out_c, pos_out, stim.fix_out_r*2.1);
    Screen('FillOval', window, stim.fix_in_c, pos_in, stim.fix_out_r*2.1);
    vbl = Screen('Flip', window);   % then fixation plus patches
    
    
    %% Loop through trials
    trigger(trigger_enc.block_start);  % trigger to mark start of block
    if strcmp(options.et,'yes'); Eyelink('message', num2str(trigger_enc.block_start)); end
    
    on_time = [];  % onset time for first trial is specified within trial script
    Behav = zeros(length(fdists),6); tRefresh = zeros(length(fdists),maxsamps);
    for t = 1:length(fdists);
        
        % Presenting trial number at the bottom of the eyetracker display
        if strcmp(options.et,'yes');
            Eyelink('command', 'record_status_message "TRIAL %d/%d"', t, length(fdists));
            Eyelink('message', 'TRIALID %d', t);
        end
        
        % Variable task parameters
        varopts.on_time = on_time;                        % controls onset time of impending trial - fed back from trial function
        varopts.xpos = XstimIn(t,~isnan(XstimIn(t,:)));   % sequence of sample X-coordinates
        varopts.ypos = YstimIn(t,~isnan(YstimIn(t,:)));   % sequence of sample Y-coordinates
        varopts.fdist = fdists(t);                        % correct answer (i.e. generative distribution @ sequence end: 1=left, 2=right)
        varopts.kbqdev = options.kb_setup;                % keyboard info
        varopts.ppd = options.ppd;                        % pixels per degree
        
        % Run trial
        [Behav(t,1:6),tRefresh(t,1:length(varopts.xpos)),on_time] = Surprise_radial_checkers_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, t, strcmp(options.et,'yes'));
    end
    
    WaitSecs(3);
    trigger(trigger_enc.block_end);  % trigger to mark end of block
    if strcmp(options.et,'yes'); Eyelink('message', num2str(trigger_enc.block_end)); end
    
    
    %% Save behavioral data
    fprintf('\nSaving behavioural data to %s\n', behavpath)
    save([behavpath,subj,'_',sess,'_',num2str(b),'.mat'],'Behav','tRefresh')
    
    
    %% Calculate accuracy from current + previous blocks and display to subject
    end_block_screen(str2double(b),[behavpath,subj,'_',sess,'_'],window,windowRect);
    
    
    %% Save Eyelink data
    if strcmp(options.et,'yes');
        fprintf('Saving EyeLink data to %s\n', ELpath)
        eyefilename = fullfile(ELpath,options.edfFile);
        Eyelink('CloseFile');
        Eyelink('WaitForModeReady', 500);
        try
            status = Eyelink('ReceiveFile', options.edfFile, eyefilename);
            disp(['File ' eyefilename ' saved to disk']);
        catch
            warning(['File ' eyefilename ' not saved to disk']);
        end
        Eyelink('StopRecording');
    end
    
    
    %% Exit
    ListenChar(0)
    PsychPortAudio('DeleteBuffer');
    PsychPortAudio('Close',stim.audio.h);
    sca
    
    
catch ME  % if any errors are encountered or user requests quit, clean up and exit
    ListenChar(0)
    PsychPortAudio('DeleteBuffer');
    PsychPortAudio('Close',stim.audio.h);
    sca
    rethrow(ME)
end


