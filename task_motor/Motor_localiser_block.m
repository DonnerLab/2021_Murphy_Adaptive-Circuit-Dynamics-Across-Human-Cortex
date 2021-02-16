% TASK: deccide which of two spatial distributions a sequence of dots was
% generated from. Distributions are Gaussian with shared SD but different
% means, and the generative distribution at a given point in time switches
% at a fixed hazard rate. Subjects task is to estimate which is the
% generative distribution at the *end* of a sequence. Sequence lengths are
% variable.


sca, clear

testing = 1;  % set to 1 for debugging version: brief stimuli, few trials


%% Path details
cd D:\Experiments\Surprise_accumulation\Task_motor_localiser
datadir = 'D:\Experiments\Surprise_accumulation\Task_motor_localiser\Data\';


%% Basic PTB setup
options = setup;
if testing, options.et = 'no'; options.do_trigger = 'no'; end
if strcmp(options.do_trigger,'yes'), addpath matlabtrigger/, else addpath faketrigger/, end
trigger_enc = setup_trigger;


%% Subject details & basic PTB setup
if ~testing
    subj = input('Subject initials? ', 's');
    sess = input('Session number? ', 's');
else
    subj = 'test';
    sess = '1';
    options.et = 'no'; options.do_trigger = 'no';
end

behavpath = [datadir,subj,filesep,'S',sess,filesep,'Behaviour',filesep];
ELpath = [datadir,subj,filesep,'S',sess,filesep,'Eyetracking',filesep];

if ~isdir([datadir,subj,filesep,'S',sess])  % making new directories for this session if they don't exist, and setting block number to 1
      mkdir(behavpath)
      mkdir(ELpath)
      b = '1';  fprintf('\nNew subject/session combination...\n')
      WaitSecs(2);
else
    b = count_blocks([behavpath,'*.mat']);  % otherwise, counting number of existing files to determine current block number
    if isempty(b), b='1'; end   % catch just in case there are no saved behavioural files (if script errored before b1 end and now rerunning)
end
fprintf('\nCurrent block number is %s...\n',b)
WaitSecs(2);


%% Fixed task settings
% Trial count
ntrials =    60;    % number of trials

% Timing (all in seconds)
timing.precue    =   0.75;    % minimum time for pre-cue active fixation
timing.precueD   =   0.75;    % maximum increment to time for pre-cue active fixation (determines upper bound on uniform distribution)
timing.cue       =    0.3;    % duration of cue presentation
timing.prep      =    1.0;    % post-cue motor preparation time
timing.rest      =    2.0;    % duration of post-response rest period

% Stimulus appearance/positioning
stim.cueYoff          =         -1.25*options.ppd;   % vertical offset (relative to screen center) of written response cue
stim.cuesize          =                        15;   % font size of response cue
stim.cuefont          =               'Trebuchet';   % font type of response cue
stim.cuecol           =             [255 255 255];   % font color of response cue

stim.fix_in_r         =   round(0.18*options.ppd);   % radius of inner fixation point
stim.fix_out_r        =   round(0.36*options.ppd);   % radius of outer fixation point
stim.fix_in_c         =             [  0   0   0];   % colour of inner fixation point
stim.fix_out_c_prep   =             [168 142 160];   % colour of outer fixation point during preparation
stim.fix_out_c_go     =             [134 171 149];   % colour of outer fixation point during response execution
stim.fix_out_c_rest   =             [128 153 169];   % colour of outer fixation point during rest period

% Make trial list [1 = left, 2 = right]
rng('default')
rng('shuffle')

resps = [ones(ceil(ntrials/2),1); ones(ceil(ntrials/2),1).*2];
resps = Shuffle(resps);

% If only testing, overwrite some of the above
if testing == 1
    resps = resps(1:4);
end


%% Initialize Psychtoolbox and create textures
setup_ptb;
timing.ifi = options.frameDur;     % inter-frame interval (1/sampling rate)


%% Present intro screen
block_instructions(b,window,windowRect);


%% Start Eyelink recording
if strcmp(options.et,'yes');
    Eyelink('StartRecording');
    WaitSecs(0.1);
    Eyelink('message', 'Start recording Eyelink');
end


%% Set response type cue appearance
Screen('TextSize', window, stim.cuesize);
Screen('TextFont', window, stim.cuefont);
Screen('TextColor', window, stim.cuecol);

%% Countdown and first fixation point
[X,Y] = RectCenter(windowRect);
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

Screen('FillOval', window, stim.fix_out_c_rest, pos_out, stim.fix_out_r*2.1);
Screen('FillOval', window, stim.fix_in_c, pos_in, stim.fix_out_r*2.1);
Screen('Flip', window);   % present fixation


%% Loop through trials
trigger(trigger_enc.block_start);  % trigger to mark start of block
if strcmp(options.et,'yes'); Eyelink('message', num2str(trigger_enc.block_start)); end

on_time = [];  % onset time for first trial is specified within trial script
Behav = zeros(length(resps),6);
for t = 1:length(resps);
    
    % Presenting trial number at the bottom of the eyetracker display
    if strcmp(options.et,'yes');
        Eyelink('command', 'record_status_message "TRIAL %d/%d"', t, length(fdists));
        Eyelink('message', 'TRIALID %d', t);
    end
    
    % Variable task parameters
    varopts.on_time = on_time;            % controls onset time of impending trial - fed back from trial function
    varopts.cresp = resps(t);             % response cue type
    varopts.kbqdev = options.kbqdev;      % keyboard info
    
    % Run trial
    [Behav(t,1:6),on_time] = Motor_localiser_trial(window, windowRect, timing, stim, varopts, trigger_enc, t, strcmp(options.et,'yes'));
end

WaitSecs(3);
trigger(trigger_enc.block_end);  % trigger to mark end of block
if strcmp(options.et,'yes'); Eyelink('message', num2str(trigger_enc.block_end)); end


%% Save behavioral data
fprintf('\nSaving behavioural data to %s\n', behavpath)
save([behavpath,subj,'_',sess,'_',num2str(b),'.mat'],'Behav')


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
sca


