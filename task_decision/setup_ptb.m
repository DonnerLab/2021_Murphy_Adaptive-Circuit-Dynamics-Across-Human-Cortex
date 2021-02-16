AssertOpenGL;
KbName('UnifyKeyNames');

screenNumber = min(Screen('Screens'));
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);


%% Open the screen
Screen('Preference', 'SkipSyncTests', 0);
[window, windowRect] = Screen('OpenWindow', screenNumber, [white/2 white/2 white/2]);  % [0 0 round(options.resolution./6)]
options.window_rect = windowRect;


%% Set the display parameters 'frameRate' and 'resolution'
options.frameDur     = Screen('GetFlipInterval', window); %duration of one frame
options.frameRate    = Screen('NominalFrameRate', window); %Hz

HideCursor(screenNumber)
Screen('Flip', window);


%% Set text appearance
Screen('Preference', 'TextRenderer', 1); % smooth text
Screen('TextFont', window, 'Helvetica'); % define text font
Screen('TextSize', window, 22); % define text font


%% Set up Eye Tracker
if strcmp(options.et,'yes')
    [el, options] = ELconfig(window, [subj,'_',sess,'_',b], options, screenNumber);
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
end


% make Kb Queue: Need to specify the device to query button box
% Find the keyboard + MEG buttons.
if strcmp(options.kb_setup,'MEG');
    [idx, names, all] = GetKeyboardIndices();
    options.kbqdev = [idx(strcmpi(names, 'ATEN USB KVMP w. OSD')), idx(strcmpi(names, 'Current Designs, Inc. 932')),...
        idx(strcmpi(names, 'Apple Internal Keyboard / Trackpad')), idx(strcmpi(names, ''))];
    
    keyList = zeros(1, 256);
    keyList(KbName({'ESCAPE','SPACE', 'LeftArrow', 'RightArrow',...
        '1', '2', '3', '4', 'b', 'g', 'y', 'r', '1!', '2@', '3#', '4$'})) = 1; % only listen to those keys!
    % first four are the buttons in mode 001, escape and space are for
    % the experimenter, rest is for testing
    for kbqdev = options.kbqdev
        PsychHID('KbQueueCreate', kbqdev, keyList);
        PsychHID('KbQueueStart', kbqdev);
        WaitSecs(.1);
        PsychHID('KbQueueFlush', kbqdev);
    end
end    


%% Maximum priority level
topPriorityLevel = MaxPriority(window);
