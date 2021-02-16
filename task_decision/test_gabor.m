sca, clear


%% Path details
cd D:\Experiments\Surprise_accumulation\Task


%% Basic PTB setup
options = setup;
setup_ptb;

% Dot appearance/positioning
stim.s_y             =      0*options.ppd;   % vertical offset of regular samples (d.v.a. relative to center of screen, converted to pixels)
stim.s_r             =   round(0.15*options.ppd);   % radius of sample gabor patches

[X,Y] = RectCenter(windowRect);

Screen('DrawLine', window, [255 255 255 0], X, Y-100, X, Y+stim.s_y*2)


cs = [255;255;255]; 
ps = [X-stim.s_r; Y+(stim.s_y-stim.s_r); X+stim.s_r; Y+(stim.s_y+stim.s_r)];  % adjust position of dot
Screen('FillOval', window, cs, ps);   % drawing fixation and dot in one call

% Gabor patch settings
g_size           =    round(6.69 *options.ppd) ;   % size of one side of square within which gabor will be drawn (d.v.a. converted to pixels)
g_sigma          =    6.69  *options.ppd;   % standard deviation of Gaussian
g_tilt           =                   180   ;   % orientation (0 = vertically oriented; positive increments yield counter-clockwise rotation)
g_contrast       =                0.2 ;   % contrast
g_aspectRatio    =                1.0;   % determines elongation (whether patch is square<--->ellipse)
g_phase          =                  0;   % phase
g_freq           =      0.65  /options.ppd;   % spatial frequency (cycles per d.v.a., converted to cycles per pixel)


backgroundOffset = [0.5 0.5 0.5 0.0];
disableNorm = 1;
preContrastMultiplier = 0.5;
[gabortex,gaborrect] = CreateProceduralGabor(window, g_size, g_size, [],...
    backgroundOffset, disableNorm, preContrastMultiplier);
gaborrect = gaborrect-(gaborrect(4)/2)-1;

% Randomise the phase of the Gabors and make a properties matrix.
propertiesMat = [g_phase, g_freq, g_sigma, g_contrast, g_aspectRatio, 0, 0, 0];

% Draw the Gabor. By default PTB will draw this in the center of the screen
% for us.
g_ps = OffsetRect(gaborrect, X, Y);
Screen('DrawTextures', window, gabortex, [], g_ps, g_tilt, [], [], [], [],...
    kPsychDontDoRotation, propertiesMat'); 

% Flip to the screen
Screen('Flip', window);

% Wait for a button press to exit
KbWait;

sca