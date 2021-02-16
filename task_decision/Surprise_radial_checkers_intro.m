% TASK: deccide which of two spatial distributions a sequence of dots was
% generated from. Distributions are Gaussian with shared SD but different
% means, and the generative distribution at a given point in time switches
% at a fixed hazard rate. Subjects task is to estimate which is the
% generative distribution at the *end* of a sequence. Sequence lengths are
% variable.


sca, clear
commandwindow;  % stops task response keys being written into task script
ListenChar(2)  % suppress keyboard input to command line


%% Path details
cd('C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Task')
datadir = 'C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Data\';
 

%% Basic PTB setup
options = setup;
if strcmp(options.do_trigger,'yes'), addpath matlabtrigger\, else addpath faketrigger\, end
trigger_enc = setup_trigger;
options.et = 'no';

% Seed random number generator
rng('default')
rng('shuffle')


%% Fixed task settings
% Generative settings
gen.mu       =         [17 -17];   % means of generative distributions (polar angles relative to downward vertical midline of zero; + is left of midline)
gen.sigma    =          [29 29];   % standard deviations
gen.range    =   [-89.99 89.99];   % range for truncation
gen.H        =           0.08;     % hazard rate of distribution switches
maxsamps     =             12;     % maximum number of samples per trial

% Trial counts
example_samps  =    30;    % number of samples to present from example distribution

ntrials_noCP   =    12;    % number of initial practice trials without any change-points (default = 12)
pshort_noCP    =   0.5;    % proportion of trials with less than max number of samples (uniformly distributed between [2,maxsamps-1] samples)

ntrials_CP     =    20;    % number of initial practice trials without any change-points (default = 20)
pshort_CP      =  0.25;

% Timing (all in seconds)
timing.s         =     0.3;    % individual sample presentation time
timing.s_ex      =     0.4;    % individual sample presentation time for initial example distributions
timing.s_refresh =     0.1;    % time between switch in checkerboard polarity within each sample (i.e. checker flicker rate)
timing.sgap      =     0.1;    % blank gap between samples
timing.post      =       1;    % minimum time between final sample and response cue
timing.postdist  =     0.5;    % maximum increment to time between final sample and response cue (determines upper bound on uniform distribution)
timing.prefb     =     0.1;    % time between response and associated auditory feedback
timing.fb        =   0.125;    % duration of auditory feedback (for each of two consecutive tones)
timing.fbtext    =     0.5;    % duration of visual feedback (training only)
timing.ITI       =     2.5;    % minimum time between feedback and next-trial onset
timing.rest      =     2.0;    % amount of fixed ITI time where blinks/rest is allowed
timing.ITIdist   =     0.5;    % maximum increment to time between feedback and next-trial onset (determines upper bound on uniform distribution)

% Stimulus appearance/positioning
stim.s_r             =    round(8.1*options.ppd);   % radial offset of regular samples (in d.v.a converted to pixels)
stim.cbar_r          =                       8.8;   % radial offset of LLR colour bar (in d.v.a)
stim.cbar_rext       =                      0.40;   % radial extent of LLR colour bar
stim.yscaling        =                         0;   % multiplicative scaling factor by which to decrease y-offset of stimuli (to mimic x/y asymmetry in human vision)
stim.fix_in_r        =   round(0.18*options.ppd);   % radius of inner fixation point
stim.fix_out_r       =   round(0.36*options.ppd);   % radius of outer fixation point
stim.fix_in_c        =             [  0   0   0];   % default colour of inner fixation point
stim.fix_out_c       =             [168 142 160];   % colour of outer fixation point for active trial period (REDISH)
stim.fix_out_c_resp  =             [134 171 149];   % colour of outer fixation point for response cueing (GREENISH)
stim.fix_out_rest    =             [128 153 169];   % colour of outer fixation point for rest period (BLUEISH)
stim.txt_yoffset     =   round(1.45*options.ppd);   % y-offset of feedback/change-point text (training only)
stim.fb_freqs        =             [350 950 nan];   % frequencies of auditory feedback tones (nan = white noise)

% Checkerboard patch settings
cb.size     =     0.8;   % size of one side of square within which checkerboard will be drawn (d.v.a.)
cb.freq     =       2;   % spatial frequency (cycles per d.v.a.)

% Make distribution PDF colour bars
load('C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Task\Teufel_rgb.mat')  % loading pre-computed Tuefel colourbar
ncols = size(rgb,1);   % number of discrete colours in LLR colourbar
rgb = rgb(size(rgb,1):-1:1,:);  % re-arranging rgbs such that green = left, red = right

[cbar_img1,cbar_alpha1] = make_gaussian_PDF_cbar_tex(gen.mu(1),gen.sigma(1),gen.range,options.ppd,options.resolution,stim.cbar_r,stim.cbar_rext,stim.yscaling);  % creating PDF colourbar image
[cbar_img2,cbar_alpha2] = make_gaussian_PDF_cbar_tex(gen.mu(2),gen.sigma(2),gen.range,options.ppd,options.resolution,stim.cbar_r,stim.cbar_rext,stim.yscaling);  % creating PDF colourbar image
[cbar_img3,cbar_alpha3] = make_gaussian_PDF_cbar_tex(0,15,gen.range,options.ppd,options.resolution,stim.cbar_r,stim.cbar_rext,stim.yscaling);  % creating PDF colourbar image
[cbar_img,cbar_alpha] = make_gaussian_LLR_cbar_tex(gen.mu,gen.sigma,gen.range,options.ppd,options.resolution,stim.cbar_r,stim.cbar_rext,stim.yscaling);  % creating LLR colourbar image

cbar_img1 = round(cbar_img1.*(ncols-1))+1;  % rescaling colorbar images for appropriate colouration
cbar_img2 = round(cbar_img2.*(ncols-1))+1;
cbar_img3 = round(cbar_img3.*(ncols-1))+1;
cbar_img = round(cbar_img.*(ncols-1))+1;

base_col = 128;  % baseline grey level for low PDFs
max_col = 180;  % maximum color level for high PDFs
cols1 = [linspace(base_col,255-max_col,ncols); linspace(base_col,max_col,ncols); linspace(base_col,255-max_col,ncols)];  % creating color maps
cols2 = [linspace(base_col,max_col,ncols); linspace(base_col,255-max_col,ncols); linspace(base_col,255-max_col,ncols)];
cols3 = [linspace(base_col,255,ncols); linspace(base_col,255,ncols); linspace(base_col,255,ncols)];
cols = rgb';

tmp = nan(size(cbar_img1)); tmp(~isnan(cbar_img1)) = cols1(1,cbar_img1(~isnan(cbar_img1))); cbar_rgbs1 = tmp;  % translating PDFs into desired colourmap
tmp = nan(size(cbar_img1)); tmp(~isnan(cbar_img1)) = cols1(2,cbar_img1(~isnan(cbar_img1))); cbar_rgbs1(:,:,2) = tmp;
tmp = nan(size(cbar_img1)); tmp(~isnan(cbar_img1)) = cols1(3,cbar_img1(~isnan(cbar_img1))); cbar_rgbs1(:,:,3) = tmp;
cbar_rgbs1(isnan(cbar_rgbs1)) = 127;
cbar_rgbs1(:,:,4) = cbar_alpha1;
tmp = nan(size(cbar_img2)); tmp(~isnan(cbar_img2)) = cols2(1,cbar_img2(~isnan(cbar_img2))); cbar_rgbs2 = tmp;  % translating PDFs into desired colourmap
tmp = nan(size(cbar_img2)); tmp(~isnan(cbar_img2)) = cols2(2,cbar_img2(~isnan(cbar_img2))); cbar_rgbs2(:,:,2) = tmp;
tmp = nan(size(cbar_img2)); tmp(~isnan(cbar_img2)) = cols2(3,cbar_img2(~isnan(cbar_img2))); cbar_rgbs2(:,:,3) = tmp;
cbar_rgbs2(isnan(cbar_rgbs2)) = 127;
cbar_rgbs2(:,:,4) = cbar_alpha2;
tmp = nan(size(cbar_img3)); tmp(~isnan(cbar_img3)) = cols3(1,cbar_img3(~isnan(cbar_img2))); cbar_rgbs3 = tmp;  % translating PDFs into desired colourmap
tmp = nan(size(cbar_img3)); tmp(~isnan(cbar_img3)) = cols3(2,cbar_img3(~isnan(cbar_img2))); cbar_rgbs3(:,:,2) = tmp;
tmp = nan(size(cbar_img3)); tmp(~isnan(cbar_img3)) = cols3(3,cbar_img3(~isnan(cbar_img2))); cbar_rgbs3(:,:,3) = tmp;
cbar_rgbs3(isnan(cbar_rgbs3)) = 127;
cbar_rgbs3(:,:,4) = cbar_alpha3;
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
stim.mid_w = 4;   % linewidth


%% Initialize Psychtoolbox and create textures
setup_ptb;

blah

stim.audio = setup_audio(stim.fb_freqs,timing.fb);
timing.ifi = options.frameDur;     % inter-frame interval (1/sampling rate)

[gabor.tex,gabor.rect] = createCircularChecker(window, cb.size, cb.freq, options.ppd, 0);
[gabor_ex.tex,gabor.rect] = createCircularChecker(window, cb.size, cb.freq, options.ppd, 0.4);

gabor.rect = gabor.rect-(gabor.rect(4)/2);    % making sure subsequent coordinates are always wrt screen center

stim.cbartex1 = Screen('MakeTexture', window, cbar_rgbs1);  % create texture for first distribution
stim.cbartex2 = Screen('MakeTexture', window, cbar_rgbs2);  % create texture for second distribution
stim.cbartex3 = Screen('MakeTexture', window, cbar_rgbs3);  % create texture for first example distribution
stim.cbartex = Screen('MakeTexture', window, cbar_rgbs);  % create texture for second distribution


%% Get some basic positioning stuff
[X,Y] = RectCenter(windowRect);  % center coords

stim.mid_xy(1,:) = stim.mid_xy(1,:)+X;  % centering reference line
stim.mid_xy(2,:) = stim.mid_xy(2,:)+Y;

pos_in = [X-stim.fix_in_r; Y-stim.fix_in_r; X+stim.fix_in_r; Y+stim.fix_in_r];  % fixation coords
pos_out = [X-stim.fix_out_r; Y-stim.fix_out_r; X+stim.fix_out_r; Y+stim.fix_out_r];

ps = [pos_out pos_in];  % positions
cs = [stim.fix_out_c' stim.fix_in_c'];  % colours


try
    %% Present initial instruction screen
    Screen('TextSize', window, 18);
    
    str1 = ['Welcome!\n \n \n',...
        'In this study, you will be shown sequences of small dots presented around\na central fixation point. You should always look at this\n'...
        'fixation point when the dots are being presented.\n \n \n',...
        'Importantly, the positions of the dots can be\ndrawn from one of two ''DECKS''. Your task will be to decide which\n',...
        'deck is being used to draw dot positions at different times.'];
    str2 = 'Press any key to see an example deck...';
    
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.26,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.9,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    KbWait;
    WaitSecs(1);
    
    
    %% Present initial example distribution with relatively low variance
    str1 = ['The curved white bar here shows the spread of different\ndot positions that could be drawn from an example deck.\n \n',...
        'Whiter places along the bar indicate the *most common* dot positions\ncontained in this deck, whereas less white places indicate\ndot positions that are quite rare in this deck.'];
    str4 = 'Press any key for more information...';
    
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    Screen('DrawTexture', window, stim.cbartex3);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.18,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str4,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('Flip', window);
    
    KbWait;
    WaitSecs(1);
    
    str2 = 'This is the central fixation point that\nyou need to always look at when\nthe dots are being shown';
    str3 = 'This is a useful reference\nline that indicates the\nmiddle of the screen';
    str4 = 'Press any key to draw some dot positions from this example deck...';
    
    Screen('DrawTexture', window, stim.cbartex3);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.18,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str4,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('TextSize', window, 15);
    DrawFormattedText(window,str2,X-370,Y-35,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str3,X+100,Y+1,[255 255 255],[],[],[],2);
    Screen('DrawLines', window, [X+15 X+92; Y+100 Y+48], 1, [255 255 255]);   % drawing pointer
    Screen('DrawLines', window, [X-160 X-3; Y Y], 1, [255 255 255]);   % drawing pointer
    Screen('Flip', window);
    
    KbWait;
    
    Screen('DrawTexture', window, stim.cbartex3);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.18,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    WaitSecs(1);
    
    pos1 = randn(1,example_samps).*15;  % sample positions
    [XstimIn,YstimIn] = pol2cart(deg2rad(pos1+90),stim.s_r);  % get XY coords relative to origin of zero
    YstimIn = YstimIn.*(1-stim.yscaling);  % decrease y-offset of dots by multiplicative scaling factor
    XstimIn = XstimIn+ X; YstimIn = YstimIn+Y;  % reference to screen center
    g_ps = OffsetRect(gabor.rect, XstimIn', YstimIn')';   % get positions of checkerboards
    vbl = GetSecs;  % get start time
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    for s = 1:length(pos1)
        % Sample on
        Screen('FillOval', window, cs, ps);   % drawing fixation point
        Screen('DrawTextures', window, gabor_ex .tex(1), [], g_ps(:,1:s));  % drawing checkerboards
        Screen('DrawTexture', window, stim.cbartex3);  % drawing colourbar
        Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
        DrawFormattedText(window,str1,'center',options.resolution(2)*0.18,[255 255 255],[],[],[],2);
        vbl = Screen('Flip', window, vbl + timing.s_ex  + timing.sgap - timing.ifi*0.5);  % set next dot to be flipped at specified dot presentation rate
    end
    
    WaitSecs(1);
    
    str3 = 'Notice how most of the dot positions came from\nthis bright region of the example deck';
    str4 = 'Press any key to continue...';
    
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('DrawTextures', window, gabor_ex .tex(1), [], g_ps);  % drawing checkerboards
    Screen('DrawTexture', window, stim.cbartex3);  % drawing colourbar
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.18,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str4,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('TextSize', window, 15);
    DrawFormattedText(window,str3,'center',options.resolution(2)*0.82,[255 255 255],[],[],[],2);
    Screen('DrawLines', window, [X-2 X; Y+20+options.ppd*stim.cbar_r*(1-stim.yscaling) options.resolution(2)*0.82], 1, [255 255 255]);   % drawing pointer
    Screen('Flip', window);  % set next dot to be flipped at specified dot presentation rate
    
    KbWait
    
    
    %% Present distribution 1 and draw a sequence of samples
    str1 = ['Now let''s have a look at the two decks that we will\nactually be using in the study.\n \n \n',...
        'As you will see, the most common dot positions from each of these decks are different;\nhowever, there is also quite a lot of overlap between the decks.'];
    str2 = 'Press any key to see DECK 1...';
    
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.35,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.9,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    str1 = 'Here''s what the distribution of positions from DECK 1 looks like.';
    str2 = 'Brighter green colors indicate the most common positions contained in this deck:';
    str4 = 'Press any key to draw some positions from DECK 1...';
    
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    Screen('DrawTexture', window, stim.cbartex1);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.21-20,[255 255 255]);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.21+20,[255 255 255]);
    DrawFormattedText(window,str4,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    Screen('DrawTexture', window, stim.cbartex1);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.21-20,[255 255 255]);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.21+20,[255 255 255]);
    Screen('Flip', window);
    
    WaitSecs(1);
    
    pos1 = gen.mu(1)+[18,-43,-5,-27,70,0,19,2,-18,-66,41,-10,48,-2,-56,-32,16,-20,9,3];  % sample positions
    [XstimIn,YstimIn] = pol2cart(deg2rad(pos1+90),stim.s_r);  % get XY coords relative to origin of zero
    YstimIn = YstimIn.*(1-stim.yscaling);  % decrease y-offset of dots by multiplicative scaling factor
    XstimIn = XstimIn+ X; YstimIn = YstimIn+Y;  % reference to screen center
    g_ps = OffsetRect(gabor.rect, XstimIn', YstimIn')';   % get positions of checkerboards
    vbl = GetSecs;  % get start time
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    for s = 1:length(pos1)
        % Sample on
        Screen('FillOval', window, cs, ps);   % drawing fixation point
        Screen('DrawTextures', window, gabor_ex .tex(1), [], g_ps(:,1:s));  % drawing checkerboards
        Screen('DrawTexture', window, stim.cbartex1);  % drawing colourbar
        Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
        DrawFormattedText(window,str1,'center',options.resolution(2)*0.21-20,[255 255 255]);
        DrawFormattedText(window,str2,'center',options.resolution(2)*0.21+20,[255 255 255]);
        vbl = Screen('Flip', window, vbl + timing.s_ex  + timing.sgap - timing.ifi*0.5);  % set next dot to be flipped at specified dot presentation rate
    end
    
    str3 = 'Note that most of the drawn positions were LEFT of the\nmid-line and the most common position from this deck is around HERE...';
    str4 = '...BUT, some of the drawn positions from this deck\ncan also be to the right of the mid-line.';
    str5 = 'It''s very unlikely that\npositions around here will be\ndrawn from this deck.';
    str6 = 'Press any key to see more information...';
    
    Screen('DrawTexture', window, stim.cbartex1);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('DrawTextures', window, gabor_ex .tex(1), [], g_ps(:,1:s));  % drawing checkerboards
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.21-20,[255 255 255]);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.21+20,[255 255 255]);
    DrawFormattedText(window,str6,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('TextSize', window, 15);
    DrawFormattedText(window,str3,X*0.46,options.resolution(2)*0.833,[255 255 255],[],[],[],2);
    Screen('DrawLines', window, [X-X*0.075 X-X*0.085; Y+Y*0.49 Y+Y*0.70], 1, [255 255 255]);   % drawing pointer
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    Screen('DrawTexture', window, stim.cbartex1);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('DrawTextures', window, gabor_ex .tex(1), [], g_ps(:,1:s));  % drawing checkerboards
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.21-20,[255 255 255]);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.21+20,[255 255 255]);
    DrawFormattedText(window,str6,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('TextSize', window, 15);
    DrawFormattedText(window,str3,X*0.49,options.resolution(2)*0.833,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str4,X*1.12,options.resolution(2)*0.833,[255 255 255],[],[],[],2);
    Screen('DrawLines', window, [X-X*0.075 X-X*0.085 X+X*0.03 X+X*0.14 X+X*0.14 X+X*0.22; Y+Y*0.49 Y+Y*0.70 Y+Y*0.50 Y+Y*0.64 Y+Y*0.64 Y+Y*0.42], 1, [255 255 255]);   % drawing pointer
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    str6 = 'Press any key to see information about DECK 2...';
    
    Screen('DrawTexture', window, stim.cbartex1);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('DrawTextures', window, gabor_ex .tex(1), [], g_ps(:,1:s));  % drawing checkerboards
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.21-20,[255 255 255]);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.21+20,[255 255 255]);
    DrawFormattedText(window,str6,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('TextSize', window, 15);
    DrawFormattedText(window,str3,X*0.49,options.resolution(2)*0.833,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str4,X*1.12,options.resolution(2)*0.833,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str5,X*1.48,options.resolution(2)*0.48,[255 255 255],[],[],[],2);
    Screen('DrawLines', window, [X-X*0.075 X-X*0.085 X+X*0.03 X+X*0.14 X+X*0.14 X+X*0.22; Y+Y*0.49 Y+Y*0.70 Y+Y*0.50 Y+Y*0.64 Y+Y*0.64 Y+Y*0.42], 1, [255 255 255]);   % drawing pointer
    Screen('DrawLines', window, [X*1.4 X*1.47; Y*1.061 Y*1.036], 1, [255 255 255]);   % drawing pointer
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    %% Present dist 2
    str1 = 'Here''s what the distribution of positions from DECK 2 looks like.';
    str2 = 'Here, brighter red colors indicate the most common positions contained in this deck:';
    str4 = 'Press any key to draw some positions from DECK 2...';
    
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    Screen('DrawTexture', window, stim.cbartex2);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.21-20,[255 255 255]);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.21+20,[255 255 255]);
    DrawFormattedText(window,str4,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    WaitSecs(0.6);
    
    pos2 = gen.mu(2)+[18,-43,-5,-27,70,0,19,2,-18,-66,41,-10,48,-2,-56,-32,16,-20,9,3].*-1;  % sample positions
    [XstimIn,YstimIn] = pol2cart(deg2rad(pos2+90),stim.s_r);  % get XY coords relative to origin of zero
    YstimIn = YstimIn.*(1-stim.yscaling);  % decrease y-offset of dots by multiplicative scaling factor
    XstimIn = XstimIn+ X; YstimIn = YstimIn+Y;  % reference to screen center
    g_ps = OffsetRect(gabor.rect, XstimIn', YstimIn')';   % get positions of checkerboards
    vbl = GetSecs;  % get start time
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    for s = 1:length(pos1)
        % Sample on
        Screen('FillOval', window, cs, ps);   % drawing fixation point
        Screen('DrawTextures', window, gabor_ex .tex(1), [], g_ps(:,1:s));  % drawing checkerboards
        Screen('DrawTexture', window, stim.cbartex2);  % drawing colourbar
        Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
        DrawFormattedText(window,str1,'center',options.resolution(2)*0.21-20,[255 255 255]);
        DrawFormattedText(window,str2,'center',options.resolution(2)*0.21+20,[255 255 255]);
        vbl = Screen('Flip', window, vbl + timing.s_ex  + timing.sgap - timing.ifi*0.5);  % set next dot to be flipped at specified dot presentation rate
    end
    
    WaitSecs(1);
    
    str3 = 'In this deck, HERE is the most common position\nand most of the drawn positions were RIGHT of the mid-line.';
    str4 = 'Some of the drawn positions from this deck\ncan also be to the left of the mid-line.';
    str5 = 'It''s very unlikely that\npositions around here will be\ndrawn from this deck.';
    str6 = 'Press any key to continue...';
    
    Screen('DrawTexture', window, stim.cbartex2);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('DrawTextures', window, gabor_ex .tex(1), [], g_ps);  % drawing checkerboards
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.21-20,[255 255 255]);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.21+20,[255 255 255]);
    DrawFormattedText(window,str6,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('TextSize', window, 15);
    DrawFormattedText(window,str3,X*1.09,options.resolution(2)*0.833,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str4,X*0.59,options.resolution(2)*0.833,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str5,X*0.31,options.resolution(2)*0.5,[255 255 255],[],[],[],2);
    Screen('DrawLines', window, [X+X*0.075 X+X*0.19 X-X*0.03 X-X*0.14 X-X*0.14 X-X*0.22; Y+Y*0.49 Y+Y*0.64 Y+Y*0.50 Y+Y*0.64 Y+Y*0.64 Y+Y*0.42], 1, [255 255 255]);   % drawing pointer
    Screen('DrawLines', window, [X*0.52 X*0.59; Y*1.061 Y*1.036], 1, [255 255 255]);   % drawing pointer
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    %% Present full colorbar indicating LLRs
    str1 = 'As you can see, there''s quite a lot of overlap between the two decks,\nwhich can make it difficult to decide which deck is providing dot\npositions at a given point in time.';
    str2 = ['To help you, we will provide a visual aid in the form of the color bar below.\n'...
        'This indicates **how likely each dot position is to come from one deck compared to the other**:\n'...
        'More GREEN indicates more likely from DECK 1, more PINK/RED indicates more likely from DECK 2.\nTake a close look.'];
    str3 = '*Very* likely from\nDECK 1';
    str4 = 'Somewhat more likely\nfrom DECK 1';
    str5 = 'Somewhat more likely\nfrom DECK 2';
    str6 = '*Very* likely from\nDECK 2';
    str7 = 'Equally likely\nfrom both decks';
    str8 = 'Press any key to continue...';
    
    Screen('DrawTexture', window, stim.cbartex);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.25-185,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.25-45,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str8,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('TextSize', window, 15);
    DrawFormattedText(window,str3,X*0.49,options.resolution(2)*0.5,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str4,X*0.61,options.resolution(2)*0.71,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str5,X*1.21,options.resolution(2)*0.71,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str6,X*1.38,options.resolution(2)*0.5,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str7,'center',options.resolution(2)*0.77,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    
    str1 = ['This color bar can really help you to decide which deck is being used.\n \n',...
        'For example, if a dot appears toward the very left (the most green location), you can\n',...
        'be confident that DECK 1 is being used to generate the dot positions.\n \n',...
        'Similarly, if a dot appears toward the very right (the most pink/red location),\n',...
        'you can be confident that the dot positions are being drawn from DECK 2.\n \n',...
        'On the other hand, dot positions close to the middle are equally likely\nto come from either deck.'];
    
    Screen('DrawTexture', window, stim.cbartex);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [180 180 180]);   % drawing vertical midline
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.1,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str8,'center',options.resolution(2)*0.93,[255 255 255]);
    Screen('TextSize', window, 15);
    DrawFormattedText(window,str3,X*0.49,options.resolution(2)*0.5,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str4,X*0.61,options.resolution(2)*0.71,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str5,X*1.21,options.resolution(2)*0.71,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str6,X*1.38,options.resolution(2)*0.5,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str7,'center',options.resolution(2)*0.77,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    %% Present trials with NO change-points and both auditory and visual feedback
    str1 = ['Now let''s do some practice trials. Remember, your task is to decide whether\nthe dot positions you see are being drawn from DECK 1 or DECK 2.\n \n \n',...
        'From now on, each dot will not stay on the screen, but will disappear shortly after it is shown.\n \n \n'...
        'When the dot stops changing position and the fixation point turns green, you should make a response.\n',...
        'If you think the positions were drawn from DECK 1, press the LEFT STRG key;\nif you think they were drawn from DECK 2, press the RIGHT STRG key.'];
    str6 = 'Press any key to continue...';
    
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.24,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str6,'center',options.resolution(2)*0.93,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    str1 = ['Remember, it''s very important that you look at the central fixation point\nat all times, and do not look at the individual dots themselves.\n \n \n',...
        'Also, you should only blink after making a response on each trial,\nand not while the sequence of dots is being presented.'];
    str6 = 'Press any key to start the practice trials...';
    
    Screen('TextSize', window, 18);
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.32,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str6,'center',options.resolution(2)*0.93,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    
    run_prac = 'y';
    while  sum(strcmp('y',run_prac))>0
        % Generate sequences of dot positions
        fdists = [ones(ceil(ntrials_noCP/2),1); ones(ceil(ntrials_noCP/2),1).*2];
        stimIn = [round(gen.mu(1)+(randn(ceil(ntrials_noCP/2),maxsamps).*gen.sigma(1))); round(gen.mu(2)+(randn(ceil(ntrials_noCP/2),maxsamps).*gen.sigma(2)))];
        stimIn(stimIn<gen.range(1)) = gen.range(1); stimIn(stimIn>gen.range(2)) = gen.range(2);  % in case drawn values exceed range limits
        
        shuforder = randperm(ntrials_noCP); % shuffle trial order
        stimIn = stimIn(shuforder,:);
        fdists = fdists(shuforder);
        
        tlengths = randsample(2:(maxsamps-1),round(ntrials_noCP*pshort_noCP),true);  % truncate sequence length of random selection of trials
        ttrials = sort(randsample(3:ntrials_noCP,round(ntrials_noCP*pshort_noCP),false))';  % picking random selection of trials to be truncated - makes sure first 2 trials are always full-length
        for t = 1:length(ttrials)
            stimIn(ttrials(t),tlengths(t)+1:end) = nan;
        end
        
        [XstimIn,YstimIn] = pol2cart(deg2rad(stimIn+90),stim.s_r);  % transform dot positions from polar to cartesian (x,y) coordinates
        YstimIn = YstimIn.*(1-stim.yscaling);  % decrease y-offset of dots by multiplicative scaling factor
        XstimIn = XstimIn+ X; YstimIn = YstimIn+Y;  % reference to screen center
        
        % Countdown and first fixation point
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
        Screen('FillOval', window, cs, ps);   % drawing fixation point
        vbl = Screen('Flip', window);   % then fixation plus patches
        
        % Loop through trials
        on_time = [];  % onset time for first trial is specified within trial script
        Behav = zeros(length(fdists),5); tRefresh = zeros(length(fdists),maxsamps);
        for t = 1:length(fdists);
            % Variable task parameters
            varopts.on_time = on_time;                        % controls onset time of impending trial - fed back from trial function
            varopts.xpos = XstimIn(t,~isnan(XstimIn(t,:)));   % sequence of sample X-coordinates
            varopts.ypos = YstimIn(t,~isnan(YstimIn(t,:)));   % sequence of sample Y-coordinates
            varopts.CP = nan(size(varopts.xpos));             % flag for change-points
            varopts.fdist = fdists(t);                        % correct answer (i.e. generative distribution @ sequence end: 1=left, 2=right)
            varopts.kbqdev = options.kb_setup;                % keyboard info
            
            % Run trial
            [Behav(t,1:5),tRefresh(t,1:length(varopts.xpos)),on_time] = Surprise_radial_checkers_trial_visFB(window, windowRect, timing, stim, gabor, varopts, trigger_enc, t, 0);
        end
        
        % Present average accuracy and ask whether practice should be repeated
        Screen('TextSize', window, 18);
        str1 = 'Well done!';
        str2 = sprintf('Your average accuracy on those practice trials was %2.1f percent.',nansum(Behav(:,3))/size(Behav,1)*100);
        str3 = 'Would you like to try another round of practice trials? (y/n)';
        
        DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
        DrawFormattedText(window,str2,'center',Y,[255 255 255]);
        DrawFormattedText(window,str3,'center',Y+120,[255 255 255]);
        Screen('Flip', window);
        WaitSecs(1);
        
        [~, keyCode] = KbWait;
        run_prac = KbName(keyCode);
        WaitSecs(1);
    end
    
    
    %% Present trials WITH change-points and both auditory and visual feedback
    str1 = ['Very good!\n \n \n',...
        'Now, we''ll move onto one of the most important features of the task. Specifically, from now on,\n',...
        'the deck that is used to draw dot positions can **change at any time within a sequence**.\n \n \n',...
        'Your task from now on is to decide which deck is being used to generate\nthe dot positions **AT THE END** of each sequence of dots. For example, if a sequence starts with\n',...
        'DECK 1, changes after some time to DECK 2, and then ends, the correct answer is DECK 2.'];
    str6 = 'Press any key for more instructions...';
    
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.28,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str6,'center',options.resolution(2)*0.93,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    str1 = ['The chance that the deck will change between any two consecutive dot positions is about 1 in 12.\n \n',...
        'This means that, on average, there will be about 1 deck change per trial for longer\n',...
        'dot sequences. However, it is also possible for there to be no deck change at all on a given trial,\n',...
        'or for there to be several deck changes.\n \n',...
        'Deck changes can occur at any point during of a sequence of dots.'];
    
    str6 = 'Press any key for more instructions...';
    
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.32,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str6,'center',options.resolution(2)*0.93,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    
    str1 = 'Let''s do some more practice trials, this time with possible changes in deck during each dot sequence.';
    str2 = ['To help you get to grips with this very important feature of the task, for these practice trials\nwe will indicate (with the word ''CHANGE!'') whenever a deck change occurs.\n'...
        '(During non-practice trials, you will not be given this obvious indication of a change.)'];
    str3 = 'Remember to always look at the central fixation point, and only blink after making a response.';
    str4 = 'Press any key to begin...';
    
    DrawFormattedText(window,str1,'center',options.resolution(2)*0.28,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str2,'center',options.resolution(2)*0.40,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str3,'center',options.resolution(2)*0.57,[255 255 255],[],[],[],2);
    DrawFormattedText(window,str4,'center',options.resolution(2)*0.93,[255 255 255],[],[],[],2);
    Screen('Flip', window);
    
    WaitSecs(1);
    KbWait;
    
    run_prac = 'y';
    while  sum(strcmp('y',run_prac))>0
        % Generate sequences of distribution switch positions
        pswitch = [ones(ntrials_CP,1) rand(ntrials_CP,maxsamps-1)];  % first draw uniformly distributed random probabilities b/w 0 and 1 that will determine switch positions (first sample is never a switch)
        pswitch(pswitch>gen.H) = 0; pswitch(pswitch~=0) = 1;  % binarize matrix to mark only switch positions
        
        pswitch(1:2,:) = 0; pswitch(1,9) = 1; pswitch(2,5) = 1;   % enforcing one, and only one, switch in first two trials
        
        % Generate sequences of which distributions will be drawn from at which times (1=left, 2=right)
        distseqs = zeros(ntrials_CP,maxsamps); dists = [1 2];
        for t = 1:ntrials_CP
            if t==2, cdist = 2;  % enforcing first two trials to have different correct responses
            elseif t<=ntrials_CP/2, cdist = 1; else cdist = 2; end  % making sure each distribution is starting distribution an equal number of times
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
        
        % Shuffle trial order (except for first two trials)
        shuforder = randperm(ntrials_CP); shuforder = [1 2 shuforder(shuforder>2)];
        pswitch = pswitch(shuforder,:);
        stimIn = stimIn(shuforder,:);
        distseqs = distseqs(shuforder,:);
        
        % Store generative distributions @ sequence end (i.e. correct choices; 1=left, 2=right)
        fdists = distseqs(:,end);
        
        % Truncate sequence length of random selection of trials
        tlengths = randsample(2:(maxsamps-1),round(ntrials_CP*pshort_CP),true);  % randomly draw lengths of truncated sequences
        ttrials = sort(randsample(3:ntrials_CP,round(ntrials_CP*pshort_CP),false))';  % randomly draw trials to be truncated (making sure 1st two trials are full-length)
        for t = 1:length(ttrials)
            pswitch(ttrials(t),tlengths(t)+1:end) = nan;
            stimIn(ttrials(t),tlengths(t)+1:end) = nan;
            distseqs(ttrials(t),tlengths(t)+1:end) = nan;
            fdists(ttrials(t)) = distseqs(ttrials(t),tlengths(t));
        end
        
        [XstimIn,YstimIn] = pol2cart(deg2rad(stimIn+90),stim.s_r);  % transform dot positions from polar to cartesian (x,y) coordinates
        YstimIn = YstimIn.*(1-stim.yscaling);  % decrease y-offset of dots by multiplicative scaling factor
        XstimIn = XstimIn+ X; YstimIn = YstimIn+Y;  % reference to screen center
        
        % Countdown and first fixation point
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
        Screen('FillOval', window, cs, ps);   % drawing fixation point
        vbl = Screen('Flip', window);   % then fixation plus patches
        
        % Loop through trials
        on_time = [];  % onset time for first trial is specified within trial script
        Behav = zeros(length(fdists),5); tRefresh = zeros(length(fdists),maxsamps);
        for t = 1:length(fdists);
            % Variable task parameters
            varopts.on_time = on_time;                        % controls onset time of impending trial - fed back from trial function
            varopts.xpos = XstimIn(t,~isnan(XstimIn(t,:)));   % sequence of sample X-coordinates
            varopts.ypos = YstimIn(t,~isnan(YstimIn(t,:)));   % sequence of sample Y-coordinates
            varopts.CP = pswitch(t,:);                        % flag for change-points
            varopts.fdist = fdists(t);                        % correct answer (i.e. generative distribution @ sequence end: 1=left, 2=right)
            varopts.kbqdev = options.kb_setup;                % keyboard info
            
            % Run trial
            [Behav(t,1:5),tRefresh(t,1:length(varopts.xpos)),on_time] = Surprise_radial_checkers_trial_visFB(window, windowRect, timing, stim, gabor, varopts, trigger_enc, t, 0);
        end
        
        % Present average accuracy and ask whether practice should be repeated
        Screen('TextSize', window, 18);
        str1 = 'Well done!';
        str2 = sprintf('Your average accuracy on those practice trials was %2.1f percent.',nansum(Behav(:,3))/size(Behav,1)*100);
        str3 = 'Would you like to try another round of practice trials? (y/n)';
        
        DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
        DrawFormattedText(window,str2,'center',Y,[255 255 255]);
        DrawFormattedText(window,str3,'center',Y+120,[255 255 255]);
        Screen('Flip', window);
        WaitSecs(1);
        
        [~, keyCode] = KbWait;
        run_prac = KbName(keyCode);
        WaitSecs(1);
    end
    
    %% Finish up
    str1 = 'Great!';
    str2 = 'The practice is now over.';
    
    DrawFormattedText(window,str1,'center',Y-60,[255 255 255]);
    DrawFormattedText(window,str2,'center',Y+60,[255 255 255]);
    Screen('Flip', window);
    
    WaitSecs(4.5);
    
    %% Exit
    FlushEvents('keyDown');
    ListenChar(0)  % reinstate keyboard input to command line
    PsychPortAudio('DeleteBuffer');
    PsychPortAudio('Close',stim.audio.h);
    sca
    
    
catch ME  % if any errors are encountered or user requests quit, clean up and exit
    sca
    ListenChar(0)
    PsychPortAudio('DeleteBuffer');
    PsychPortAudio('Close',stim.audio.h);
    rethrow(ME)
end



