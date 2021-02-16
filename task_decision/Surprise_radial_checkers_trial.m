% Presents one trial of luminance patch signal detection task.
%
% Two circular patches are presented in the lower hemifield with luminances
% that fluctuate independently over time. At some point during the trial,
% the luminance of one patch may increase (on average) for some
% pre-specified duration. Trial ends with written feedback on accuracy.
%
% Outputs:
%    Behav = [signalPos signalDur signalOn signalStr RTabs RTrel ACCabs ACCdet start_time]
%    tRefresh = timing of all flips during trial
%    onset_nplus1 = specifies time for next stimulus presentation

function [Behav,tRefresh,onset_nplus1] = Surprise_radial_checkers_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, t, EL)

% Response parameters
respL = '1';
respR = '4';
quitkey = 'ESCAPE';

% flush_kbqueues(varopts.kbqdev);
FlushEvents('keyDown');    %%%%%%%%%%% PM: this may mess things up......

% Get pixel coordinate of screen center
[X,Y] = RectCenter(windowRect);

% Getting eye that's being measured (matters for retrieving online gaze data)
if EL, eyeused = Eyelink('EyeAvailable'); end
sacc=[]; % and initializing saccade counter

% Fixation point positions and colours
pos_fix_in = [X-stim.fix_in_r; Y-stim.fix_in_r; X+stim.fix_in_r; Y+stim.fix_in_r];
pos_fix_out = [X-stim.fix_out_r; Y-stim.fix_out_r; X+stim.fix_out_r; Y+stim.fix_out_r];
ps = [pos_fix_out pos_fix_in];  % positions
cs = [stim.fix_out_c' stim.fix_in_c'];  % colours

% Centering vertical midline
stim.mid_xy(1,:) = stim.mid_xy(1,:)+X;
stim.mid_xy(2,:) = stim.mid_xy(2,:)+Y;

% Present initial static checkerboard & rest fixation before the first trial of a block (not constructed for any other trial because it will carry over from previous iteration)
if t==1
    cs = [stim.fix_out_c_rest' stim.fix_in_c'];  % colours
    g_ps = OffsetRect(gabor.rect, X, Y+(stim.s_r.*(1-stim.yscaling)));   % specify position of central gabor @ block onset
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');   % these blend functions are necessary for transparency in the colourbar/checkerboard textures
    Screen('DrawTexture', window, gabor.tex(2), [], g_ps);  % drawing checkerboard
    Screen('DrawTexture', window, stim.cbartex);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);   % drawing vertical midline
    vbl = Screen('Flip', window);   % present fixation plus mask
    trigger(trigger_enc.block_start);  % trigger to mark onset of baseline period
    if EL, Eyelink('message', num2str(trigger_enc.block_start)); [~,~,~] = Eyelink('GetQueuedData'); end  % sending EL message for block start, and flushing queue
    
    cs = [stim.fix_out_c' stim.fix_in_c'];  % switch to active fixation
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('DrawTexture', window, gabor.tex(2), [], g_ps);  % drawing checkerboard
    Screen('DrawTexture', window, stim.cbartex);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);   % drawing vertical midline
    vbl = Screen('Flip', window, vbl+4);   % present fixation plus mask
    
    varopts.on_time = vbl + timing.ITI - timing.rest + rand(1)*timing.ITIdist - timing.fb*2;  % drawing random jitter for timing of first-trial onset relative to active fixation onset
end

if EL, sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end

% Present forward mask: central checkerboard with flipped tiling, then briefly remove
g_ps = OffsetRect(gabor.rect, X, Y+(stim.s_r*(1-stim.yscaling)));   % specify position of forward mask gabor
Screen('FillOval', window, cs, ps);   % drawing fixation point
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawTexture', window, gabor.tex(1), [], g_ps);  % drawing checkerboard
Screen('DrawTexture', window, stim.cbartex);
Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);
start_time = Screen('Flip', window, varopts.on_time);  % set forward mask to be flipped at specified trial onset time
trigger(trigger_enc.premask_on);
if EL, Eyelink('message', num2str(trigger_enc.premask_on)); sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end

cmask = 2;
for f = 1:ceil((timing.s/timing.s_refresh)-1)  % flicker forward mask
    Screen('FillOval', window, cs, ps);
    Screen('DrawTexture', window, gabor.tex(cmask), [], g_ps);
    Screen('DrawTexture', window, stim.cbartex);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);
    Screen('Flip', window, start_time+(timing.s_refresh*f)-timing.ifi*0.5);
    cmask = find([1 2]~=cmask);
    if EL, sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end
end

Screen('FillOval', window, cs, ps);   % drawing fixation/color bar/midline only
Screen('DrawTexture', window, stim.cbartex);
Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);
vbl = Screen('Flip', window, start_time + timing.s - timing.ifi*0.5);
trigger(trigger_enc.premask_off);
if EL, Eyelink('message', num2str(trigger_enc.premask_off)); end

% Present sequence of evidence samples
tRefresh = nan(1,length(varopts.xpos));

for s = 1:length(varopts.xpos)
    % Sample on
    g_ps = OffsetRect(gabor.rect, varopts.xpos(s), varopts.ypos(s));   % set new position of Gabor
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    Screen('DrawTexture', window, gabor.tex(1), [], g_ps);  % drawing checkerboard
    Screen('DrawTexture', window, stim.cbartex);  % drawing colourbar
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);   % drawing vertical midline
    vbl = Screen('Flip', window, vbl + timing.sgap - timing.ifi*0.5);  % set next dot to be flipped at specified dot presentation rate
    trigger(trigger_enc.sample_on);
    if EL, Eyelink('message', num2str(trigger_enc.sample_on)); sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end
    tRefresh(s) = vbl-start_time;   % logging sample onset times to verify timing consistency
    
    cmask = 2;
    for f = 1:ceil((timing.s/timing.s_refresh)-1)  % flicker forward mask
        Screen('FillOval', window, cs, ps);
        Screen('DrawTexture', window, gabor.tex(cmask), [], g_ps);
        Screen('DrawTexture', window, stim.cbartex);
        Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);
        Screen('Flip', window, vbl+(timing.s_refresh*f)-timing.ifi*0.5);
        cmask = find([1 2]~=cmask);
        if EL, sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end
    end
    
    % Sample off
    Screen('FillOval', window, cs, ps);   % drawing fixation/color bar/midline only
    Screen('DrawTexture', window, stim.cbartex);
    Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);
    vbl = Screen('Flip', window, vbl + timing.s - timing.ifi*0.5);  % remove gabor at desired time
    trigger(trigger_enc.sample_off);
    if EL, Eyelink('message', num2str(trigger_enc.sample_off)); end
end

% Present backward mask: central checkerboard with inverse tiling
g_ps = OffsetRect(gabor.rect, X, Y+(stim.s_r*(1-stim.yscaling)));   % specify position of backward mask gabor
Screen('FillOval', window, cs, ps);
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawTexture', window, gabor.tex(2), [], g_ps);  % drawing checkerboard
Screen('DrawTexture', window, stim.cbartex);
Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);
vbl = Screen('Flip', window, vbl + timing.sgap - timing.ifi*0.5);  % set backward mask to be flipped at desired time
trigger(trigger_enc.postmask);
if EL, Eyelink('message', num2str(trigger_enc.postmask)); sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end

% Prepare for response checking
%flush_kbqueues(varopts.kbqdev);
FlushEvents('keyDown')
post_int = timing.post + rand(1)*timing.postdist;  % drawing pre-response cue interval from specified uniform distribution
keyIsDown = false;

% Present response cue: change in outer fixation point color
cs(:,1) = stim.fix_out_c_resp';   % change colour of outer fixation point
Screen('FillOval', window, cs, ps);
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawTexture', window, gabor.tex(2), [], g_ps);  % drawing checkerboard
Screen('DrawTexture', window, stim.cbartex);
Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);
vbl = Screen('Flip', window, vbl + post_int - timing.ifi*0.5);  % flip response cue after jittered interval
trigger(trigger_enc.resp_cue);
if EL, Eyelink('message', num2str(trigger_enc.resp_cue)); sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end

% Wait for response
while ~keyIsDown
    %[keyIsDown, firstPress] = check_kbqueues(kbqdev);
    [keyIsDown, secs, firstPress] = KbCheck;
    
    if keyIsDown  % logging response type
        
        RT = secs-vbl;  % logging RT relative to response cue onset
        keys = KbName(firstPress);  % retrieving string variable containing currently pressed key(s)
        
        if iscell(keys)
            resp = 99;  % in case of a double-press...having this as first 'if' test means it takes absolute precedence
            trigger(trigger_enc.resp_bad);  % trigger to mark a bad response
            if EL, Eyelink('message', num2str(trigger_enc.resp_bad)); end
        else
            switch keys
                case quitkey  % user requests quit
                    throw(MException('EXP:Quit', 'User request quit'));
                case {respL, 'LeftArrow', '1!', 'LeftControl'}
                    resp = 1;
                    trigger(trigger_enc.resp_left);  % trigger to mark a left response
                    if EL, Eyelink('message', num2str(trigger_enc.resp_left)); end
                case {respR, 'RightArrow', '4$', 'RightControl'}
                    resp = 2;
                    trigger(trigger_enc.resp_right);  % trigger to mark a right response
                    if EL, Eyelink('message', num2str(trigger_enc.resp_right)); end
                otherwise
                    resp = 99;  % in case any button other than task relevant ones is pressed
                    trigger(trigger_enc.resp_bad);  % trigger to mark a bad response
                    if EL, Eyelink('message', num2str(trigger_enc.resp_bad)); end
            end
        end
    end
end

% Classify response accuracy and associated feedback
if resp == 99                            % bad press
    ACC = nan;
    fb_trig = trigger_enc.fb_bad;
elseif resp == varopts.fdist             % correct response
    ACC = 1;
    fb_trig = trigger_enc.fb_correct;
else                                     % incorrect response
    ACC = 0;
    fb_trig = trigger_enc.fb_error;
end

% Present auditory feedback
WaitSecs(timing.prefb-GetSecs-secs);   % waiting short interval before feedback presentation
trigger(fb_trig);  % trigger to mark feedback onset
if EL, Eyelink('message', num2str(fb_trig)); end
if ACC == 1                     % ascending tone for correct response
    PsychPortAudio('SetLoop',stim.audio.h, stim.audio.tonepos(1,1), stim.audio.tonepos(1,2));
    fbtime = PsychPortAudio('Start', stim.audio.h, 1, 0, 1);
    PsychPortAudio('Stop', stim.audio.h, 1);
    PsychPortAudio('SetLoop',stim.audio.h, stim.audio.tonepos(2,1), stim.audio.tonepos(2,2));
    PsychPortAudio('Start', stim.audio.h, 1, fbtime+timing.fb, 1);
elseif ACC == 0                 % descending tone for incorrect response
    PsychPortAudio('SetLoop',stim.audio.h, stim.audio.tonepos(2,1), stim.audio.tonepos(2,2));
    fbtime = PsychPortAudio('Start', stim.audio.h, 1, 0, 1);
    PsychPortAudio('Stop', stim.audio.h, 1);
    PsychPortAudio('SetLoop',stim.audio.h, stim.audio.tonepos(1,1), stim.audio.tonepos(1,2));
    PsychPortAudio('Start', stim.audio.h, 1, fbtime+timing.fb, 1);
else                            % white noise tone for bad response
    PsychPortAudio('SetLoop',stim.audio.h, stim.audio.tonepos(3,1), stim.audio.tonepos(3,2));
    fbtime = PsychPortAudio('Start', stim.audio.h, 1, 0, 1);
end
PsychPortAudio('Stop', stim.audio.h, 1);

% Present rest-period fixation
cs(:,1) = stim.fix_out_rest';   % change colour of outer fixation point
Screen('FillOval', window, cs, ps);
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawTexture', window, gabor.tex(2), [], g_ps);  % drawing checkerboard
Screen('DrawTexture', window, stim.cbartex);
Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);
vbl = Screen('Flip', window, fbtime + (timing.fb*2) - timing.ifi*0.5);

% Flush Eyelink queue
[~,~,~] = Eyelink('GetQueuedData');

% Revert to regular active fixation
cs(:,1) = stim.fix_out_c';   % change colour of outer fixation point
Screen('FillOval', window, cs, ps);
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawTexture', window, gabor.tex(2), [], g_ps);  % drawing checkerboard
Screen('DrawTexture', window, stim.cbartex);
Screen('DrawLines', window, stim.mid_xy, stim.mid_w, [200 200 200]);
vbl = Screen('Flip', window, vbl + timing.rest - timing.ifi*0.5);

% Draw onset time of next trial
onset_nplus1 = vbl + timing.ITI - timing.rest + rand(1)*timing.ITIdist - timing.fb*2;

% Count number of broken fixations
if isempty(sacc), sacc = nan; else sacc = sum(sacc); end

% Concatenate final output variable
Behav = [varopts.fdist resp ACC RT start_time sacc];

end
