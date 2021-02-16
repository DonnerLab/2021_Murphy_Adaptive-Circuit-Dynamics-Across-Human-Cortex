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

function [Behav,onset_nplus1] = Motor_localiser_trial(window, windowRect, timing, stim, varopts, trigger_enc, t, EL)

% Response parameters
resptypes = {'LEFT','RIGHT'};
respkeys = {'1','4'};
quitkey = 'ESCAPE';
RT = nan;
kbqdev = varopts.kbqdev;

%flush_kbqueues(kbqdev);
FlushEvents('keyDown')

% Get pixel coordinate of screen center
[X,Y] = RectCenter(windowRect);

% Fixation cross positions and colours
pos_fix_in = [X-stim.fix_in_r; Y-stim.fix_in_r; X+stim.fix_in_r; Y+stim.fix_in_r];
pos_fix_out = [X-stim.fix_out_r; Y-stim.fix_out_r; X+stim.fix_out_r; Y+stim.fix_out_r];
ps = [pos_fix_out pos_fix_in];  % positions
cs = [stim.fix_out_c_rest' stim.fix_in_c'];  % colours

% Present initial fixation point for the first trial of a block (not constructed for any other trial because it will carry over from previous iteration)
if t==1
    Screen('FillOval', window, cs, ps);   % drawing fixation point
    vbl = Screen('Flip', window);   % present fixation plus mask
    varopts.on_time = vbl+2;
end

% Switch to active fixation
cs(:,1) = stim.fix_out_c_prep';   % change colour of outer fixation point
Screen('FillOval', window, cs, ps);   % drawing fixation point
start_time = Screen('Flip', window, varopts.on_time);  % set stimulus to be flipped at specified trial onset time
trigger(trigger_enc.fix_on);
if EL, Eyelink('message', num2str(trigger_enc.fix_on)); end

% Present response type cue after jittered duration
Screen('FillOval', window, cs, ps);
DrawFormattedText(window, resptypes{varopts.cresp}, 'center', Y+stim.cueYoff);
vbl = Screen('Flip', window, start_time + timing.precue + timing.precueD*rand(1) - timing.ifi*0.5);
trigger(trigger_enc.resp_cue_on);
if EL, Eyelink('message', num2str(trigger_enc.resp_cue_on)); end
cue_time = vbl-start_time;

% Remove response type cue
Screen('FillOval', window, cs, ps);
vbl = Screen('Flip', window, vbl + timing.cue - timing.ifi*0.5);

% Prepare for response checking
%flush_kbqueues(kbqdev);
FlushEvents('keyDown')
keyIsDown = false;

% Present go cue
cs(:,1) = stim.fix_out_c_go';   % change colour of outer fixation point
Screen('FillOval', window, cs, ps);
vbl = Screen('Flip', window, vbl + timing.prep - timing.ifi*0.5);
trigger(trigger_enc.go_cue_on);
if EL, Eyelink('message', num2str(trigger_enc.go_cue_on)); end
go_time = vbl-start_time;

% Wait for response
while ~keyIsDown
    %[keyIsDown, firstPress] = check_kbqueues(kbqdev);
    [keyIsDown, secs, firstPress] = KbCheck;
    
    if keyIsDown  % logging response
        RT = GetSecs-start_time-go_time;  % logging RT relative to go cue onset
        keys = KbName(firstPress);  % retrieving string variable containing currently pressed key(s)
        
        if iscell(keys)
            resp = 99;  % in case of a double-press...having this as first 'if' test means it takes absolute precedence
            trigger(trigger_enc.resp_bad);  % trigger to mark a bad response
            if EL, Eyelink('message', num2str(trigger_enc.resp_bad)); end
        else
            switch keys
                case quitkey  % user requests quit
                    sca
                    throw(MException('EXP:Quit', 'User request quit'));
                case {respkeys{1}, 'LeftArrow', '1!', 'LeftControl'}
                    resp = 1;
                    trigger(trigger_enc.respL);  % trigger to mark a left response
                    if EL, Eyelink('message', num2str(trigger_enc.respL)); end
                case {respkeys{2}, 'RightArrow', '4$', 'RightControl'}
                    resp = 2;
                    trigger(trigger_enc.respR);  % trigger to mark a right response
                    if EL, Eyelink('message', num2str(trigger_enc.respR)); end
                otherwise
                    resp = 99;  % in case any button other than task relevant ones is pressed
                    trigger(trigger_enc.resp_bad);  % trigger to mark a bad response
                    if EL, Eyelink('message', num2str(trigger_enc.resp_bad)); end
            end
        end
        
    end
end

%flush_kbqueues(kbqdev);
FlushEvents('keyDown')

% Switch back to rest fixation point
cs(:,1) = stim.fix_out_c_rest';   % change colour of outer fixation point
Screen('FillOval', window, cs, ps);   % drawing fixation point
vbl = Screen('Flip', window);
trigger(trigger_enc.rest_on);
if EL, Eyelink('message', num2str(trigger_enc.rest_on)); end

% Draw onset time of next trial
onset_nplus1 = vbl + timing.rest - timing.ifi*0.5;

% Concatenate final output variable: [desired_response executed_response timing(1:3) RT]
Behav = [varopts.cresp resp start_time cue_time go_time RT];

end
