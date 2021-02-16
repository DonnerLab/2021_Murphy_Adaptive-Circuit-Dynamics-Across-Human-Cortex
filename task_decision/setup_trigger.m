function trig = setup_trigger()

trig.zero = 0;
trig.width = 0.005; 

trig.block_start = 1;  % start of block
trig.block_end = 2;    % end of block

trig.premask_on = 11;   % onset of pre-sequence mask
trig.premask_off = 12;  % offset of pre-sequence mask

trig.sample_on = 21;   % onset of individual sample
trig.sample_off = 22;  % offset of individual sample

trig.postmask = 31;    % onset of individual sample

trig.resp_cue = 41;     % cue for response
trig.resp_left = 42;    % 'left' response
trig.resp_right = 43;   % 'right' response
trig.resp_bad = 44;     % bad response (either double press, or task-irrelevant button)

trig.fb_correct = 51;   % feedback for hit or correct rejection
trig.fb_error = 52;     % feedback for false alarm, miss, mislocalization, or premature response
trig.fb_bad = 53;       % feedback for bad responses

trig.trial_end = 61;    % offset of feedback/onset of break period

end
