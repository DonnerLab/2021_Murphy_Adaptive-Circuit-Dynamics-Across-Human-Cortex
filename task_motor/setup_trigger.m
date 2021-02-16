function trig = setup_trigger()

trig.zero = 0;
trig.width = 0.005; 

trig.block_start = 1;     % start of block
trig.block_end = 2;       % end of block

trig.rest_on = 11;        % onset of rest fixation
trig.fix_on = 21;         % onset of active fixation
trig.resp_cue_on = 31;    % onset of response type cue
trig.go_cue_on = 41;      % onset of go cue

trig.respL = 41;          % 'left' response
trig.respR = 42;          % 'right' response
trig.resp_bad = 99;       % bad response

end
