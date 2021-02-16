function audio = setup_audio(tone_freqs,tone_dur)

% Initialize audio structure
audio           = [];

% Initial setup
InitializePsychSound(1);  % request low latency mode

audio.freq      = 44100;
audio.h         = PsychPortAudio('Open', [], 1, [], audio.freq, 2, [], []);  % open default soundport, in stereo (to match the sound matrix we create)

% Create tones & load into audio buffer
[audio.tonebuf, audio.tonepos] = CreateAudioBuffer(CreateTone(tone_freqs(1), tone_dur, audio.freq), ...
    CreateTone(tone_freqs(2), tone_dur, audio.freq), ...
    CreateTone(tone_freqs(3), tone_dur*2, audio.freq));

PsychPortAudio('FillBuffer', audio.h, audio.tonebuf);
end


% Creates an audio buffer
function [buffer,loop] = CreateAudioBuffer(varargin)

try
    buffer = cat(2,varargin{:});

    loop = zeros(length(varargin),2);
    loop(1,:) = [0,size(varargin{1},2)-1];
    
    for i = 2:length(varargin)
        loop(i,:) = loop(i-1,2)+[1,size(varargin{i},2)];
    end % end of the for loop
    
catch % if something went wrong with the CreateAudioBuffer function, do the following
    % ShowCursor;
    disp('CreateAudioBuffer failed!');
    Screen('CloseAll');
    commandwindow;
    
    
end % end of the try-catch function
end % end of the CreateAudioBuffer function
