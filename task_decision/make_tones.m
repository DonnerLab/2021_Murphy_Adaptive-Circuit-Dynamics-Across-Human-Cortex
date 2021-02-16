function tones = make_tones(freqs,durations)

Fs = 8192;   % sampling frequency

for t = 1:length(freqs)
    if ~isnan(freqs(t))   % for sine waves
        values = 0:1/Fs:durations(t);
        tones{t} = sin(2*pi*freqs(t)*values);
    else               % for white noise
        tones{t} = (rand(Fs,1)*2)-1;
    end
end