function data = asc2dat(asc)
% takes asc data from EyeLink file and converts this into events and data structure

% make data struct
data                = [];

% create event structure for messages
evcell = cell(length(asc.msg),1);
event = struct('type', evcell, 'sample', evcell, 'value', evcell, 'offset', evcell, 'duration', evcell );

for i=1:length(asc.msg),
    
    strtok = tokenize(asc.msg{i});
    event(i).type = strtok{3};
    str2double(strtok{2});
    
    % match the message to its sample
    smpstamp = dsearchn(asc.dat(1,:)', str2double(strtok{2}));
    % find closest sample index of trigger in ascii dat
    
    if ~isempty(smpstamp)
        event(i).sample = smpstamp(1);
    else % if no exact sample was found
        warning('no sample found');
    end
    event(i).value = asc.msg{i};
end

data.event = event;  % add events into data structure

% add data into data structure
% important: match the right data chans to their corresponding labels...
data.Xgaze          = asc.dat(2, :);
data.Ygaze          = asc.dat(3, :);
data.pupil          = asc.dat(4, :);
data.fsample        = asc.fsample;
data.times          = 0:1/data.fsample:length(asc.dat(1,:))/data.fsample-1/data.fsample;
data.sampleinfo     = [1 length(asc.dat(1,:))];

% if data.fsample ~= 1000,
%     warning('pupil not sampled with 1000Hz');
% end

% parse blinks
if ~isempty(asc.eblink)
    blinktimes = cellfun(@regexp, asc.eblink, ...
        repmat({'\d*'}, length(asc.eblink), 1), repmat({'match'}, length(asc.eblink), 1), ...
        'UniformOutput', false); % parse blinktimes from ascdat
    blinktimes2 = nan(length(blinktimes), 2);
    for s = 1:length(blinktimes), a = blinktimes{s};
        for j = 1:2, blinktimes2(s, j) = str2double(a{j}); end
    end
    timestamps = asc.dat(1,:); % get the time info
    try
        blinksmp = arrayfun(@(x) find(timestamps == x, 1,'first'), blinktimes2, 'UniformOutput', true ); %find sample indices of blinktimes in timestamps
    catch
        blinksmp = arrayfun(@(x) dsearchn(timestamps', x), blinktimes2, 'UniformOutput', true ); %find sample indices of blinktimes in timestamps
    end
    data.blinksmp = blinksmp;   % add blink samples into data structure
else
    data.blinksmp = [];
end

% parse saccades
if ~isempty(asc.esacc)
    sacctimes = cellfun(@regexp, asc.esacc, ...
        repmat({'\d*'}, length(asc.esacc), 1), repmat({'match'}, length(asc.esacc), 1), ...
        'UniformOutput', false); % parse blinktimes from ascdat
    sacctimes2 = nan(length(sacctimes), 2);
    for s = 1:length(sacctimes), a = sacctimes{s};
        for j = 1:2,
            if str2double(a{j}) ~= 0,
                sacctimes2(s, j) = str2double(a{j});
            else
                sacctimes2(s, j) = str2double(a{j+1});
            end
        end
    end
    timestamps = asc.dat(1,:); % get the time info
    try
        saccsmp = arrayfun(@(x) find(timestamps == x, 1,'first'), sacctimes2, 'UniformOutput', true ); %find sample indices
    catch
        saccsmp = arrayfun(@(x) dsearchn(timestamps', x), sacctimes2, 'UniformOutput', true );
    end
    data.saccsmp = saccsmp;   % add saccade  samples into data structure
else
    data.saccsmp = [];
end

end