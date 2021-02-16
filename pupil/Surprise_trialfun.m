 % Creates pupil trial matrix of the following form:
 %    [nsamps stimonsmp responsmp respsmp respdir ACC RT[s] block_num trial_num]
 % and another matrix containing sample onset times.

function data = Surprise_trialfun(data,b)

event = data.event;
respdirs = [1 2 nan];  % left, right, bad resp
accs = [1 0 nan];  % correct, error, bad resp

% Pull indices of stimulus onset events
trialinds = find(cellfun(@(x) strcmp(x,'TRIALID'),{event.type}));  % TRIALID marks start of each trial

% Loop through each trial
trl = zeros(length(trialinds),9);  % initializing trial content matrix
smp = nan(length(trialinds),12);  % initializing sample times matrix
for i = 1:length(trialinds);
    
    j = trialinds(i);
    
    deets = tokenize(event(j).value);   % parsing event details
    trial_num = str2double(deets{find(strcmp(deets,'TRIALID'))+1});  % TRIAL NUMBER
    
    k = j+1;  % event counter
    s = 1;    % sample counter
    while sum(strcmp(event(k).type,{'TRIALID','2'}))==0
        if strcmp(event(k).type,'11')
            stimonsmp = event(k).sample;  % pre-mask onset time
        elseif strcmp(event(k).type,'21')
            smp(i,s) = event(k).sample;  % sample onset time
            s = s+1;  % increment sample counter
        elseif strcmp(event(k).type,'41')
            responsmp = event(k).sample;  % response cue onset time
        elseif sum(strcmp(event(k).type,{'42','43','44'}))>0
            respsmp = event(k).sample;  % response time
            respdir = respdirs(strcmp(event(k).type,{'42','43','44'}));  % response type
        elseif sum(strcmp(event(k).type,{'51','52','53'}))>0
            ACC = accs(strcmp(event(k).type,{'51','52','53'}));
        end
        k = k+1;
    end
    
    % Add to trial matrix
    RT = (respsmp-responsmp)/data.fsample;
    trl(i,:) = [s-1 stimonsmp responsmp respsmp respdir ACC RT b trial_num];
end

data.event = trl;  % add new event matrix into data structure (replacing old one)
data.eventsmp = smp;  % add sample times matrix into data structure

% Discard all data after block end marker
endind = find(cellfun(@(x) strcmp(x,'2'),{event.type}));  % '2' marks the end of each block
endind = event(endind).sample;

data.Xgaze          = data.Xgaze(1:endind);
data.Ygaze          = data.Ygaze(1:endind);
data.pupil          = data.pupil(1:endind);
data.times          = data.times(1:endind);
data.sampleinfo     = [1 length(data.times)];
if ~isempty(data.blinksmp), data.blinksmp = data.blinksmp(data.blinksmp(:,2)<=length(data.times),:); end
if ~isempty(data.saccsmp), data.saccsmp   = data.saccsmp(data.saccsmp(:,2)<=length(data.times),:); end



