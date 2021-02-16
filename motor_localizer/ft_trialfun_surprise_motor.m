function [trl, event] = ft_trialfun_surprise_motor(cfg)
% dir_behav = '/mnt/homes/home024/pmurphy/Surprise_accumulation/DataSensLoc/';

%triggers of interest
cue_tr = 36;             %onset of response direction cue
go_tr = 46;              %offset of dot
respL_tr = 56;
respR_tr = 57;
badresp_tr = 58;
startblock_tr = 5;       %start block
endblock_tr = 6;         %end block

% read the header, contains the sampling frequency
hdr = ft_read_header(cfg.dataset);

subject = cfg.subj;
session = cfg.sess;
recording = cfg.rec;

% % Load dot angles
% dir_behav = [dir_behav,subject,'/S',session,'/Behaviour/'];
% cfile = dir([dir_behav,'*.mat']);
% load([dir_behav, cfile.name]);
% angles = Behav(:,1);

% read MEG events
if ~(strcmp(subject,'QNV') && strcmp(session,'3'))
    if isfield(cfg, 'event')
       event = cfg.event;
    else
       event = ft_read_event(cfg.dataset);
    end
else
    % Load reconstructed triggers for session 3 subject QNV
    % Build and save artificial event structure that can by used by trial function
    load(['/mnt/homes/home024/chernandezsav/meg_data/surprise/preprocessed/Data/qnv_corr_coef/' 'qnv_triggersET_3_' recording '.mat']);
    event = [];
    for i = 1:length(triggersET)
        event(i,1).type = 'UPPT001';
        event(i,1).sample = triggersET(i,2);
        event(i,1).value = triggersET(i,1);
        event(i,1).offset = [];
        event(i,1).duration = [];
    end
end

% Extract real RTs from UPPT002 (some subjects don't do task correctly)
events = event;  % load events

events = events(ismember({events.type},{'UPPT001','UPPT002'})); % isolate only stim/response triggers
events = events([events.sample]>=events([events.value]==5).sample & [events.sample]<=events([events.value]==6).sample);  % pull motor loc segment

if ~isempty(find([events.value]==36))
    RTs = [];
    for stim = find([events.value]==36)  % loop through each trial
        i = stim;
        while ~strcmp(events(i).type,'UPPT002')
            i = i+1;
            if strcmp(events(i).type,'UPPT002')
                RTs(end+1,1) = (events(i).sample - events(stim).sample) / 1200;
                break
            elseif i == length(events) | events(i).type==36  % if there's not response trigger before end of trial or block
                RTs(end+1,1) = nan;
                break
            end
        end
    end
end

% pull only motor localizer events
event = event(strcmp({event.type},'UPPT001'));
event = event([event.value]==cue_tr | [event.value]==go_tr | [event.value]==respL_tr | [event.value]==respR_tr | [event.value]==badresp_tr | [event.value]==endblock_tr);

assert(isequal(length(find([event.value]==cue_tr)),length(RTs)),'ERROR: Mismatch between onsets and RTs')

% Loop through events and add epochs to trl matrix
trl = [];
if ~isfield(cfg.trialdef, 'eventvalue')
    ctype = cue_tr;  % default epoching around stim onset if not explicitly specified
else ctype = cfg.trialdef.eventvalue;
end

ticker=0;
for i = 1:length(event)
    if event(i).value == ctype
        
        ticker = ticker+1;
        if ~isfield(cfg.trialdef, 'prestim')
            trloff = event(i).offset;
            trlbeg = event(i).sample;  % default prestimulus period of 0 seconds
        else
            trloff = round(-cfg.trialdef.prestim*hdr.Fs);
            trlbeg = event(i).sample + round(-cfg.trialdef.prestim*hdr.Fs);
        end
        
        if ~isfield(cfg.trialdef, 'poststim')
            trlend = event(i).sample + 1*hdr.Fs;  % default trial length of 1 second
        else
            trlend = event(i).sample + cfg.trialdef.poststim*hdr.Fs;
        end
        
        if i<length(event)-1
            if event(i+2).value==respL_tr || event(i+2).value==respR_tr
                resp = event(i+2).value;
            else resp = nan;
            end
        end
        
        new_trial = [trlbeg,trlend,trloff,resp,RTs(ticker)];
        if isempty(trl)
          trl = new_trial;
        else
          trl = [trl; new_trial];
        end
        
    end
end
cfg.trl = trl;
cfg.full_event = event;

end

