% Computes time of first button press on each trial of motor localizer task
% based on *response triggers in the MEG data* - a viridical readout of
% the timing of button presses (as opposed to task behavioural data, which
% only records motor responses during designated response period). This is
% used to check whether some participants didn't wait until go cue to make
% their responses, thus rendering motor localizer data unuseable.

function [] = motor_loc_check_resp_chan(subject, session, recording)
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
% %addpath '/Users/caro/Documents/UKE/Matlab/Toolboxes/fieldtrip-20170322'
ft_defaults
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts/PreprocScripts/
cd /mnt/homes/home024/pmurphy/meg_data/surprise/  % specify path where dataset is located
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/PreprocessedMotor/Data/Resp_channel/';

subjfile = dir([cd,'/',subject,'-',session,'_Surprise_','*',recording,'.ds']);

% ==================================================================
% 1. READ EVENTS FROM ENTIRE RECORDING & DEFINE TRIALS
% ==================================================================
events = ft_read_event(subjfile(1).name);  % load events

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
    
    % save
    namefile = strcat(subject,'-',session,'_trueRTs.mat');
    save([savepath,namefile],'RTs','events','-v7.3');
end
