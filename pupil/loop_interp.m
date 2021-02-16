% Interpolate pupil time-series

clear; close all; clc

user_input = true;

loadpath = 'D:\Experiments\Surprise_accumulation\Analysis\Pupil\2.matConverted\';
savepath = 'D:\Experiments\Surprise_accumulation\Analysis\Pupil\3.Interpolated\';
addpath(genpath('D:\Experiments\Surprise_accumulation\Analysis\Pupil'));

allsubj = {'DHB','TFD','EXF','JTB','TNB','QNV','PDP','GSB','OMF','NIF','ECB','TSJ','KSV','HBC','EMB','DCB','EXG'};

% Loop through subjects
for_checking = {};
% Loop through individual files
for subj = 1:length(allsubj);
    
    bnames = dir([loadpath,allsubj{subj},'*.mat']);
    
    for b = 1:length(bnames)
        fprintf('Subj %s, %s\n',allsubj{subj},bnames(b).name)  % print progress
        
        % Load current asc conversion
        inFile = [loadpath,bnames(b).name];
        outFile = [savepath,bnames(b).name(1:end-4),'_interp.mat'];
        
        load(inFile)
        
        % Run interpolation routine & add results to data structure
        [newpupil, newXgaze, newYgaze, newblinksmp, badsmp] = blink_interpolatePM(data,bnames(b).name);
        data.pupil = newpupil;
        data.Xgaze = newXgaze;
        data.Ygaze = newYgaze;
        data.newblinksmp = newblinksmp;
        data.badsmp = badsmp;
        
        pause(2)  % wait for a second to give user a better look when user input is specified to false
        
        if user_input
            % Get user input on interpolation quality if specified
            instr = input('Happy with interpolation? (y/n) ...','s');
            
            if strcmp(instr,'y')
                save(outFile,'data')
            else
                for_checking{end+1} = bnames(b).name(1:end-4);
            end
        else % Otherwise just save
            save(outFile,'data')
        end
        close all
    end
end

% Reminder of files to be checked further
if ~isempty(for_checking)
    fprintf('\nThe following files need to be inspected more closely:\n')
    for f = 1:length(for_checking);
        fprintf('%s\n',for_checking{f})
    end
end