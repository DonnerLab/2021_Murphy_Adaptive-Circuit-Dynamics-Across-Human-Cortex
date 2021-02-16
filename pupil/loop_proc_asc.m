loadpath = 'D:\Experiments\Surprise_accumulation\Data\';
savepath = 'D:\Experiments\Surprise_accumulation\Analysis\Pupil\2.matConverted\';

addpath D:\Experiments\Surprise_accumulation\Analysis\Pupil
addpath 'C:\Program Files\MATLAB\fieldtrip-20160221'  % tell Matlab where FieldTrip is
ft_defaults

allsubj = {'DHB','TFD','EXF','JTB','TNB','QNV','PDP','GSB','OMF','NIF','ECB','TSJ','KSV','HBC','EMB','DCB','EXG'};

% Loop through subjects
for subj = 1:length(allsubj)
    
    sdirs = dir([loadpath,allsubj{subj}]); sdirs = sdirs(3:end);
    
    for s = 1:length(sdirs)
                
        sesspath = [loadpath,allsubj{subj},filesep,sdirs(s).name,filesep];
        bnames = dir([sesspath,'Eyetracking',filesep,'*.asc']);
        
        for b = 1:length(bnames)
            
            if strcmp(sdirs(s).name(1),'S')   % making sure this isn't a training block
                fprintf('Subj %s, %s\n',allsubj{subj},bnames(b).name)
                
                % Define file names
                ascFile = [sesspath,'Eyetracking',filesep,bnames(b).name];
                matFile = [savepath,bnames(b).name(1:end-4),'.mat'];
                
                if ~exist(matFile)
                    % Read the asc file into matlab
                    asc = read_eyelink_ascNK_AU(ascFile);
                    
                    % Create data and events structures, and matrices of blink/saccade times
                    data = asc2dat(asc); clear asc
                    
                    % Create refined event structure
                    data = Surprise_trialfun(data,b);
                                        
                    % Save output
                    save(matFile,'data')
                end
            end
        end
    end
end