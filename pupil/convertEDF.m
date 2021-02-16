% Convert from Eyelink .edf files to .asc files using edf2asc converter

loadpath = 'D:\Experiments\Surprise_accumulation\Data\';
savepath = 'D:\Experiments\Surprise_accumulation\Analysis\Pupil\1.ascConverted\';

addpath D:\Experiments\Surprise_accumulation\Analysis\Pupil
addpath 'C:\Program Files\MATLAB\fieldtrip-20160221'  % tell Matlab where FieldTrip is
ft_defaults

subj = 'DHB';
sess = [1 2 3];
nblocks = [8 8 9];

% Loop through individual files
for s = 1:length(sess);
    for b = 1:nblocks(s);
        
        fprintf('Subj %s, session % d, block %d\n',subj,s,b)
        
        % Define file names
        edfFile = [loadpath,subj,'\S',num2str(sess(s)),'\Eyetracking\',subj,'_',num2str(sess(s)),'_',num2str(b),'.edf'];
        ascFile = [loadpath,subj,'\S',num2str(sess(s)),'\Eyetracking\',subj,'_',num2str(sess(s)),'_',num2str(b),'.asc'];
%             edfFile = [loadpath,subj,'\Training\Eyetracking\',subj,'_',num2str(s),'_',num2str(b),'.edf'];
%             ascFile = [loadpath,subj,'\Training\Eyetracking\',subj,'_',num2str(s),'_',num2str(b),'.asc'];
        matFile = [savepath,subj,'_',num2str(sess(s)),'_',num2str(b),'.mat'];
        
            % edf2asc
            if ~exist(ascFile)
                system(sprintf('%s %s -failsafe -input', [loadpath,'edf2asc-linux'], edfFile));
                assert(exist(ascFile, 'file') > 1, 'Edf not properly converted...');  % check that asc has actually been created
            end
        
        % Read the asc file into matlab & save
        %if ~exist(matFile)
            asc = read_eyelink_ascNK_AU(ascFile);
            save(matFile,'asc')
        %end
    end
end

