% Script for calculating existing number of blocks within a specified
% subject/session directory

function b = count_blocks(dir_in)

files=dir(dir_in);   % get list of all files in directory
filenames = {files.name};   % pull only file names

bs=[];
for f = 1:length(filenames);
    
    current_name = filenames{f};  % current filename
    
    dotmark=[];
    l = length(current_name);
    while ~strcmp(current_name(l),'_')
        l = l-1;
        if strcmp(current_name(l),'.')
            dotmark = l;  % marking position of dot
        end
    end
    undermark = l;  % marking position of last underscore
    
    bs(end+1) = str2num(current_name(undermark+1:dotmark-1));  % pulling block number from current filename
end

b = num2str(max(bs)+1);  % calculate number of new block to be presented