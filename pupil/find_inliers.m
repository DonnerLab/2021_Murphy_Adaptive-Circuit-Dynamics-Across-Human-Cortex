% Function for finding elements of input vector/matrix that lie within
% outlier discrimination threshold. If mat_in is matrix, variables are
% assumed to be represented in columns. Cutoff is a z-threshold, and
% function will normalize mat_in appropriately.

function inliers = find_inliers(mat_in,cutoff)

% Z-score input variables
mat_in = nanzscore(mat_in,0,1);

% Iterate through input variables and find inliers
ins=[];
for i = 1:size(mat_in,2)
    ins(:,i) = mat_in(:,i)>(mean(mat_in(:,i))-std(mat_in(:,i))*cutoff) & mat_in(:,i)<(mean(mat_in(:,i))+std(mat_in(:,i))*cutoff);
end

% Find rows satisfing all thresholds
inliers = find(sum(ins,2)==size(mat_in,2));

