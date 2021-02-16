% Performs two-tailed, within-subjects cluster-based permutation test on
% one-dimensional time-series data. Clusters are formed via thresholded
% paired-samples T-scores; permutation distribution is formed from the peak
% absolute summed T-score across clusters formed within each permutation.
% Permutations consist of randomly shuffling condition labels independently
% per subject.
%
% Inputs:  data_in = subj*time*cond matrix where cond indexes task condition,
%                    subj indexes subject and time indexes the dimension of
%                    the data being input (e.g. time, space, frequency)
%          nperm   = number of permutations
%          clustalpha = alpha level applied to first-level statistics -
%                    determines how clusters are defined (usually 0.05)
%          alpha   = alpha level applied to cluster test statistic (usually 0.05)
%
% Outputs: sig_ts  = vector of length time, with 1 marking significant
%                    cluster-corrected time points and nan elsewhere
%          sig_ts_uncorr = vector of length time, with 1 marking significant
%                    time points without correction and nan elsewhere
%          cP      = p-values for each cluster
%          permdist = permutation distribution (for possible sanity checks)
%
% Peter Murphy, 16/05/17
% 
% Dropped some for-loops (decreased computation time) 
% KD, November 18

function [sig_ts,sig_ts_uncorr,cP,permdist] = cluster_permWS_fast(data,nperm,clustalpha,alpha)

% Isolate clusters in observed data and calculate associated test statistics
[Hobs, ~, ~, ts] = ttest(data(:,:,1),data(:,:,2),'alpha',clustalpha);
Hobs = squeeze(Hobs);
Tobs = squeeze(ts.tstat);
sig_ts_uncorr = Hobs; sig_ts_uncorr(sig_ts_uncorr==0) = nan;

if ~isempty(find(Hobs==1))  % making sure there's a significant effect at at least one time point in observed data
    [Cobs,nclust] = bwlabel(Hobs);
    for c = 1:nclust
        cTobs(c) = abs(sum(Tobs(Cobs==c)));  % take absolute sum of T-scores within each cluster
    end
else
    sprintf('No significant effects in observed data... :-/')
    Cobs = 1;
    nclust = 1;
    cTobs = max(abs(Tobs));
end

% Run successive permutations
clabels = [1 2];
for p = 1:nperm
    % Shuffle condition labels
    cdata=[];
    corder = randi(clabels,size(data,1),1);
    for s = 1:size(data,1)
        cdata(s,1:size(data,2),1:2) = cat(3,data(s,:,clabels==corder(s)),data(s,:,clabels~=corder(s)));
    end
    
    % Isolate clusters in permuted data and store peak test statistic
    Hperm=[]; Tperm=[]; cTperm=[];
    [Hperm,~,~,ts] = ttest(cdata(:,:,1),cdata(:,:,2),'alpha',clustalpha);
    Hperm = squeeze(Hperm);
    Tperm = squeeze(ts.tstat);
    
    if isempty(find(Hperm==1))  % if there's no significant effect at any time point in this permuation, just take max standalone T-score
        TpermMax(p) = max(abs(Tperm));
    else
        [Cperm,nclustperm] = bwlabel(Hperm);
        %[Cperm,nclustperm] = get_clusts(Hperm);
        for c = 1:nclustperm
            cTperm(c) = abs(sum(Tperm(Cperm==c)));  % take absolute sum of T-scores within each cluster
        end
        TpermMax(p) = max(cTperm);
    end
end

% Calculate corrected p-value for each observed cluster and form vector of significant clusters
sig_ts = sig_ts_uncorr;
for c = 1:nclust
    cP(c) = 1-(length(find(TpermMax<=cTobs(c)))/nperm);
    if cP(c)>alpha
        sig_ts(Cobs==c) = nan;
    end
end

permdist = TpermMax;


end


% Subfunction for isolating temporal clusters from significance flags
function [c,n] = get_clusts(H)

c = nan(1,length(H)); % initializing vector of cluster IDs

n = 1; ts = find(H==1,1,'first'); sigperiod = 1;
while ts<=length(H)
    if sigperiod
        if H(ts)==1
            c(ts) = n;
        else sigperiod = 0;
            n = n+1;
        end
        ts = ts+1;
    else
        if H(ts)==1
            c(ts) = n;
            sigperiod = 1;
        end
        ts = ts+1;
    end
end
end
