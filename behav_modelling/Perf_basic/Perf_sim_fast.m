% Creates sequence of updated beliefs given some input sequence of log
% likelihood ratios, according to perfect accumulation.

% B = gain term applied to sample-wise LLRs
% noise = divisive scaling factor applied to final belief

function CP = Perf_sim_fast(LLRin,nsamps,B,noise)

% Apply gain term to LLRs (represents subjective component of generative variance)
LLRin = LLRin.*B;

% Increment belief for each sample
LPRout = zeros(size(LLRin,1),1);
for s = 1:size(LLRin,2)
    LPRout(:,end+1) = LPRout(:,end)+LLRin(:,s);
end

% Retrieve final LPRs for each trial (accounting for variable sequence lengths)
subinds = sub2ind(size(LPRout),1:size(LPRout,1),nsamps'+1);
LPRfinal = LPRout(subinds)';

% Calculate choice probabilities based on LPR @ end of each sequence, scaled by noise
CP = 0.5+0.5.*erf(LPRfinal./noise);
