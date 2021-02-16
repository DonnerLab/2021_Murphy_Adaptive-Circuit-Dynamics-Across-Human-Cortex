% Script for simulating behaviour generated by Glaze et al. normative model
% during two-alternative sequential evidence accumulation task where
% samples are drawn from one of 2 distributions with same std but different
% means, and mean can switch at fixed hazard rate.

function CP = FitPsi_sim_fast_npLLR_InconUp(LLR,nsamps,L_nm1,rLLR,psi_n,gains,noise,incon)

% Apply full non-parametric gain function to LLRs
gains = cumsum(gains);
gains = [fliplr(gains) 0 gains];  % creating full vector of gain terms, symmetric around LLR=0
rLLR = [sort(-rLLR) 0 rLLR];      % creating full vector of matching original LLRs

for s = 1:size(LLR,2)
    LLR(:,s) = interp1(rLLR,rLLR.*gains,LLR(:,s),'spline');
end

% Create full non-parametric Glaze transfer function
psi_n = [fliplr(-psi_n) 0 psi_n];  % creating full vector of psi terms, symmetric around L_nm1(0)=psi_n=0
L_nm1 = [sort(-L_nm1) 0 L_nm1];  % creating full vector of matching L_nm1s

% Increment belief for each sample, scaling increments by interpolated phi function
LPRout = zeros(size(LLR,1),1);
for s = 1:size(LLR,2)
    LPRout(:,end+1) = nan(size(LPRout,1),1);
    
    LLR(abs(sign(LPRout(:,end-1))-sign(LLR(:,s)))==2,s) = LLR(abs(sign(LPRout(:,end-1))-sign(LLR(:,s)))==2,s).*incon;  % apply gain factor to new samples that are inconsistent with current belief
    
    LPRout(LPRout(:,end-1)>=min(L_nm1) & LPRout(:,end-1)<=max(L_nm1),end) = ...
        LLR(LPRout(:,end-1)>=min(L_nm1) & LPRout(:,end-1)<=max(L_nm1),s) + interp1(L_nm1,psi_n,LPRout(LPRout(:,end-1)>=min(L_nm1) & LPRout(:,end-1)<=max(L_nm1),end-1),'spline');
    LPRout(LPRout(:,end-1)<min(L_nm1),end) = LLR(LPRout(:,end-1)<min(L_nm1),s) + psi_n(1);   % setting any posterior beliefs below lowest fitted phi to current parameter estimate for lowest fitted phi
    LPRout(LPRout(:,end-1)>max(L_nm1),end) = LLR(LPRout(:,end-1)>max(L_nm1),s) + psi_n(end);   % setting any posterior beliefs above highest fitted phi to current parameter estimate for highest fitted phi
end

% Retrieve final LPRs for each trial (accounting for variable sequence lengths)
subinds = sub2ind(size(LPRout),1:size(LPRout,1),nsamps'+1);
LPRfinal = LPRout(subinds)';

% Calculate choice probabilities based on LPR @ end of each sequence, scaled by noise
CP = 0.5+0.5.*erf(LPRfinal./noise);