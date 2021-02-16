% Takes sequences of samples, computes trial-wise choice probabilities, and
% computes cross-entropy error term (e) w.r.t. observed responses.

% pm = [phi_n, LLR_n, noise]

function e = FitPsi_npLLR_InconUp_cross_entropy_fitting(pm)

% Retrieving global variables from initializing script
global LLRin choices nsamps L_nm1 rLLR

% Looping through PSO particles
for p = 1:size(pm,1)
    % Generate choice probabilities given current parameter set
    CPs = FitPsi_sim_fast_npLLR_InconUp(LLRin,nsamps,L_nm1,rLLR,pm(p,1:length(L_nm1)),pm(p,length(L_nm1)+1:length(L_nm1)+length(rLLR)),pm(p,end-1),pm(p,end));
    CPs(CPs==1) = 0.9999; CPs(CPs==0) = 0.0001;  % adjusting extreme CPs that lead to inf cross-entropy values
    
    % Calculate cross-entropy between model choice probabilities and observed responses,
    % adding Tikhonov regularization penalty terms on derivatives of phi & LLR transfer functions to promote smooth function fits (long-time settings: 3/5 & 1/20)
    e(p,1) = -sum((ones(size(choices))-choices).*log(1-CPs)+(choices.*log(CPs))) + (1/2)*sum(diff([0 pm(p,1:length(L_nm1))]).^2) + (1/20)*sum(pm(p,length(L_nm1)+1:length(L_nm1)+length(rLLR)).^2);
    
    % assert(~isnan(e(p)),'ERF=nan for param set H=%1.4f, B=%1.4f, noise=%1.4f...',pm(p,1),pm(p,2),pm(p,3))
end