% Takes sequences of samples, computes trial-wise choice probabilities
% using Glaze model with specified hazard rate (H), B-slope (exponent
% applied to LLR transfer function), B-scale (gain factor applied to
% normalized non-linear LLR transfer function) and noise, and computes
% cross-entropy error term (e) with respect to observed responses.

% pm(1:4) = [H, Bslope, Bscale, noise]

function e = Glaze_cross_entropy_fitting_npLLR_InconUp(pm)

% Retrieving global variables from initializing script
global LLRin choices nsamps rLLR

% Looping through PSO particles
for p = 1:size(pm,1)
    % Generate choice probabilities given current parameter set
    CPs = Glaze_sim_fast_npLLR_InconUp(LLRin,nsamps,rLLR,pm(p,1),pm(p,2),pm(p,3),pm(p,4:end));
    CPs(CPs==1) = 0.9999; CPs(CPs==0) = 0.0001;  % adjusting extreme CPs that lead to inf cross-entropy values
    
    % Calculate cross-entropy between model choice probabilities and observed responses,
    % add Tikhonov regularization penalty term on derivative of additive LLR gains (pm(*5*:end)) to promote smooth function fits
    e(p,1) = -sum((ones(size(choices))-choices).*log(1-CPs)+(choices.*log(CPs))) + 1/20*sum(diff(pm(p,5:end)).^2);
    
    % assert(~isnan(e(p)),'ERF=nan for param set H=%1.4f, B=%1.4f, noise=%1.4f...',pm(p,1),pm(p,2),pm(p,3))
end