function err = fit_pstim_from_LLR(log_pstim,xy)

stim = xy(:,1);
LLR = xy(:,2);

% Interpolate pdf to full resolution and normalize so integral = 1
fullpstim = exp(interp1(stim',log_pstim,[min(stim) -90:1:90 max(stim)],'spline'));
fullpstim = log(fullpstim./sum(fullpstim));  % will still be in log space afer normalization

% Calculate predicted LLR for this fit iteration
LLRpred = fullpstim-fliplr(fullpstim);

% Interpolate input LLR function to matched resolution
LLRin = interp1(stim,LLR,[min(stim) -90:1:90 max(stim)],'spline');

% Calculate SSR between input and predicted LLR functions
err = sum((LLRin-LLRpred).^2);
