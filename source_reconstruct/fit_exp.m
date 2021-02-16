function ssr = fit_exp(params,data)


% Function for fitting one-parameter exponential to data. Only appropriate
% for input data with amplitude of 1 @ times=0 and decaying thereafter
% See also get_exp_fit.m
%
% Input:
% params(1) = timescale (tau)
%
% data(:,1) = times (should be a vector, in s, whose first element is zero)
% data(:,2) = values to be fit;

tau = params(1);  % timescale

times = data(:,1);  % transforming times vector before passing through exponential
if times(1)~=0, times = times-times(1); end  % aligning times vector to zero if not already the case

pred = exp(-times./tau);  % calculating exponential

ssr = sum((data(:,2)-pred).^2);  % sum of squared residuals




% % Function for fitting two-parameter exponential to data. Only appropriate
% % for input data with amplitude of 1 @ times=0 and decaying thereafter
% % See also get_exp_fit.m
% %
% % Input:
% % params(1) = offset constant (B)
% % params(2) = timescale (tau)
% %
% % data(:,1) = times (should be a vector, in s, whose first element is zero)
% % data(:,2) = values to be fit;
% 
% B = params(1);  % offset
% tau = params(2);  % timescale
% 
% times = data(:,1);  % transforming times vector before passing through exponential
% if times(1)~=0, times = times-times(1); end  % aligning times vector to zero if not already the case
% 
% pred = (exp(-times./tau)+B)./(1+B);  % calculating exponential
% 
% ssr = sum((data(:,2)-pred).^2);  % sum of squared residuals




% % Function for fitting three-parameter Murray exponential to data.
% % See also get_exp_fit.m
% %
% % Input:
% % params(1) = amplitude constant (A)
% % params(2) = offset constant (B)
% % params(3) = timescale (tau)
% %
% % data(:,1) = times (should be a vector, in s, whose first element is zero)
% % data(:,2) = values to be fit;
% 
% A = params(1);  % amplitude
% B = params(2);  % offset
% tau = params(3);  % timescale
% 
% times = data(:,1);  % transforming times vector before passing through exponential
% if times(1)~=0, times = times-times(1); end  % aligning times vector to zero if not already the case
% 
% pred = A.*(exp(-times./tau)+B);  % calculating exponential
% 
% ssr = sum((data(:,2)-pred).^2);  % sum of squared residuals




