function [p_out,p_nans] = interp_nans(p_in,badsmps)

% make the pupil NaN within all identified blink windows
for b = 1:size(badsmps,1),
    p_in(badsmps(b,1):badsmps(b,2)) = NaN;
end

% store indices of nans
p_nans = find(isnan(p_in));

% interpolate linearly
p_in(isnan(p_in)) = interp1(find(~isnan(p_in)), ...
    p_in(~isnan(p_in)), find(isnan(p_in)), 'linear');

% also extrapolate ends
p_in(isnan(p_in)) = interp1(find(~isnan(p_in)), ...
    p_in(~isnan(p_in)), find(isnan(p_in)), 'nearest', 'extrap');

% make output
p_out = p_in;