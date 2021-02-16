function [imgout,alphaout] = make_gaussian_PDF_cbar_tex(mean,sigma,lims,ppd,res,r_in,r_ext,yscaling)

% mean     =  mean angle (in degrees) of evidence distribution (Gaussian) - 0 degrees = vertical midline
% sigma    =  std of evidence distribution
% ppd      =  pixels per d.v.a
% res      =  resolution of display [horizontal, vertical]
% r_in     =  inner radius (d.v.a.) of elliptical LLR colour bar along x-axis
% r_ext    =  outer radius (when added to r_in) of elliptical LLR colour bar along x-axis
% yscaling =  multiplicative scaling factor by which to decrease eccentricity along yhe y-axis (to mimic x/y asymmetry in human vision)

mean = mean+90;   % rereference by 90 degrees to translate into polar coordinates
lims = lims+90;

dva_res = res./2./ppd;  % screen resolution expressed as d.v.a.

[x, y] = meshgrid(linspace(-dva_res(1), dva_res(1), res(1)), linspace(-dva_res(2), dva_res(2), res(2)));  % create matrix of pixels
t = rad2deg(atan2(y, x));        % calculate polar angle coordinate for each pixel
r = sqrt(y.^2 + x.^2);  % calculate polar radius coordinate for each pixel
imgout = normpdf(t, mean, sigma);    % draw distribution likelihood values for each pixel from Gaussian pdf

y_rescaled = y.*(1/(1-yscaling));    % rescaling y coords to create apropriate-looking scaled colourbar (not sure why this works, but looks better than rescaling method below)
tm = rad2deg(atan2(y_rescaled, x));
rm = sqrt(y_rescaled.^2 + x.^2);
mask = rm<r_in | rm>(r_in+r_ext) | (tm<lims(1) | tm>lims(2));
imgout(mask==1) = nan;   % only keep colour bar in pixels of interest; set everything else to nan
imgout = (imgout-nanmin(nanmin(imgout)))./(nanmax(nanmax(imgout))-nanmin(nanmin(imgout)));  % rescaling colorbar such that range = [0,1]

% imgout(r<r_in | r>(r_in+r_ext) | (t<lims(1) | t>lims(2))) = nan;   % only keep colour bar in pixels of interest; set everything else to nan
% imgout = (imgout-nanmin(nanmin(imgout)))./(nanmax(nanmax(imgout))-nanmin(nanmin(imgout)));  % rescaling colorbar such that range = [0,1]
% imgout = [nan(round(res(2)*yscaling/2),res(1)); imresize(imgout,[round(res(2)*(1-yscaling)) res(1)])];  % rescaling colorbar along y-axis
% imgout(end+1:res(2),:) = nan(res(2)-size(imgout,1),res(1));

alphaout = ones(size(imgout)).*255;   % intializing matrix of transparency values
alphaout(isnan(imgout)) = 0;  % setting everything but colorbar to be completely transparent
