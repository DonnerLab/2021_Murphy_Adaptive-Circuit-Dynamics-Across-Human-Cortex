function [checkerid,checkerrect] = createCircularChecker(win, diam, freq, ppd, transparency)

% inputs:
%     win   =  PTB window
%     diam  =  diameter of checkerboard circle (in d.v.a.)
%     freq  =  spatial frequency of checkerboard (in cycles per d.v.a.)
%     ppd   =  pixels per d.v.a.


% get number of pixels per cycle - forcing this to be integer, so spatial freq might be ever-so-slightly off as a result)
ppc = round(ppd/freq);

% if non-integer number of cycles fit within desired diameter, make initial
% checkerboard slightly larger (and trim later)
if mod(diam,1/freq)==0
    c_diam = ppc*freq*diam;
else
    c_diam = ppc*freq*diam + ppc*(1-(mod(diam,1/freq)*freq));
end

if mod(c_diam,2)~=0, uneven_number_of_pixels, end  % making sure there's an even number of pixels to pass on to next steps

% specifying boundary of texture to be drawn
checkerrect = [0 0 c_diam c_diam];

% compute number of complete cycles to be drawn 
numCheckers =  ceil(diam*freq);

% make an atomic checkerboard (one complete cycle)
miniboard = eye(2,'uint8') .* 255;

% repeat it by number of desired cycles
cb1 = repmat(miniboard, numCheckers, numCheckers)';
%cb2 = repmat(ones(2,2).*255-miniboard, numCheckers, numCheckers)';
cb2 = 255-cb1;

% scale the images up
cb1_full = imresize(cb1,[c_diam c_diam],'box');
cb2_full = imresize(cb2,[c_diam c_diam],'box');

% turn all but aperture of desired diam transparent
dva_res = c_diam/2/ppd;  % screen resolution expressed as d.v.a.
[x, y] = meshgrid(linspace(-dva_res, dva_res, ceil(c_diam)), linspace(-dva_res, dva_res, ceil(c_diam)));  % create matrix of pixels
r = sqrt(y.^2 + x.^2);  % calculate polar radius coordinate for each pixel
cb_alpha = ones(size(r)).*255;   % intializing matrix of transparency values
cb_alpha(r>diam/2) = 0;  % setting all but central circle to be transparent

% make final checkerboard texture
cb1_full(:,:,2) = cb_alpha.*(1-transparency);
cb2_full(:,:,2) = cb_alpha.*(1-transparency);
checkerid(1) = Screen('MakeTexture', win, cb1_full);
checkerid(2) = Screen('MakeTexture', win, cb2_full);

% 
% % convert diameter from d.v.a. to pixels
% diam = diam*ppd;
% 
% % if non-integer number of cycles fit within desired diameter, make initial
% % checkerboard slightly larger (and trim later)
% if mod(diam,ppd/freq)==0
%     c_diam = diam;
% else
%     c_diam = diam + (ppd/freq - mod(diam,ppd/freq));
% end
% checkerrect = [0 0 ceil(c_diam) ceil(c_diam)];  % specify boundary of texture to be drawn
% 
% % compute number of complete cycles to be drawn 
% numCheckers =  ceil(diam/(ppd/freq));
% 
% % make an atomic checkerboard (one complete cycle)
% miniboard = eye(2,'uint8') .* 255;
% 
% % repeat it by number of desired cycles
% cb = repmat(miniboard, numCheckers, numCheckers)';
% 
% % scale the images up
% cb_full = imresize(cb,[ceil(c_diam) ceil(c_diam)],'box');
% 
% % turn all but aperture of desired diam transparent
% dva_res = ceil(c_diam)./2./ppd;  % screen resolution expressed as d.v.a.
% [x, y] = meshgrid(linspace(-dva_res, dva_res, ceil(c_diam)), linspace(-dva_res, dva_res, ceil(c_diam)));  % create matrix of pixels
% r = sqrt(y.^2 + x.^2);  % calculate polar radius coordinate for each pixel
% cb_alpha = ones(size(r)).*255;   % intializing matrix of transparency values
% cb_alpha(r>diam/2/ppd) = 0;  % setting all but central circle to be transparent
% 
% % make final checkerboard texture
% cb_full(:,:,2) = cb_alpha;
% checkerid = Screen('MakeTexture', win, cb_full);


