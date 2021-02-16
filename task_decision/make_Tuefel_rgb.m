clear

% Teufel CIE coordinates
x = [0.2858, 0.2779, 0.2731, 0.2721, 0.2756, 0.2841, 0.2977, 0.3149, 0.332, 0.3445, 0.3489, 0.3449, 0.335, 0.3221, 0.3087, 0.2963];  % CIE x coordinates of Teufel colours
y = [0.2587, 0.259, 0.2667, 0.2818, 0.304, 0.3321, 0.3629, 0.3907, 0.4075, 0.4074, 0.3904, 0.3622, 0.331, 0.3027, 0.2806, 0.2659];  % CIE y coordinates of Teufel colours
Y = [42.675, 42.792, 42.923, 43.048, 43.149, 43.209, 43.221, 43.181, 43.096, 42.979, 42.848, 42.722, 42.622, 42.561, 42.550, 42.590];  % CIE Y coordinates of Teufel colours

% Pulling range of colours to be used for colourbar
x = x(7:14);
y = y(7:14);
Y = Y(7:14);

% Get n linearly spaced X coordinates
xfull = [];
for i = 1:length(x)-1
    xfull = [xfull linspace(x(i),x(i+1),50)];
end

% Flipping segment of color plane, after max x value, before interpolation
flippos = find(x==max(x),1,'last')+1:length(x);
temp_x = x;
temp_x(flippos) = temp_x(flippos)+(max(x)-temp_x(flippos)).*2;

flippos_full = find(xfull==max(x),1,'last')+1:length(xfull);
temp_xfull = xfull;
temp_xfull(flippos_full) = temp_xfull(flippos_full)+(max(x)-temp_xfull(flippos_full)).*2;

% Interpolate associated y/Y coordinates
yfull = interp1(temp_x,y,temp_xfull,'pchip');
Yfull = interp1(temp_x,Y,temp_xfull,'pchip');

% Convert from CIE to RGB
rgb = cie2rgb(xfull,yfull,Yfull);

% Save
save('D:\Experiments\Surprise_accumulation\Task\Teufel_rgb.mat','rgb')