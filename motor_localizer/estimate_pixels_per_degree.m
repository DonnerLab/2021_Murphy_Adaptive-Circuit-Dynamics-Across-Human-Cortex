function ppd = estimate_pixels_per_degree()
    dist=60;
    width=50;
    height=32;
    resolution=[1920,1080];
    o = tan(0.5*pi/180)*dist;
    ppd = 2*o*resolution(1)/width;
end

