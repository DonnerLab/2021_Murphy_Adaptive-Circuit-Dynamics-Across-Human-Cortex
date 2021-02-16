function sacc = check_saccade(xgaze, ygaze, xcenter, ycenter, ppd)
    sacc = 0;
    x = (xgaze-xcenter)/ppd;
    y = (ygaze-ycenter)/ppd;
    d = (x.^2 + y.^2).^.5;
    a=d(2:length(d));
    if any(a>4)
            sacc = 1;
end

