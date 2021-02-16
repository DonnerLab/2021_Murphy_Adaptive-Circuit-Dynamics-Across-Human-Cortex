function sacc = check_saccade(eye, xc, yc, ppd)
pause(0.002)
sacc = 0;
[samples, ~, ~] = Eyelink('GetQueuedData');
%sprintf('samples')
%samples
if eye==0
    x = (samples(14,:)-xc)/ppd;
    y = (samples(16,:)-yc)/ppd;
else
    x = (samples(15,:)-xc)/ppd;
    y = (samples(17,:)-yc)/ppd;
end
   
d = (x.^2 + y.^2).^.5;
a=d(2:length(d));
if any(a>4)
    sacc = 1;
end
