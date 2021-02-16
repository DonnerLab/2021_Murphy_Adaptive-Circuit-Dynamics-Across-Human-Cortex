function [x, y, p] = eye_voltage2gaze(raw, ranges, screen_x, screen_y, ch_mapping)  
    %Converts analog output of EyeLink 1000+ to gaze coordinates
    %Based in Niklas python script
    minvoltage = ranges(1);
    maxvoltage = ranges(2);
    minrange=0;
    maxrange = 1;
    screenright = screen_x(1);
    screenleft = screen_x(2);
    screenbottom = screen_y(1);
    screentop = screen_y(2);
    
    %obtain the idx of the channels of interest
    %TODO: obtain idx based on names, for now, i got the fixed idx 380,381,382    
    idx = 380;%UADC002
    idy = 382;%UADC004
    idp = 381;%UADC003
        
    R = (raw(idx,:)-minvoltage)/(maxvoltage-minvoltage);
    S = R*(maxrange-minrange)+minrange;
    x = S*(screenright-screenleft+1)+screenleft;
    
    R = (raw(idy,:)-minvoltage)/(maxvoltage-minvoltage);
    S = R*(maxrange-minrange)+minrange;
    y = S*(screenbottom-screentop+1)+screentop;

    p = raw(idp,:);  
end