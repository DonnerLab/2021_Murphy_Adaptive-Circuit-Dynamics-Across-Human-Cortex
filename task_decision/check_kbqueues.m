function [keyIsDown firstPress] = check_kbqueues(devices)
firstPress = boolean(zeros(1, 256)); 
keyIsDown = false;
for device = devices
    [kD, fP] = PsychHID('KbQueueCheck', device);
    keyIsDown = keyIsDown | kD;
    firstPress = firstPress | boolean(fP);
end

end
