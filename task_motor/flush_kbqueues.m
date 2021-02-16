function flush_kbqueues(devices)
for device = devices
    PsychHID('KbQueueFlush', device);
end
