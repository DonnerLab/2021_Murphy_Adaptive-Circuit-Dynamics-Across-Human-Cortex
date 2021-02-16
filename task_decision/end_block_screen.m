function end_block_screen(block,loadstr,window,windowRect)

% Calculate current block and total earnings
acc = 0; ntrials = 0;
for b = 1:block;
    load([loadstr,num2str(b),'.mat']);
    acc = acc+nansum(Behav(:,3));
    ntrials = ntrials+size(Behav,1);
end

% Get window center
[~,Y] = RectCenter(windowRect);

% Specify text
str1 = 'Block complete!';
str2 = sprintf('Accuracy for that block = %2.1f percent',nansum(Behav(:,3))/size(Behav,1)*100);
str3 = sprintf('Total accuracy for this session = %2.1f percent',acc/ntrials*100);
str4 = sprintf('You broke fixation during a trial on %d trials.',length(find(Behav(:,6)>0)));
str5 = 'Remember to try as hard as you can to only do this during the rest period between trials!';
str6 = 'Please wait while we save your data...';

% Display text
DrawFormattedText(window,str1,'center',Y-140,[255 255 255]);
DrawFormattedText(window,str2,'center',Y-60,[255 255 255]);
DrawFormattedText(window,str3,'center',Y-20,[255 255 255]);
DrawFormattedText(window,str4,'center',Y+60,[255 255 255]);
DrawFormattedText(window,str5,'center',Y+10,[255 255 255]);
DrawFormattedText(window,str6,'center',Y+180,[255 255 255]);

Screen('Flip', window);

WaitSecs(5);