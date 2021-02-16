function run_trials(i)

cd('/home/pmurphy/Surprise_acc/wang2002/c_wplus')
rootdir = '/home/pmurphy/Surprise_acc/wang2002/c_wplus/';

ntrials = 50;  % number of trials to simulate

decvar.wplus       = 1.68;  %recurrent connectivity strength
decvar.c           = 19;    %LLR multiplicative scaling

decvar.meanR       = 40;    %input rate for LLR = 0
decvar.dt          = .2;    %time step
decvar.ms          = 5600;  %decision length
decvar.N           = 2000;  %Ncells
decvar.f           = 0.15;  %fraction selective
decvar.vext        = 2.4;   %external noise input
decvar.stim_start  = 400;   %stim on time (ms)      
decvar.samp_dur    = 400;   %sample duration (ms)     

seedval = 1;

if ~exist([rootdir,'output_c',num2str(decvar.c),'_wplus',num2str(decvar.wplus)],'dir')
    mkdir([rootdir,'output_c',num2str(decvar.c),'_wplus',num2str(decvar.wplus)])
end
savepath = [rootdir,'output_c',num2str(decvar.c),'_wplus',num2str(decvar.wplus),filesep];

load(['w',num2str(decvar.wplus),'.mat']);     % load predefined connectivity matrix
load('norm_vars.mat','allLLR','allCPP_M','allUncertainty_M') ;  % load LLR/pCP/uncertainty sequences

i = str2double(i);  % converting string input to number
LLRs = allLLR(ntrials*(i-1)+1:min([ntrials*(i) size(allLLR,1)]),:);  % pull LLR sequences for current iteration
surp = allCPP_M(ntrials*(i-1)+1:min([ntrials*(i) size(allLLR,1)]),:);
uncert = allUncertainty_M(ntrials*(i-1)+1:min([ntrials*(i) size(allLLR,1)]),:);

clear allLLR allCPP_M allUncertainty_M

for t = 1:size(LLRs,1)
    fprintf('\n\nSimulating trial %d of %d... ',t,ntrials) ;
    decvar.LLR = LLRs(t,:);  % define LLR sequence for this trial
    [choice(t,1),sw_relFR(t,:),events] = wang2002(decvar,seedval,w);
end
clear w t

save([savepath,'sim',num2str(i),'.mat'])