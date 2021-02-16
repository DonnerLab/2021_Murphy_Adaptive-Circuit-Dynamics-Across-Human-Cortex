% Uses Glaze model to generate finite number of observed choices for a
% given combination of H, noise and inconsistency upweighting parameters.
% In this case the stimulus->LLR transfer function is estimated
% non-parametrically, dividing the stimulus space into sections, estimating
% a weight for each section, and interpolating the full function from
% these estimated weights.

function [pm_fit,err] = Glaze_fit_npLLR_InconUp(subj)

% Set global variables to pass to PSO algorithm
global LLRin choices nsamps rLLR

% Define discrete LLRs for which corresponding gain factor should be estimated
rLLR = [0.5:0.5:3.0 3.35 3.7];  % assumes dot->LLR transfer function is symmetric around 0; max |LLR| given task is ~3.64

% Path stuff
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Modelling/Glaze_npLLR_InconUp
addpath(genpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Simulations/particle_swarm_Glaze'))
loadpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Data/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Modelling/Glaze_npLLR_InconUp/Fits/';

% Define parameter estimation settings
range_p.H = [0.0001 0.5];  % parameter bounds
range_p.noise = [0.0001 10];
range_p.incon = [0.01 10];
range_p.LLRn = repmat([0 10],length(rLLR),1);  % parameter bounds - each gain factor represents an *increment* to the last - lower bound of zero constrains function to be monotonically increasing

mv = [0.005;...    % maximum particle velocities
    0.2;...
    0.15;...
    ones(length(rLLR),1).*0.2;]';

seeds.H = [0.08 0.03];  % good seed distributions for a selection of particles - [mean sd]
seeds.noise = [1 0.5];
seeds.incon = [1 0.5];
seeds.LLRn = [[1 0.5]; repmat([0.4 2],length(rLLR)-1,1)];  % these seeds will slightly bias towards linear LLR transfer, but some instances will have more extreme gain factors

% Seed random number generator
rng('default')
rng('shuffle')

% Load data for current subject
allsubj = {'DHB','TFD','EXF','JTB','TNB','QNV','PDP','GSB','OMF','NIF','ECB','TSJ','KSV','HBC','EMB','DCB','EXG'};

subji = find(strcmp(allsubj,subj));
fprintf('Pulling behavioural data for subject %d of %d: %s\n',subji,length(allsubj),subj)
LLRin=[]; choices=[];

sdirs = dir([loadpath,subj]); sdirs = sdirs(3:end);
for s = 1:length(sdirs)
    
    sesspath = [loadpath,subj,filesep,sdirs(s).name,filesep];
    bnames = dir([sesspath,'Behaviour',filesep,'*.mat']);
    
    for b = 1:length(bnames)
        
        load([sesspath,'Behaviour',filesep,bnames(b).name])  % load behaviour
        load([sesspath,'Sample_seqs',filesep,bnames(b).name])  % load stimulus sequences
        
        % Converting sample and choice values to appropriate signs for choice regressions
        stimIn = round(stimIn.*-1);
        Cchoices = Behav(:,2)-1;
        
        % Convert stimulus values to LLRs
        if size(stimIn,1)>length(Cchoices), stimIn = stimIn(1:length(Cchoices),:); end  % trimming stimIn trials in case they exceed .mat trials (happens if block was terminated prematurely)
        CLLRin = log(normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2))./normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
        
        % Concatenate
        LLRin = [LLRin; CLLRin(~isnan(Behav(:,3)),:)];
        choices = [choices; Cchoices(~isnan(Behav(:,3)),:)];
    end
end

for t = 1:size(LLRin,1)  % Calculating number of samples per trial, to be passed on to objective function calculation for fast indexing
    nsamps(t,1) = length(find(~isnan(LLRin(t,:))));
end

% Defining PSO options (see pso_Trelea_vectorized.m for details)
P(1)=0;  P(2)=1500;     P(3)=300;    P(4:13)=[1.6 1.9 0.9 0.4 400 1e-25 250 NaN 0 1];
% display  n_iterations  n_particles       acceleration, inertia, tolerance, etc

% Seeding first n particles with parameters drawn from realistic distributions
n_seeded = 200;
PSOseedValue=[];

PSOseedValue(1:n_seeded,end+1) = seeds.H(1)+(randn(n_seeded,1).*seeds.H(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.H(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.H(1)),end) = range_p.H(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.H(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.H(2)),end) = range_p.H(2); end

PSOseedValue(1:n_seeded,end+1) = seeds.noise(1)+(randn(n_seeded,1).*seeds.noise(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.noise(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.noise(1)),end) = range_p.noise(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.noise(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.noise(2)),end) = range_p.noise(2); end

PSOseedValue(1:n_seeded,end+1) = seeds.incon(1)+(randn(n_seeded,1).*seeds.incon(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.incon(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.incon(1)),end) = range_p.incon(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.incon(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.incon(2)),end) = range_p.incon(2); end

for p = 1:length(rLLR)
    PSOseedValue(1:n_seeded,end+1) = seeds.LLRn(p,1)+(randn(n_seeded,1).*seeds.LLRn(p,2));
    if ~isempty(find(PSOseedValue(:,end)<range_p.LLRn(p,1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.LLRn(p,1)),end) = range_p.LLRn(p,1); end
    if ~isempty(find(PSOseedValue(:,end)>range_p.LLRn(p,2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.LLRn(p,2)),end) = range_p.LLRn(p,2); end
end

% Concatenating parameter ranges
par_range = [range_p.H; range_p.noise; range_p.incon; range_p.LLRn];

% Randomly seeding remaining particles within prespecified bounds
PSOseedValue(size(PSOseedValue,1)+1:P(3),1:size(PSOseedValue,2)) = normmat(rand([P(3)-n_seeded,size(PSOseedValue,2)]),par_range',1);

% Running PSO routine
[output,te,tr] = pso_Trelea_vectorized_Glaze('Glaze_cross_entropy_fitting_npLLR_InconUp',length(mv),mv,par_range,0,P,'goplotpso',PSOseedValue);

% Store parameter estimates
pm_fit = output(1:end-1);
err = output(end);

% Save generative stats and fitted parameters for this iteration
save([savepath,subj,'_fixed_fit.mat'],'pm_fit','rLLR','err','gen','P','PSOseedValue','seeds','range_p','te','tr')

clear global

