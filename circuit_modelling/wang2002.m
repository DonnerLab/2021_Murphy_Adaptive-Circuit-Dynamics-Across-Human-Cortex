function [choice,sw_relFR,events] = wang2002(decvar,seedval,w)

% [winning_pop,reaction_time,Firing_Rates,time_list] = wang2002(decvar,seedval,w)
% 
% decvar contains decision parameters - meanA and meanB must be specified
%   (Wang2002 has meanA and meanB averaging to 40Hz)
% seedval (optional) is seed for random number generator, to exactly
%   reproduce model results 
% w is a pre-specified weight matrix - otherwise this gets built (a little
%   slowly) in the script

%% input arguments
if nargin<1
    error('Not enough input arguments');
end

if nargin>1 %set random seed using input
    if ~isempty(seedval)
        if ischar(seedval); seedval = str2num(seedval);end
        seedval = uint32(seedval);
        st = RandStream('mt19937ar','Seed',seedval);
        RandStream.setGlobalStream(st);
    end
else
    seedval = nan;
end

meanR = decvar.meanR;  %input when LLR = 0 (i.e. symmetric point around which non-zero LLRs vary)
LLR = decvar.LLR; %vector of LLRs
c = decvar.c; %multiplicative scaling factor applied to LLRs to determine effective input
ms = decvar.ms; %decision length
dt = decvar.dt; %time step
N = decvar.N; %Ncells
f = decvar.f; %fraction selective
vext = decvar.vext; %external noise input
stim_start = decvar.stim_start; %stim on time (ms)                
samp_dur = decvar.samp_dur; %stim off time (ms)                

%% 
tic
fprintf(1,'\nInitializing Parameters... ') ;

time_list = dt:dt:ms;
events = [stim_start:samp_dur:(samp_dur*length(LLR))];  % vector marking sample onsets

stim_start = stim_start / dt;
samp_dur = samp_dur / dt;
stim_end =  stim_start + samp_dur*length(LLR);

tstop = round(ms / dt) ;
Ne = .8 * N ;   %number excitatory
Ni = .2 * N ;   %number inhibitory

%indexing:
ind.Pool1 = 1:(Ne*f) ;
ind.Pool2 = (Ne*f+1):(Ne*f*2) ;
ind.Ns = (Ne*f*2+1):Ne ;
ind.Inh = (Ne+1):N ;

toc

%%
if nargin<3
    fprintf(1, '\nBuilding Weight Matrix... ') ;
    w = zeros(N) ;  %weights matrix - all to all
    wplus = 1.7 ;
    wminus = 1 - f * (wplus - 1)/(1 - f) ;
    
    for i=1:N
        for j=1:N
            if i==j % No Self Connections
                w(i,j) = 0 ;
            elseif any(i == ind.Pool1) && any(j == ind.Pool1) % Same group
                w(i,j) = wplus ;
            elseif any(i == ind.Pool2) && any(j == ind.Pool2) % Same group
                w(i,j) = wplus ;
            elseif any(i == ind.Pool1) && any(j == ind.Pool2)  %% Different Group
                w(i,j) = wminus ;
            elseif any(i == ind.Pool2) && any(j == ind.Pool1) %% Different Group
                w(i,j) = wminus ;
            elseif any(i == union(ind.Pool1,ind.Pool2)) && any(j == ind.Ns) %Nonselective
                w(i,j) = wminus ;
            else %all other connections
                w(i,j) = 1;
            end
        end
    end
    toc
    if 0 %this is how to create the sparse matrix
        indr = randsample(N^2,N^2/2,false); %randomly sample half the connections
        w(indr) = 0; %set these to 0
        for i = 1:N; 
            w(:,i) = w(:,i)./mean(w(:,i)); %normalise so mean input is still 1
        end
    end
        
else
    fprintf(1, '\nUsing pre-loaded weight matrix... ') ;
end

if ~all(size(w)==N)||~ndims(w)==2 %check weight matrix size
    error('w must be Ncells*Ncells matrix')
end

%%
fprintf(1, '\nGenerating Stimulus... ') ; 

Stim = zeros(N,tstop) ;
breaks = [1 stim_start stim_start+samp_dur*(1:length(LLR)) tstop];

LLRin1 = (meanR + LLR*c)./1000.*dt; LLRin1(LLRin1<0)=0;
LLRin2 = (meanR - LLR*c)./1000.*dt; LLRin2(LLRin2<0)=0;

Stim(:,1:breaks(2)) = poissrnd(vext*dt,N,stim_start) ;  % pre-onset input
for smp = 1:length(LLR)
    Stim(1:floor(Ne*f),breaks(smp+1)+1:breaks(smp+2)) = poissrnd(LLRin1(smp)+vext*dt,floor(Ne*f),samp_dur) ; %pop 1 input
    Stim(floor(Ne*f)+1:floor(2*Ne*f),breaks(smp+1)+1:breaks(smp+2)) = poissrnd(LLRin2(smp)+vext*dt,floor(Ne*f),samp_dur) ; %pop 2 input
    Stim(floor(2*Ne*f)+1:end,breaks(smp+1)+1:breaks(smp+2)) = poissrnd(vext*dt,N-floor(2*Ne*f),samp_dur) ;  % non-selective pop input
end
Stim(:,breaks(smp+2):tstop) = poissrnd(vext*dt,N,length(breaks(smp+2):tstop)) ;  % post-stim input
toc


%% Run simulation
fprintf(1, '\nStarting Simulation... ') ;
spikes = spikeloop(dt, tstop, Ne, Ni, w, Stim);
toc
clear Stim

%% Average
fprintf(1, '\nAveraging Firing Rates... ') ;

window = 50 / dt  ;
fn = fieldnames(ind) ;
Firing_Rates = zeros(tstop,length(fn)) ;

for i = 1:length(fn)
    Firing_Rates(:,i) = filter(1/window*ones(1,window),1,mean(spikes(ind.(fn{i}),:),1)) * 1000 / dt;
end
toc

%% Determine winning population and compute sample-wise relative activation
t = 400; % latency (ms) post-sample onset at which to measure relative activation

relFR = Firing_Rates(:,1)-Firing_Rates(:,2);
for smp = 1:length(events)
    sw_relFR(1,smp) = relFR(find(round(time_list,2)==(events(smp)+t)));
end

if relFR(stim_end)>0 % at stim_offset, firing rate of 1>2, so population A won
    choice = 1;
elseif relFR(stim_end)<0 % population B won
    choice = 0;
end
