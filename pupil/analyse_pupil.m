clear; clc
close all

% Path stuff
loadpath = 'C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Analysis\Pupil\3.Interpolated\';
behavpath = 'C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Data\';
addpath(genpath('C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Analysis\Pupil'));
addpath('C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Analysis\Behaviour')
addpath('C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Analysis\Gen_fun')
addpath('C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Simulations')
addpath('C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Analysis\MEG\Scripts')

allsubj = {'DHB','TFD','EXF','JTB','TNB','QNV','PDP','GSB','OMF','NIF','ECB','TSJ','KSV','HBC','EMB','DCB','EXG'};
sess    = [3 3 3 3 3 4 3 3 3 3 3  3 3 3 3 3 3];
nblocks = [8 8 8 7 8 8 8 9 8 9 8  9 8 8 8 9 8;
           8 9 9 7 9 9 9 9 8 9 9 10 9 9 8 9 9;
           9 9 9 9 9 9 9 9 9 8 9  9 8 8 9 9 9;
           0 0 0 0 0 7 0 0 0 0 0  0 0 0 0 0 0];
               
nsubj = length(allsubj);

% Analysis stuff
maxsamps = 12;
model_type = 'Glaze_basic_InconUp';

freqs = [0.06 6];  % filter cutoffs [lo hi]
newFs = 50;  % new sampling rate if desired (set to empty if no resampling required)

basewin = [-0.05 0.05];  % window for baselining, relative to eliciting event (default = [-0.05 0.05])
fullwin = [-2 7];  % window for plotting full dilation response aligned to trial onset
sampwin = [-0.2 1.6];  % window for plotting dilation response to indiviudal samples (default = [0 1.0])
fbwin = [-0.5 4];    % window for plotting dilation response to feeback

first_deriv = 1;      % take first derivative of pupil signal rather than original time course (will not impose any baseline)
singlesamp_base = 0;  % baseline at trialwise or samplewise levels for regressions
outliercut = inf;     % |z|-threshold for throwing out observations from regression models

plot_type = 'av';     % plot either single-subjects ('ss') or grand-averages ('av')
plot_ss_scatters = 1; % generate per-subject, per-sample scatterplots

if ~first_deriv
    dilwin = [1.2 1.5];  % window for measuring mean/peak sample-wise pupil dilation
else dilwin = [0.3 0.9];
end

% Loop through subjects
for subj = 1:nsubj
    dil_full=[]; dil_samp=[]; dil_samp_peak=[]; dil_fbC=[]; dil_fbE=[];
    X_samp=[]; Y_samp=[]; base_samp=[]; dil_delay=[];
    oLLR_full=[]; LLR_full=[]; LPR_full=[]; psi_full=[]; LPR_final=[]; RT_full=[]; acc_full=[]; accN_full=[];
    surprise_full=[]; surpriseS_full=[]; absLsurprise_full=[]; Gsurprise_full=[]; pCP_full=[]; PE_full=[]; choices_full=[]; switch_full=[];
    
    % loading parameter estimates if requested
    if strcmp(model_type,'Glaze_basic')
        fitpath = 'C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Analysis\Modelling\Glaze_basic\Fits\';
        load([fitpath,allsubj{subj},'_fixed_fit.mat']);
        GA_H(subj,1) = pm_fit(1);
        GA_B(subj,1) = pm_fit(2);
    elseif strcmp(model_type,'Glaze_basic_InconUp')
        fitpath = 'C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation\Analysis\Modelling\Glaze_basic_InconUp\Fits\';
        load([fitpath,allsubj{subj},'_fixed_fit.mat']);
        GA_H(subj,1) = pm_fit(1);
        GA_B(subj,1) = pm_fit(2);
        GA_incon(subj,1) = pm_fit(4);
    end
    
    % loop through sessions
    for s = 1:sess(subj)
        fprintf('Subj %s, session %d\n',allsubj{subj},s)  % print progress
        if s==1 GApupilPSD{subj}=[]; end
        for b = 1:nblocks(s,subj)
            % Load behaviour and sample sequences
            load([behavpath,allsubj{subj},filesep,'S',num2str(s),filesep,'Behaviour',filesep,allsubj{subj},'_',num2str(s),'_',num2str(b),'.mat'])
            load([behavpath,allsubj{subj},filesep,'S',num2str(s),filesep,'Sample_seqs',filesep,allsubj{subj},'_',num2str(s),'_',num2str(b),'.mat'])
            
            % Converting sample and choice values to appropriate signs for choice regressions
            stimIn = round(stimIn.*-1);
            choices = Behav(:,2)-1;
            
            % Convert stimulus values to LLRs
            LLRinO = log(normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2))./normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
            LLRin = LLRinO.*pm_fit(2);
            
            % Calculate p(stim)
            if strcmp(model_type,'normative')
                pIn = cat(3,normpdf(stimIn,gen.mu(1),gen.sigma(1)),normpdf(stimIn,gen.mu(2),gen.sigma(2)));  % make matrix of sample probabilities for both right (:,:,1) and left (:,:,2) distributions
            else
                csmp=1; while LLRin(1,csmp)==0, csmp = csmp+1; end
                cLLR = LLRinM(1,csmp); cStim = stimIn(1,csmp);
                sigma = sqrt((-(cStim-17).^2 + (cStim+17).^2)./(2.*cLLR));
                pIn = cat(3,normpdf(stimIn,17,sigma),normpdf(stimIn,-17,sigma));
            end
            
            if strcmp(model_type,'normative')
                [LPRout,surprise,scaled_prior] = accGlaze_fast(LLRin,gen.H,0,'DY_prior_weighted',pIn);
                [~,Gsurprise] = accGlaze_fast(LLRin,gen.H,0,'conditional',pIn);
                [~,pCP] = accGlaze_fast(LLRin,gen.H,0,'pCP',pIn);
            elseif strcmp(model_type,'Glaze_basic_InconUp')
                [LPRout,surprise,scaled_prior] = accGlaze_InconUp_fast(LLRin,pm_fit(1),GA_incon(subj),0,'DY_prior_weighted',pIn);
                [~,Gsurprise] = accGlaze_InconUp_fast(LLRin,pm_fit(1),GA_incon(subj),0,'conditional',pIn);
                [~,pCP] = accGlaze_InconUp_fast(LLRin,pm_fit(1),GA_incon(subj),0,'pCP',pIn);
            else
                [LPRout,surprise,scaled_prior] = accGlaze_fast(LLRin,pm_fit(1),0,'DY_prior_weighted',pIn);
                [~,Gsurprise] = accGlaze_fast(LLRin,pm_fit(1),0,'conditional',pIn);
                [~,pCP] = accGlaze_fast(LLRin,pm_fit(1),0,'pCP',pIn);
            end
            pCP = log(pCP);
                        
            % Load eyetracker file
            load([loadpath,allsubj{subj},'_',num2str(s),'_',num2str(b),'_interp.mat']);
            pupil = zscore(data.pupil);  % z-scoring within-block
            times = data.times;
            if length(data.event(:,1))>length(choices), data.event = data.event(1:length(choices),:); end  % trimming EL trials in case they exceed .mat trials (can happen if block was terminated prematurely)
            
            % Downsampling EL data (speeds processing and aligns all datasets to same Fs - some were not recorded @ desired 1000Hz)
            pupil = resample(pupil,newFs,data.fsample)';
            data.Xgaze = resample(data.Xgaze,newFs,data.fsample)';    % X-GAZE regressor
            data.Ygaze = resample(data.Ygaze,newFs,data.fsample)';    % Y-GAZE regressor
            times = (0:(length(pupil)-1))./newFs; data.times = times;  % manually creating new times vector
            
            data.event(:,2:4) = round(data.event(:,2:4).*(newFs/data.fsample));
            data.eventsmp = round(data.eventsmp.*(newFs/data.fsample));
            data.badsmp = unique(round(data.badsmp.*(newFs/data.fsample)));  % log of all bad samples that were previously interpolated
            data.badsmp(data.badsmp==0) = [];  % in case a sample index was rounded to zero
            
            data.fsample = newFs;  % replacing stored sampling rate in data structure
            
            % Ensuring behavioural and pupil datasets contain same trials
            nsampstest=[];
            for t = 1:length(choices)
                nsampstest(t,1) = sum(~isnan(stimIn(t,:)));
            end
            assert(sum(nsampstest-data.event(:,1))==0,'Sample mismatch for subject %s, session %d, block %d.',allsubj{subj},s,b);
            if sum(nsampstest-data.event(:,1))>0, blah, end
            
            % Initialize times vectors
            if subj==1
                fulltimes = times(times<=diff(fullwin))+fullwin(1);
                samptimes = times(times<=diff(sampwin))+sampwin(1);
                fbtimes = times(times<=diff(fbwin))+fbwin(1);
            end
            
            % Isolating useable full-sequence trials
            ts=[];
            for t = 1:length(choices)
                if sum(isnan(stimIn(t,:)))==0 && choices(t)<2, ts(end+1) = t; end   %  && Behav(t,6)==0
            end
            
            % Collating useable single trials
            if ~strcmp(model_type,'normative'), oLLR_full = [oLLR_full; LLRinO(ts,:)]; end  % storing actual LLRs for PPI analyses if using model-estimated subjective LLRs elsewhere
            LLR_full = [LLR_full; LLRin(ts,:)];
            LPR_full = [LPR_full; LPRout(ts,:)];
            psi_full = [psi_full; scaled_prior(ts,:)];
            LPR_final = [LPR_final; LPRout(ts,end)];
            surprise_full = [surprise_full; surprise(ts,:)];
            Gsurprise_full = [Gsurprise_full; Gsurprise(ts,:)];
            pCP_full = [pCP_full; pCP(ts,:)];
            choices_full = [choices_full; choices(ts)];
            
            % Filter
            [bfilt,afilt] = butter(3, freqs(1)*2/data.fsample, 'high');   % hi-pass
            pupil = filtfilt(bfilt,afilt, pupil);
            
            [bfilt,afilt] = butter(3, freqs(2)*2/data.fsample, 'low');   % lo-pass
            pupil = filtfilt(bfilt,afilt, pupil);
            
            pupil = zscore(pupil); rawpupil = pupil;
            if first_deriv, pupil = diff(pupil).*data.fsample; end  % taking 1st derivative of pupil signal if desired, in z/s
            
            % Loop through trials
            for t = ts
                % Pull full trial response
                smp1 = find(times>=times(data.event(t,2))+fullwin(1),1,'first');
                if ~first_deriv
                    dil_full(end+1,:) = pupil(smp1:smp1+length(fulltimes)-1)-mean(pupil(times>=times(data.event(t,2))+basewin(1) & times<=times(data.event(t,2))+basewin(2)));
                else dil_full(end+1,:) = pupil(smp1:smp1+length(fulltimes)-1);  % no baseline if first derivative
                end
                % Pull individual sample response
                csmp = size(dil_samp,1)+1;
                if ~singlesamp_base, cbase = mean(pupil(times>=times(data.event(t,2))+basewin(1) & times<=times(data.event(t,2))+basewin(2))); end
                for smp = 1:size(data.eventsmp,2)
                    smp1 = find(times>=times(data.eventsmp(t,smp))+sampwin(1),1,'first');
                    if ~first_deriv
                        if singlesamp_base
                            dil_samp(csmp,smp,:) = pupil(smp1:smp1+length(samptimes)-1)-mean(pupil(times>=times(data.eventsmp(t,smp))+basewin(1) & times<=times(data.eventsmp(t,smp))+basewin(2)));
                        else dil_samp(csmp,smp,:) = pupil(smp1:smp1+length(samptimes)-1)-cbase;
                        end
                    else dil_samp(csmp,smp,:) = pupil(smp1:smp1+length(samptimes)-1);
                    end
                end
                for smp = 1:size(data.eventsmp,2)
                    smp1 = find(times>=times(data.eventsmp(t,smp))+sampwin(1),1,'first');
                    X_samp(csmp,smp,:) = data.Xgaze(smp1:smp1+length(samptimes)-1);
                    Y_samp(csmp,smp,:) = data.Ygaze(smp1:smp1+length(samptimes)-1);
                    
                    base_samp(csmp,smp) = mean(rawpupil(times>=times(data.eventsmp(t,smp))+basewin(1) & times<=times(data.eventsmp(t,smp))+basewin(2)));  % sample-wise measure of baseline pupil
                end
                % Pull feedback response
                smp1 = find(times>=times(data.event(t,4))+fbwin(1),1,'first');
                if smp1+length(fbtimes)-1 <= length(pupil)  % making sure epoch lies within range
                    if data.event(t,6)==1
                        dil_fbC(end+1,:) = pupil(smp1:smp1+length(fbtimes)-1)-mean(pupil(times>=times(data.event(t,4))+basewin(1) & times<=times(data.event(t,4))+basewin(2)));
                    elseif data.event(t,6)==0
                        dil_fbE(end+1,:) = pupil(smp1:smp1+length(fbtimes)-1)-mean(pupil(times>=times(data.event(t,4))+basewin(1) & times<=times(data.event(t,4))+basewin(2)));
                    end
                end
            end
        end
    end
    % Collating trial-averaged responses
    GAdil_full(subj,:) = mean(dil_full,1);
    GAdil_samp(subj,:,:) = mean(dil_samp,1);
    GAdil_fbC(subj,:) = mean(dil_fbC,1);
    GAdil_fbE(subj,:) = mean(dil_fbE,1);
            
    % Running sample-wise multiple regressions with pCP, |LLR| & -|psi| as predictors
    pupil_r_full = {}; pupil_r_reg1={}; pupil_r_reg1_nobase={}; r_full_ts={};
    for i = 1:length(samptimes)
        ts = find_inliers([squeeze(dil_samp(:,:,i)) pCP_full(:,2:end) psi_full(:,2:end)],outliercut);
        for smp = 1:maxsamps-1
            Rs = [nanzscore(pCP_full(:,smp+1)) nanzscore(abs(LLR_full(:,smp+1))) nanzscore(-abs(psi_full(:,smp+1)))];
            if smp-1<1
                Rs_reg1 = [nanzscore(pCP_full(:,smp+1)) nanzscore(abs(LLR_full(:,smp+1))) nanzscore(-abs(psi_full(:,smp+1)))];
            else
                Rs_reg1 = [nanzscore(pCP_full(:,smp+1)) nanzscore(abs(LLR_full(:,smp+1))) nanzscore(-abs(psi_full(:,smp+1))) nanzscore(pCP_full(:,smp)) nanzscore(abs(LLR_full(:,smp))) nanzscore(-abs(psi_full(:,smp)))];
            end
            smpmodel = regstats(nanzscore(dil_samp(ts,smp+1,i)),[Rs(ts,:) X_samp(ts,smp+1,i) Y_samp(ts,smp+1,i)],'linear',{'beta','r','tstat'});
            GA_dil_B_LPR_LLR_full(subj,smp,i,1:3) = smpmodel.beta(2:4);
            pupil_r_full{i}(:,smp) = nanzscore(smpmodel.r);
            
            smpmodel = regstats(nanzscore(dil_samp(ts,smp+1,i)),[Rs_reg1(ts,:) nanzscore(X_samp(ts,smp+1,i)) nanzscore(Y_samp(ts,smp+1,i))],'linear',{'beta','r','tstat'});
            GA_dil_B_LPR_LLR_full_reg1(subj,smp,i,1:3) = smpmodel.beta(2:4);
            pupil_r_reg1{i}(:,smp) = nanzscore(smpmodel.r);
            
            smpmodel = regstats(nanzscore(dil_samp(ts,smp+1,i)),[Rs_reg1(ts,:) nanzscore(X_samp(ts,smp+1,i)) nanzscore(Y_samp(ts,smp+1,i)) nanzscore(base_samp(ts,smp+1))],'linear',{'beta','r','tstat'}); % BASELINE-COVARIATE model
            GA_dil_B_LPR_LLR_full_reg1_nobase(subj,smp,i,1:3) = smpmodel.beta(2:4);
            pupil_r_reg1_nobase{i}(:,smp) = nanzscore(smpmodel.r);
        end
        r_full_ts{i} = ts;
    end
        
    % Running sample-wise choice regressions with LLR+pupil/surprise residual interaction terms (i.e. PPI) (ALWAYS USING TRUE, NOT MODEL-ESTIMATED, LLRs)
    for i = 1:length(samptimes)  % PPIs with LLR*surprise interaction terms also included in models
        B = glmfit([nanzscore(oLLR_full(r_full_ts{i},:)) nanzscore(oLLR_full(r_full_ts{i},2:end).*nanzscore(pCP_full(r_full_ts{i},2:end))) nanzscore(oLLR_full(r_full_ts{i},2:end).*nanzscore(-abs(psi_full(r_full_ts{i},2:end)))) ...
            nanzscore(oLLR_full(r_full_ts{i},2:end).*squeeze(pupil_r_full{i}),0,1)],[choices_full(r_full_ts{i}) ones(length(r_full_ts{i}),1)],'binomial');
        GA_dil_surp_PPI_Bfull(subj,:,i) = B((end-size(pupil_r_full{i},2)+1):end);
    end
    
    for i = 1:length(samptimes)  % PPIs with LLR*surprise interaction terms also included in models
        B = glmfit([nanzscore(oLLR_full(r_full_ts{i},:)) nanzscore(oLLR_full(r_full_ts{i},2:end).*nanzscore(pCP_full(r_full_ts{i},2:end))) nanzscore(oLLR_full(r_full_ts{i},2:end).*nanzscore(-abs(psi_full(r_full_ts{i},2:end)))) ...
            nanzscore(oLLR_full(r_full_ts{i},2:end).*squeeze(pupil_r_reg1{i}),0,1)],[choices_full(r_full_ts{i}) ones(length(r_full_ts{i}),1)],'binomial');
        GA_dil_surp_PPI_Breg1(subj,:,i) = B((end-size(pupil_r_reg1{i},2)+1):end);
    end
    
    for i = 1:length(samptimes)  % PPIs with LLR*surprise interaction terms also included in models, using residuals from BASELINE-COVARIATE model
        B = glmfit([nanzscore(oLLR_full(r_full_ts{i},:)) nanzscore(oLLR_full(r_full_ts{i},2:end).*nanzscore(pCP_full(r_full_ts{i},2:end))) nanzscore(oLLR_full(r_full_ts{i},2:end).*nanzscore(-abs(psi_full(r_full_ts{i},2:end)))) ...
            nanzscore(oLLR_full(r_full_ts{i},2:end).*squeeze(pupil_r_reg1_nobase{i}),0,1)],[choices_full(r_full_ts{i}) ones(length(r_full_ts{i}),1)],'binomial');
        GA_dil_surp_PPI_Breg1_nobase(subj,:,i) = B((end-size(pupil_r_reg1_nobase{i},2)+1):end);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PLOTTING %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting average trial-evoked and sample-evoked pupil responses
smarks = repmat(0.4:0.4:0.4*12,2,1);
sampcols = [linspace(0.8,0,size(GAdil_samp,2))' linspace(0.8,0,size(GAdil_samp,2))' linspace(0.8,0,size(GAdil_samp,2))'];

figure,
subplot(1,2,1), hold on,
lx = line(fullwin,[0 0]); set(lx,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])

se = std(GAdil_full,[],1)./sqrt(size(GAdil_full,1));
shadedErrorBar(fulltimes,mean(GAdil_full,1),se,{'Color',[0 0 1],'LineWidth',1.5},0);
for subj = 1:size(GAdil_full,1)
    plot(fulltimes,GAdil_full(subj,:),'Color',[0.7 0.7 0.7],'LineWidth',0.75)
end
plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
l = line([smarks],repmat(get(gca, 'ylim')',1,size(smarks,2))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
xlim(fullwin), xlabel('Time relative to trial onset (s)'), ylabel('Pupil dilation (z)')
set(gca,'TickDir','out','box','off')

subplot(1,2,2), hold on,
lx = line(sampwin,[0 0]); set(lx,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
for samp = 1:size(GAdil_samp,2)
    plot(samptimes,squeeze(mean(GAdil_samp(:,samp,:),1)),'Color',sampcols(samp,:));
end
plot([0 0],get(gca, 'ylim')','Color',[0 0 0],'LineWidth',0.5);
l = line(repmat(0.4:0.4:2.0,2,1),repmat(get(gca, 'ylim')',1,length(0.4:0.4:2.0))); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
xlim(sampwin), xlabel('Time relative to sample onset (s)'), ylabel('Pupil dilation (z)')
set(gca,'TickDir','out','box','off')


% Plotting sample-averaged encoding of pCP, |LLR| & -|psi| from single multiple regression
figure,
subplot(1,2,1), hold on,
seS = std(squeeze(mean(GA_dil_B_LPR_LLR_full(:,:,:,1),2)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full,1));
seL = std(squeeze(mean(GA_dil_B_LPR_LLR_full(:,:,:,2),2)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full,1));
seP = std(squeeze(mean(GA_dil_B_LPR_LLR_full(:,:,:,3),2)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full,1));
seSpeak = std(squeeze(mean(GA_dilPeak_B_LPR_LLR_full(:,:,1),2)),[],1)./sqrt(size(GA_dilPeak_B_LPR_LLR_full,1));
seLpeak = std(squeeze(mean(GA_dilPeak_B_LPR_LLR_full(:,:,2),2)),[],1)./sqrt(size(GA_dilPeak_B_LPR_LLR_full,1));
sePpeak = std(squeeze(mean(GA_dilPeak_B_LPR_LLR_full(:,:,3),2)),[],1)./sqrt(size(GA_dilPeak_B_LPR_LLR_full,1));

shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full(:,:,:,1),2),1)),seS,{'Color',[1 0 0],'LineWidth',1.5},0);
shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full(:,:,:,2),2),1)),seL,{'Color',[0 1 0],'LineWidth',1.5},0);
shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full(:,:,:,3),2),1)),seP,{'Color',[0 0 0],'LineWidth',1.5},0);
for samp = 1:size(GA_dil_absLsurp_B,2)
    p1=plot(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full(:,:,:,1),2),1)),'Color',[1 0 0],'LineWidth',1.5);
    p2=plot(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full(:,:,:,2),2),1)),'Color',[0 1 0],'LineWidth',1.5);
    p3=plot(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full(:,:,:,3),2),1)),'Color',[0 0 0],'LineWidth',1.5);
end

plot([0.9 0.9],[squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full(:,:,1),2),1))-seSpeak squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full(:,:,1),2),1))+seSpeak],'Color',[0.6 0.6 0.6])
plot([0.95 0.95],[squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full(:,:,2),2),1))-seLpeak squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full(:,:,2),2),1))+seLpeak],'Color',[0.6 0.6 0.6])
plot([1.0 1.0],[squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full(:,:,3),2),1))-sePpeak squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full(:,:,3),2),1))+sePpeak],'Color',[0.6 0.6 0.6])
S=scatter(0.9,squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full(:,:,1),2),1)),60); set(S,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0])
S=scatter(0.95,squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full(:,:,2),2),1)),60); set(S,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0])
S=scatter(1.0,squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full(:,:,3),2),1)),60); set(S,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])

l = line(sampwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--')
xlim([0 sampwin(2)]), ylabel('\beta'), xlabel('Time relative to sample onset (s)'), if first_deriv, ylim([-0.02 0.09]), else ylim([-0.04 0.04]), end
legend([p1 p2 p3],{'Surprise_k','|LLR|_k','|LPR|_k_-_1'}), set(gca,'TickDir','out')

if first_deriv, avwin = [0.45 0.75]; else avwin = [0.5 1.6]; end
av_dil_surp_B = squeeze(mean(mean(GA_dil_B_LPR_LLR_full(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3),1));
seS = std(squeeze(mean(GA_dil_B_LPR_LLR_full(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full,1));
if first_deriv, avwin = [0.5 0.8]; else avwin = [0.7 1.1]; end
av_dil_deltaL_B = squeeze(mean(mean(GA_dil_B_LPR_LLR_full(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),2),3),1));
seL = std(squeeze(mean(GA_dil_B_LPR_LLR_full(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),2),3)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full,1));
if first_deriv, avwin = [0.05 0.35]; else avwin = [0.1 0.8]; end
av_dil_prior_B = squeeze(mean(mean(GA_dil_B_LPR_LLR_full(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),3),3),1));
seP = std(squeeze(mean(GA_dil_B_LPR_LLR_full(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),3),3)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full,1));
subplot(1,2,2), hold on,
shadedErrorBar(2:12,av_dil_prior_B,seP,{'Color',[0 0 0],'LineStyle','--'},0)
shadedErrorBar(2:12,av_dil_deltaL_B,seL,{'Color',[0 1 0],'LineStyle','--'},0)
shadedErrorBar(2:12,av_dil_surp_B,seS,{'Color',[1 0 0],'LineStyle','--'},0)
plot(2:12,av_dil_deltaL_B,'Color',[0 1 0],'LineStyle','--');
plot(2:12,av_dil_surp_B,'Color',[1 0 0],'LineStyle','--');
l = line([0 12],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
for samp = 2:12
    S2=scatter(samp,av_dil_prior_B(samp-1)); set(S2,'MarkerEdgeColor',sampcols(samp-1,:));
    S=scatter(samp,av_dil_deltaL_B(samp-1)); set(S,'MarkerEdgeColor',sampcols(samp-1,:));
    S1=scatter(samp,av_dil_surp_B(samp-1)); set(S1,'MarkerEdgeColor',sampcols(samp-1,:));
end
xlim([1 12.5]),  if first_deriv, ylim([-0.04 0.15]), else ylim([-0.04 0.04]), end, ylabel('\beta'), xlabel('Sample position')
set(gca,'TickDir','out')

% Plotting sample-averaged encoding of pCP, |LLR| & -|psi| from single multiple regression WITH PREVIOUS-SAMPLE VARS AS COVARIATES
nperm = 10000; clustalpha = 0.05; alpha = 0.05;
disp('Running cluster-based permutation test #1')
[sig_ts1,sig_ts_uncorr1,cP1,~] = cluster_permWS_fast(cat(3,squeeze(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,1),2)),zeros(size(GA_dil_B_LPR_LLR_full_reg1,1),size(GA_dil_B_LPR_LLR_full_reg1,3))),nperm,clustalpha,alpha);
disp('Running cluster-based permutation test #2')
[sig_ts2,sig_ts_uncorr2,cP2,~] = cluster_permWS_fast(cat(3,squeeze(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,2),2)),zeros(size(GA_dil_B_LPR_LLR_full_reg1,1),size(GA_dil_B_LPR_LLR_full_reg1,3))),nperm,clustalpha,alpha);
disp('Running cluster-based permutation test #3')
[sig_ts3,sig_ts_uncorr3,cP3,~] = cluster_permWS_fast(cat(3,squeeze(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,3),2)),zeros(size(GA_dil_B_LPR_LLR_full_reg1,1),size(GA_dil_B_LPR_LLR_full_reg1,3))),nperm,clustalpha,alpha);

figure,
subplot(1,2,1), hold on,
seS = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,1),2)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1,1));
seL = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,2),2)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1,1));
seP = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,3),2)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1,1));
seSpeak = std(squeeze(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,1),2)),[],1)./sqrt(size(GA_dilPeak_B_LPR_LLR_full_reg1,1));
seLpeak = std(squeeze(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,2),2)),[],1)./sqrt(size(GA_dilPeak_B_LPR_LLR_full_reg1,1));
sePpeak = std(squeeze(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,3),2)),[],1)./sqrt(size(GA_dilPeak_B_LPR_LLR_full_reg1,1));

shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,1),2),1)),seS,{'Color',[1 0 0],'LineWidth',1},0);
shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,2),2),1)),seL,{'Color',[0 1 0],'LineWidth',1},0);
shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,3),2),1)),seP,{'Color',[0 0 0],'LineWidth',1},0);
for samp = 1:size(GA_dil_absLsurp_B,2)
    p1=plot(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,1),2),1)),'Color',[1 0 0],'LineWidth',1);
    p2=plot(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,2),2),1)),'Color',[0 1 0],'LineWidth',1);
    p3=plot(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,:,3),2),1)),'Color',[0 0 0],'LineWidth',1);
end

plot([0.9 0.9],[squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,1),2),1))-seSpeak squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,1),2),1))+seSpeak],'Color',[0.6 0.6 0.6])
plot([0.95 0.95],[squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,2),2),1))-seLpeak squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,2),2),1))+seLpeak],'Color',[0.6 0.6 0.6])
plot([1.0 1.0],[squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,3),2),1))-sePpeak squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,3),2),1))+sePpeak],'Color',[0.6 0.6 0.6])
S=scatter(0.9,squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,1),2),1)),60); set(S,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0])
S=scatter(0.95,squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,2),2),1)),60); set(S,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0])
S=scatter(1.0,squeeze(mean(mean(GA_dilPeak_B_LPR_LLR_full_reg1(:,:,3),2),1)),60); set(S,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])

plot(samptimes,sig_ts1.*-0.02,'Color',[1 0 0],'LineWidth',3.5)
plot(samptimes,sig_ts2.*-0.024,'Color',[0 1 0],'LineWidth',3.5)
plot(samptimes,sig_ts3.*-0.028,'Color',[0 0 0],'LineWidth',3.5)
l = line(sampwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
xlim([0 sampwin(2)]), ylabel('\beta'), xlabel('Time relative to sample onset (s)'), if first_deriv, ylim([-0.03 0.1]), else ylim([-0.04 0.04]), end
legend([p1 p2 p3],{'surprise_k','|LLR|_k','-|psi|_k'}), set(gca,'TickDir','out')

if first_deriv, avwin = [0.45 0.75]; else avwin = [0.5 1.6]; end
av_dil_surp_B = squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3),1));
seS = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1,1));
if first_deriv, avwin = [0.3 0.5]; else avwin = [0.7 1.1]; end
av_dil_deltaL_B = squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),2),3),1));
seL = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),2),3)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1,1));
if first_deriv, avwin = [0.2 0.4]; else avwin = [0.1 0.8]; end
av_dil_prior_B = squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),3),3),1));
seP = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),3),3)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1,1));
subplot(1,2,2), hold on,
shadedErrorBar(2:12,av_dil_prior_B,seP,{'Color',[0 0 0],'LineStyle','--'},0)
shadedErrorBar(2:12,av_dil_deltaL_B,seL,{'Color',[0 1 0],'LineStyle','--'},0)
shadedErrorBar(2:12,av_dil_surp_B,seS,{'Color',[1 0 0],'LineStyle','--'},0)
plot(2:12,av_dil_deltaL_B,'Color',[0 1 0],'LineStyle','--');
plot(2:12,av_dil_surp_B,'Color',[1 0 0],'LineStyle','--');
l = line([0 12],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
for samp = 2:12
    S2=scatter(samp,av_dil_prior_B(samp-1)); set(S2,'MarkerEdgeColor',sampcols(samp-1,:));
    S=scatter(samp,av_dil_deltaL_B(samp-1)); set(S,'MarkerEdgeColor',sampcols(samp-1,:));
    S1=scatter(samp,av_dil_surp_B(samp-1)); set(S1,'MarkerEdgeColor',sampcols(samp-1,:));
end
xlim([1 12.5]),  if first_deriv, ylim([-0.04 0.15]), else ylim([-0.04 0.04]), end, ylabel('\beta'), xlabel('Sample position')
set(gca,'TickDir','out')


% Plotting sample-averaged encoding of pCP, |LLR| & -|psi| from single multiple regression with previous sample vars AND BASELINE PUPIL as covariates
disp('Running cluster-based permutation test #1 - nobase')
[sig_ts1_nobase,sig_ts_uncorr1_nobase,cP1_nobase,~] = cluster_permWS_fast(cat(3,squeeze(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,1),2)),zeros(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1),size(GA_dil_B_LPR_LLR_full_reg1_nobase,3))),nperm,clustalpha,alpha);
disp('Running cluster-based permutation test #2 - nobase')
[sig_ts2_nobase,sig_ts_uncorr2_nobase,cP2_nobase,~] = cluster_permWS_fast(cat(3,squeeze(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,2),2)),zeros(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1),size(GA_dil_B_LPR_LLR_full_reg1_nobase,3))),nperm,clustalpha,alpha);
disp('Running cluster-based permutation test #3 - nobase')
[sig_ts3_nobase,sig_ts_uncorr3_nobase,cP3_nobase,~] = cluster_permWS_fast(cat(3,squeeze(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,3),2)),zeros(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1),size(GA_dil_B_LPR_LLR_full_reg1_nobase,3))),nperm,clustalpha,alpha);

figure,
subplot(1,2,1), hold on,
seS = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,1),2)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1));
seL = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,2),2)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1));
seP = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,3),2)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1));

shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,1),2),1)),seS,{'Color',[1 0 0],'LineWidth',1},0);
shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,2),2),1)),seL,{'Color',[0 1 0],'LineWidth',1},0);
shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,3),2),1)),seP,{'Color',[0 0 0],'LineWidth',1},0);
for samp = 1:size(GA_dil_absLsurp_B,2)
    p1=plot(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,1),2),1)),'Color',[1 0 0],'LineWidth',1);
    p2=plot(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,2),2),1)),'Color',[0 1 0],'LineWidth',1);
    p3=plot(samptimes,squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,:,3),2),1)),'Color',[0 0 0],'LineWidth',1);
end

plot(samptimes,sig_ts1_nobase.*-0.02,'Color',[1 0 0],'LineWidth',3.5)
plot(samptimes,sig_ts2_nobase.*-0.024,'Color',[0 1 0],'LineWidth',3.5)
plot(samptimes,sig_ts3_nobase.*-0.028,'Color',[0 0 0],'LineWidth',3.5)
l = line(sampwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0.6 0.6 0.6])
xlim([0 sampwin(2)]), ylabel('\beta'), xlabel('Time relative to sample onset (s)'), if first_deriv, ylim([-0.03 0.1]), else ylim([-0.04 0.04]), end
legend([p1 p2 p3],{'surprise_k','|LLR|_k','-|psi|_k'}), set(gca,'TickDir','out')

if first_deriv, avwin = [0.45 0.75]; else avwin = [0.5 1.6]; end
av_dil_surp_B = squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3),1));
seS = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),1),3)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1));
if first_deriv, avwin = [0.3 0.5]; else avwin = [0.7 1.1]; end
av_dil_deltaL_B = squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),2),3),1));
seL = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),2),3)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1));
if first_deriv, avwin = [0.2 0.4]; else avwin = [0.1 0.8]; end
av_dil_prior_B = squeeze(mean(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),3),3),1));
seP = std(squeeze(mean(GA_dil_B_LPR_LLR_full_reg1_nobase(:,:,samptimes>=avwin(1) & samptimes<=avwin(2),3),3)),[],1)./sqrt(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1));
subplot(1,2,2), hold on,
shadedErrorBar(2:12,av_dil_prior_B,seP,{'Color',[0 0 0],'LineStyle','--'},0)
shadedErrorBar(2:12,av_dil_deltaL_B,seL,{'Color',[0 1 0],'LineStyle','--'},0)
shadedErrorBar(2:12,av_dil_surp_B,seS,{'Color',[1 0 0],'LineStyle','--'},0)
plot(2:12,av_dil_deltaL_B,'Color',[0 1 0],'LineStyle','--');
plot(2:12,av_dil_surp_B,'Color',[1 0 0],'LineStyle','--');
l = line([0 12],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
for samp = 2:12
    S2=scatter(samp,av_dil_prior_B(samp-1)); set(S2,'MarkerEdgeColor',sampcols(samp-1,:));
    S=scatter(samp,av_dil_deltaL_B(samp-1)); set(S,'MarkerEdgeColor',sampcols(samp-1,:));
    S1=scatter(samp,av_dil_surp_B(samp-1)); set(S1,'MarkerEdgeColor',sampcols(samp-1,:));
end
xlim([1 12.5]),  if first_deriv, ylim([-0.04 0.15]), else ylim([-0.04 0.04]), end, ylabel('\beta'), xlabel('Sample position')
set(gca,'TickDir','out')


% Plotting results of PPI analyses
if first_deriv, avwin = [0.35 0.55]; else avwin = [0.5 1.2]; end
smps = 2:12;

nperm = 10000; clustalpha = 0.05; alpha = 0.05;
testtimes = [0 1.0];  % times within which to run cluster-based permutation test
disp('Running PPI cluster-based permutation test #1')
[sig_ts1,sig_ts_uncorr1,cP1,~] = cluster_permWS_fast(cat(3,squeeze(mean(GA_dil_surp_PPI_Breg1(:,smps-1,samptimes>=testtimes(1) & samptimes<=testtimes(2)),2)),...
                                                zeros(size(GA_dil_B_LPR_LLR_full_reg1,1),length(find(samptimes>=0 & samptimes<=1)))),nperm,clustalpha,alpha);
disp('Running PPI cluster-based permutation test #1 - nobase')
[sig_ts1_nobase,sig_ts_uncorr1_nobase,cP1_nobase,~] = cluster_permWS_fast(cat(3,squeeze(mean(GA_dil_surp_PPI_Breg1_nobase(:,smps-1,samptimes>=testtimes(1) & samptimes<=testtimes(2)),2)),...
                                                zeros(size(GA_dil_B_LPR_LLR_full_reg1_nobase,1),length(find(samptimes>=0 & samptimes<=1)))),nperm,clustalpha,alpha);

sigtimes = samptimes(samptimes>=testtimes(1) & samptimes<=testtimes(2));


figure,  %  RESIDUALS FROM UNIVARIATE REGRESSIONS

subplot(2,3,2), hold on, %  RESIDUALS FROM MULTIVARIATE REGRESSIONS
se = std(squeeze(mean(GA_dil_surp_PPI_Bfull(:,smps-1,:),2)),[],1)./sqrt(size(GA_dil_surp_PPI_Bfull,1));
sePeak = std(squeeze(mean(GA_dilPeak_surp_PPI_Bfull,2)),[],1)./sqrt(size(GA_dilPeak_surp_PPI_Bfull,1));
l = line(sampwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--')
shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_surp_PPI_Bfull(:,smps-1,:),1),2)),se,{'Color',[1 0 0],'LineWidth',1.5},0);
for samp = 1:size(GA_dil_surp_PPI_Bfull,2)
    plot(samptimes,squeeze(mean(GA_dil_surp_PPI_Bfull(:,samp,:),1)),'Color',sampcols(samp,:));
end
plot([0.9 0.9],[mean(mean(GA_dilPeak_surp_PPI_Bfull))-sePeak mean(mean(GA_dilPeak_surp_PPI_Bfull))+sePeak],'Color',[0.6 0.6 0.6])
S=scatter(0.9,mean(mean(GA_dilPeak_surp_PPI_Bfull)),60); set(S,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0])
xlim([0 1]), ylim([-0.055 0.13]), ylabel('Regression coefficient (a.u.)'), xlabel('Time relative to sample onset (s)')
set(gca,'TickDir','out')

subplot(2,3,6), hold on,
av_dil_surp_B = squeeze(mean(mean(GA_dil_surp_PPI_Bfull(:,:,samptimes>=avwin(1) & samptimes<=avwin(2)),3),1));
se = std(squeeze(mean(GA_dil_surp_PPI_Bfull(:,:,samptimes>=avwin(1) & samptimes<=avwin(2)),3)),[],1)./sqrt(size(GA_dil_surp_PPI_Bfull,1));
shadedErrorBar(2:12,av_dil_surp_B,se,{'Color',[1 0 0],'LineStyle','--'},0)
l = line([0 12],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
for samp = 2:12
    S=scatter(samp,av_dil_surp_B(samp-1)); set(S,'MarkerEdgeColor',sampcols(samp-1,:))
end
ylabel('Regression coefficient (a.u.)'), xlabel('Sample position'), xlim([1 12.5])%, ylim([-0.005 0.08])
title('Residuals from model including |LLR| & |LPR|-1')
set(gca,'TickDir','out')

subplot(2,3,3), hold on, %  RESIDUALS FROM MULTIVARIATE REGRESSIONS WITH PREVIOUS-SAMPLE VARS AS COVARIATES                                            
se = std(squeeze(mean(GA_dil_surp_PPI_Breg1(:,smps-1,:),2)),[],1)./sqrt(size(GA_dil_surp_PPI_Breg1,1));
sePeak = std(squeeze(mean(GA_dilPeak_surp_PPI_Breg1,2)),[],1)./sqrt(size(GA_dilPeak_surp_PPI_Breg1,1));
l = line(sampwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--')
shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_surp_PPI_Breg1(:,smps-1,:),1),2)),se,{'Color',[1 0 0],'LineWidth',1.5},0);
for samp = 1:size(GA_dil_surp_PPI_Breg1,2)
    plot(samptimes,squeeze(mean(GA_dil_surp_PPI_Breg1(:,samp,:),1)),'Color',sampcols(samp,:));
end
plot([0.9 0.9],[mean(mean(GA_dilPeak_surp_PPI_Breg1))-sePeak mean(mean(GA_dilPeak_surp_PPI_Breg1))+sePeak],'Color',[0.6 0.6 0.6])
S=scatter(0.9,mean(mean(GA_dilPeak_surp_PPI_Breg1)),60); set(S,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0])
plot(sigtimes,sig_ts1.*-0.04,'LineWidth',3,'Color',[0 0 0])
xlim([0 1]), ylim([-0.055 0.13]), ylabel('Regression coefficient (a.u.)'), xlabel('Time relative to sample onset (s)')
set(gca,'TickDir','out')

subplot(2,3,7), hold on,
av_dil_surp_B = squeeze(mean(mean(GA_dil_surp_PPI_Breg1(:,:,samptimes>=avwin(1) & samptimes<=avwin(2)),3),1));
se = std(squeeze(mean(GA_dil_surp_PPI_Breg1(:,:,samptimes>=avwin(1) & samptimes<=avwin(2)),3)),[],1)./sqrt(size(GA_dil_surp_PPI_Breg1,1));
shadedErrorBar(2:12,av_dil_surp_B,se,{'Color',[1 0 0],'LineStyle','--'},0)
l = line([0 12],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
for samp = 2:12
    S=scatter(samp,av_dil_surp_B(samp-1)); set(S,'MarkerEdgeColor',sampcols(samp-1,:))
end
ylabel('Regression coefficient (a.u.)'), xlabel('Sample position'), xlim([1 12.5])%, ylim([-0.005 0.08])
title('Residuals from model including |LLR| & |LPR|-1')
set(gca,'TickDir','out')

subplot(2,3,4), hold on, %  RESIDUALS FROM MULTIVARIATE REGRESSIONS WITH PREVIOUS-SAMPLE VARS & BASELINE PUPIL AS COVARIATES                                            
se = std(squeeze(mean(GA_dil_surp_PPI_Breg1_nobase(:,smps-1,:),2)),[],1)./sqrt(size(GA_dil_surp_PPI_Breg1_nobase,1));
l = line(sampwin,[0 0]); set(l,'LineWidth',0.5,'LineStyle','--')
shadedErrorBar(samptimes,squeeze(mean(mean(GA_dil_surp_PPI_Breg1_nobase(:,smps-1,:),1),2)),se,{'Color',[1 0 0],'LineWidth',1.5},0);
for samp = 1:size(GA_dil_surp_PPI_Breg1_nobase,2)
    plot(samptimes,squeeze(mean(GA_dil_surp_PPI_Breg1_nobase(:,samp,:),1)),'Color',sampcols(samp,:));
end
plot(sigtimes,sig_ts1_nobase.*-0.04,'LineWidth',3,'Color',[0 0 0])
xlim([0 1]), ylim([-0.055 0.13]), ylabel('Regression coefficient (a.u.)'), xlabel('Time relative to sample onset (s)')
set(gca,'TickDir','out')

subplot(2,3,8), hold on,
av_dil_surp_B = squeeze(mean(mean(GA_dil_surp_PPI_Breg1_nobase(:,:,samptimes>=avwin(1) & samptimes<=avwin(2)),3),1));
se = std(squeeze(mean(GA_dil_surp_PPI_Breg1_nobase(:,:,samptimes>=avwin(1) & samptimes<=avwin(2)),3)),[],1)./sqrt(size(GA_dil_surp_PPI_Breg1_nobase,1));
shadedErrorBar(2:12,av_dil_surp_B,se,{'Color',[1 0 0],'LineStyle','--'},0)
l = line([0 12],[0 0]); set(l,'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
for samp = 2:12
    S=scatter(samp,av_dil_surp_B(samp-1)); set(S,'MarkerEdgeColor',sampcols(samp-1,:))
end
ylabel('Regression coefficient (a.u.)'), xlabel('Sample position'), xlim([1 12.5])%, ylim([-0.005 0.08])
title('Residuals from model including |LLR| & |LPR|-1 & RawBase')
set(gca,'TickDir','out')



