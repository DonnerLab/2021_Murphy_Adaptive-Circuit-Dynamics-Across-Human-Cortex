% Load behavioural data and variety of model fits, and estimate desired
% behavioural measures (accuracies, psychophysical kernels, etc.).
% For illustration, this example script deals with some of the more complex
% model fits from paper, including models in which 'non-parametric' L->psi
% and stimulus->LLR transfer functions were estimated. Many analyses in the
% paper deal with simpler models. See paper for details.

clear, close all

dirstr = 'C:\Users\Peter UKE\Desktop\Experiments\Surprise_accumulation';  % depending on which pc

datpath = [dirstr,'\Data\'];   % behavioural data
npfitpath = [dirstr,'\Analysis\Modelling\FitPsi_data_npLLR_InconUp\Fits\'];  % Model fit w/ non-parametric transfer function, non-param LLR, and inconsistency weighting
normfitpath = [dirstr,'\Analysis\Modelling\Glaze_npLLR_InconUp\Fits\'];  % Model fit w/ fitted H, non-param LLR, and inconsistency weighting
perffitpath = [dirstr,'\Analysis\Modelling\Perf_basic\Fits\'];   % Model fit w/ perfect integration

addpath([dirstr,'\Simulations'])
addpath(genpath([dirstr,'\Analysis\Gen_fun']))
addpath ('C:\Users\Peter UKE\Desktop\Experiments\Tools')

allsubj = {'DHB','EXF','TFD','JTB','TNB','QNV','PDP','GSB','OMF','NIF','ECB','KSV','TSJ','HBC','EMB','DCB','EXG'};

Htrue = 0.08;
nsmps = 12;

% Loop through subjects
for subj = 1:length(allsubj)
    oLLR_full=[]; LPR_full=[];
    surprise_full=[]; choices_full=[]; fdist_full=[]; acc_full=[]; fsmps_full=[]; choices_m1_full=[]; fdist_m1_full=[]; choicesCE=[];
    LPR_full_norm=[]; LPR_full_ideal=[]; LPR_full_perf=[]; surprise_full_norm=[]; surprise_full_ideal=[]; surprise_full_perf=[]; psi_full=[]; psi_full_norm=[]; psi_full_ideal=[]; psi_full_perf=[];
    LPR_final_fsmps=[]; LPR_final_norm_fsmps=[]; LPR_final_ideal_fsmps=[]; LPR_final_perf_fsmps=[]; LLRfinal_fsmps=[]; LPR_perfacc_fsmps=[];
    
    fprintf('Processing subject %s...\n',allsubj{subj})
    
    % Load GLAZE model fits for comparison of fitted psi function with corresponding H
    load([normfitpath,allsubj{subj},'_fixed_fit.mat'])
    GA_H(subj,1) = pm_fit(1);
    GA_noise_norm(subj,1) = pm_fit(2);
    GA_incon_norm(subj,1) = pm_fit(3);
        
    LLRs = log(normpdf(floor(gen.range(1)):1:ceil(gen.range(2)),gen.mu(1),gen.sigma(1))./normpdf(floor(gen.range(1)):1:ceil(gen.range(2)),gen.mu(2),gen.sigma(2)));
    gains = cumsum(pm_fit(4:end)');        % reconstructing fitted LLR gain terms
    gains_norm = [fliplr(gains) 0 gains];  % creating full vector of gain terms, symmetric around LLR=0
    rLLRfull_norm = [sort(-rLLR) 0 rLLR];  % creating full vector of matching original LLRs
    
    GA_LLR_norm(subj,1:length(rLLRfull_norm)) = rLLRfull_norm.*gains_norm;   % storing actual parameter estimates, translated from point-wise gain parameters to point-wise LLRs
    tLLRs_norm{subj} = interp1(rLLRfull_norm,rLLRfull_norm.*gains_norm,LLRs,'spline');
    % One awkward side-effect of 'non-parametrically estimated' stimulus->LLR
    % functions is that, unlike other model variants with a simple linear scaling
    % of the function, they don't allow us to easily compute the probability
    % of a stimulus each location, p(stim). Since this quantity is needed
    % for computation of CPP, we here fit the p(stim) distribution under
    % the hard constraint that it's integral over the full stimulus space equals 1.
    LLRext = log(normpdf(-100:1:100,gen.mu(1),gen.sigma(1))./normpdf(-100:1:100,gen.mu(2),gen.sigma(2)));  % go beyond stim range to avoid edge artifacts
    xy = [interp1(LLRext,-100:1:100,rLLRfull_norm,'linear')' (rLLRfull_norm.*gains_norm)'];  % [x=true stimulus positions, y=fitted LLRs]
    log_pstim = fminsearchbnd(@(params) fit_pstim_from_LLR(params,xy), log(normpdf(xy(:,1),gen.mu(1),gen.sigma(1))), ones(1,size(xy,1)).*-60, ones(1,size(xy,1)).*-1);  % fit p(stim) in log space
    stim=xy(:,1)'; pstim = exp(interp1(stim,log_pstim,[min(stim) -90:1:90 max(stim)],'spline')); pstim = pstim./sum(pstim);  % calculate normalized hi-res p(stim)
    GA_pstim_norm(subj,:) = pstim;
        
    % Load perfect accumulator model fits
    load([perffitpath,allsubj{subj},'_fixed_fit.mat'])
    
    GA_B_perf(subj,1) = pm_fit(1);
    GA_noise_perf(subj,1) = pm_fit(2);
    
    tLLRs_perf{subj} = LLRs.*GA_B_perf(subj);
            
    % Load model fit with non-parametrically estimated L->psi transfer function & store fitted parameters
    load([npfitpath,allsubj{subj},'_fixed_fit.mat'])
    
    GA_psi_n(subj,1:(length(L_nm1)*2 + 1)) = [fliplr(-pm_fit(1:length(L_nm1))') 0 pm_fit(1:length(L_nm1))'];
    GA_noise(subj,1) = pm_fit(end-1);
    GA_incon(subj,1) = pm_fit(end);
    GA_err(subj,1) = err;
    
    gains = cumsum(pm_fit(length(L_nm1)+1:length(L_nm1)+length(rLLR))');
    gains = [fliplr(gains) 0 gains];  % creating full vector of gain terms, symmetric around LLR=0
    rLLRfull = [sort(-rLLR) 0 rLLR];  % creating full vector of matching original LLRs
    
    GA_LLR_n(subj,1:length(rLLRfull)) = rLLRfull.*gains;   % storing actual parameter estimates, translated from point-wise gain parameters to point-wise LLRs
    tLLRs{subj} = interp1(rLLRfull,rLLRfull.*gains,LLRs,'spline');
    if subj==1, xLLR = floor(gen.range(1)):1:ceil(gen.range(2)); end
    
    LLRext = log(normpdf(-100:1:100,gen.mu(1),gen.sigma(1))./normpdf(-100:1:100,gen.mu(2),gen.sigma(2))); % Fit p(stim) from fitted LLRs under constraint that integral of p(stim) within stimulus bounds = 1
    xy = [interp1(LLRext,-100:1:100,rLLRfull,'linear')' (rLLRfull.*gains)'];  % [x=true stimulus positions, y=fitted LLRs]
    log_pstim = fminsearchbnd(@(params) fit_pstim_from_LLR(params,xy), log(normpdf(xy(:,1),gen.mu(1),gen.sigma(1))), ones(1,size(xy,1)).*-60, ones(1,size(xy,1)).*-1);  % fit p(stim) in log space
    stim=xy(:,1)'; pstim = exp(interp1(stim,log_pstim,[min(stim) -90:1:90 max(stim)],'spline')); pstim = pstim./sum(pstim);  % calculate normalized hi-res p(stim)
    GA_pstim(subj,:) = pstim;
    
    % Load behavioural data & calculate model-derived belief updates
    sdirs = dir([datpath,allsubj{subj},filesep,'S*']);
    for s = 1:length(sdirs)
        
        sesspath = [datpath,allsubj{subj},filesep,sdirs(s).name,filesep];
        bnames = dir([sesspath,'Behaviour',filesep,'*.mat']);
        
        for b = 1:length(bnames)
            % Load behaviour and sample sequences
            load([sesspath,'Behaviour',filesep,bnames(b).name])
            load([sesspath,'Sample_seqs',filesep,bnames(b).name])
            
            % Converting sample and choice values to appropriate signs for choice regressions
            stimIn = round(stimIn.*-1);
            choices = Behav(:,2)-1;
            acc = Behav(:,3);
            
            % Convert stimulus polar angles to LLRs & run belief updating
            LLRinO = log(normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2))./normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
            pIn_ideal = cat(3,normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2)),normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
            
            LLRinM = nan(size(LLRinO)); pIn = nan(size(pIn_ideal));
            for samp = 1:size(LLRinO,2)
                LLRinM(:,samp) = interp1(rLLRfull,rLLRfull.*gains,LLRinO(:,samp),'spline');
                pIn(:,samp,1) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(GA_pstim(subj,:)),stimIn(:,samp),'spline'));  % always interpolating in log space
                pIn(:,samp,2) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],fliplr(log(GA_pstim(subj,:))),stimIn(:,samp),'spline'));
            end
            
            LLRinM_norm = nan(size(LLRinO)); pIn_norm = nan(size(pIn_ideal));
            for samp = 1:size(LLRinO,2)
                LLRinM_norm(:,samp) = interp1(rLLRfull_norm,rLLRfull_norm.*gains_norm,LLRinO(:,samp),'spline');
                pIn_norm(:,samp,1) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],log(GA_pstim_norm(subj,:)),stimIn(:,samp),'spline'));  % always interpolating in log space
                pIn_norm(:,samp,2) = exp(interp1([min(xy(:,1)) -90:1:90 max(xy(:,1))],fliplr(log(GA_pstim_norm(subj,:))),stimIn(:,samp),'spline'));
            end
                        
            LLRinM_perf = LLRinO.*GA_B_perf(subj);
            cLLR = LLRinM_perf(1,5); cStim = stimIn(1,5);
            sigma = sqrt((-(cStim-17).^2 + (cStim+17).^2)./(2.*cLLR));
            pIn_perf = cat(3,normpdf(stimIn,17,sigma),normpdf(stimIn,-17,sigma));
            
            [LPRout,surprise,psi] = accPhi_InconUp_fast(LLRinM,[sort(-L_nm1) 0 L_nm1],GA_psi_n(subj,:),'pCP',pIn,GA_H(subj),GA_incon(subj));
            psi = psi(:,1:end-1);  % dumping transformed prior after last sample
            [LPRout_norm,surprise_norm,psi_norm] = accGlaze_InconUp_fast(LLRinM_norm,GA_H(subj),GA_incon_norm(subj),0,'pCP',pIn_norm);
            [LPRout_ideal,surprise_ideal,psi_ideal] = accGlaze_fast(LLRinO,Htrue,0,'pCP',pIn_ideal);
            [LPRout_perf,surprise_perf,psi_perf] = accPerf_fast(LLRinM_perf,0,'pCP',pIn_perf,GA_H(subj));
                        
            surprise = log(surprise);
            surprise_norm = log(surprise_norm);
            surprise_ideal = log(surprise_ideal);
            surprise_perf = log(surprise_perf);
                        
            % Pulling samples per trial and position of final change-point
            nsamps=[]; fCPpos=[];
            for t = 1:size(LLRinO,1)
                nsamps(t,1) = length(find(~isnan(LLRinO(t,:))));
                if isempty(find(pswitch(t,:)==1, 1))
                    fCPpos(t,1) = 1;
                else fCPpos(t,1) = find(pswitch(t,:)==1,1,'last');
                end
            end
            subinds = sub2ind(size(LPRout_ideal),1:size(LPRout_ideal,1),nsamps');
            LPRfinal = LPRout(subinds)';
            LPRfinal_norm = LPRout_norm(subinds)';
            LPRfinal_ideal = LPRout_ideal(subinds)';
            LPRfinal_perf = LPRout_perf(subinds)';
            LLRfinal = LLRinO(subinds)';
            
            % Pulling # samples from last change-point
            fsmps = (nsamps-fCPpos)+1;  % full-length no-CP trials will take value of 13
            fsmps(fsmps==13) = 12;
                        
            % Isolating useable full-sequence trials
            ts=[]; tsCE=[];
            for t = 1:length(choices)
                if sum(isnan(stimIn(t,:)))==0 && choices(t)<2, ts(end+1) = t; end   % full sequence
                if choices(t)<2, tsCE(end+1) = t; end   % all trials, regardless of sequence length
            end
            
            % Collating useable single trials
            oLLR_full = [oLLR_full; LLRinO(ts,:)];
            
            LPR_full = [LPR_full; LPRout(ts,end)];
            LPR_full_norm = [LPR_full_norm; LPRout_norm(ts,end)];
            LPR_full_ideal = [LPR_full_ideal; LPRout_ideal(ts,end)];
            LPR_full_perf = [LPR_full_perf; LPRout_perf(ts,end)];
            
            LPR_final_fsmps = [LPR_final_fsmps; LPRfinal(tsCE)];
            LPR_final_norm_fsmps = [LPR_final_norm_fsmps; LPRfinal_norm(tsCE)];
            LPR_final_ideal_fsmps = [LPR_final_ideal_fsmps; LPRfinal_ideal(tsCE)];
            LPR_final_perf_fsmps = [LPR_final_perf_fsmps; LPRfinal_perf(tsCE)];
            LLRfinal_fsmps = [LLRfinal_fsmps; LLRfinal(tsCE)];
            LPR_perfacc_fsmps = [LPR_perfacc_fsmps; nansum(LLRinO(tsCE,:),2)];
            
            surprise_full = [surprise_full; surprise(ts,:)];
            surprise_full_norm = [surprise_full_norm; surprise_norm(ts,:)];
            surprise_full_ideal = [surprise_full_ideal; surprise_ideal(ts,:)];
            surprise_full_perf = [surprise_full_perf; surprise_perf(ts,:)];
            
            psi_full = [psi_full; psi(ts,:)];
            psi_full_norm = [psi_full_norm; psi_norm(ts,:)];
            psi_full_ideal = [psi_full_ideal; psi_ideal(ts,:)];
            psi_full_perf = [psi_full_perf; psi_perf(ts,:)];
            
            choices_full = [choices_full; choices(ts)];
            fdist_full = [fdist_full; Behav(tsCE,1)-1];
            acc_full = [acc_full; acc(tsCE)];
            fsmps_full = [fsmps_full; fsmps(tsCE)];
            if ~isempty(ts)
                if ts(1)==1
                    choices_m1_full = [choices_m1_full; [nan; choices(ts(2:end)-1)]];
                    fdist_m1_full = [fdist_m1_full; [nan; Behav(ts(2:end)-1,1)]];
                else choices_m1_full = [choices_m1_full; choices(ts-1)];
                    fdist_m1_full = [fdist_m1_full; Behav(ts-1,1)];
                end
            end
            
       end
    end
        
    % Calculate average accuracy
    maxsamps = size(oLLR_full,2);
        
    % Running LLR/SURPRISE/UNCERTAINTY regressions of choice
    [B,~,stats] = glmfit([nanzscore(oLLR_full,0,1) nanzscore(oLLR_full(:,2:end).*nanzscore(surprise_full(:,2:end)),0,1) nanzscore(oLLR_full(:,2:end).*nanzscore(-abs(psi_full(:,2:end))),0,1)],[choices_full ones(length(choices_full),1)],'binomial');
    GA_ssB_regU(subj,1:maxsamps) = B(2:maxsamps+1);
    GA_ssB_surpU(subj,:) = B(maxsamps+2:end-length(2:maxsamps));
    GA_ssB_uncertU(subj,:) = B(end-length(2:maxsamps)+1:end);
            
    % Calculating accuracy per subject & as a function of # samples from final change-point
    GAacc(subj,1) = nansum(acc_full)/length(find(~isnan(acc_full)));
    
    bnds_fsmp = [1 3; 4 6; 7 9; 10 12];  % # final sample bounds within which to average
    for b = 1:size(bnds_fsmp,1)
        GA_fsmp_acc(subj,b) = nanmean(acc_full(fsmps_full>=bnds_fsmp(b,1) & fsmps_full<=bnds_fsmp(b,2)));
        GA_fsmp_n(subj,b) = length(find(fsmps_full>=bnds_fsmp(b,1) & fsmps_full<=bnds_fsmp(b,2)));
    end
    
    fdistc = fdist_full; fdistc(fdistc==0) = -1;
    choices_c = choicesCE; choices_c(choices_c==0) = -1;
    
    cacc = LPR_final_norm_fsmps./GA_noise_norm(subj);  %%% Normative model fit
    cacc(sign(fdistc)==-1) = cacc(sign(fdistc)==-1).*-1; % conditioning LPRs on accuracy, not choice
    cacc = 0.5+0.5.*erf(cacc);
    acc_norm(subj,1) = mean(cacc);
    for b = 1:size(bnds_fsmp,1)
        acc_fsmp_norm(subj,b) = mean(cacc(fsmps_full>=bnds_fsmp(b,1) & fsmps_full<=bnds_fsmp(b,2)));
    end
    cacc = sign(LPR_final_norm_fsmps);
    acc_norm_binary(subj,1) = length(find(cacc-fdistc==0))./length(cacc);
    
    cacc = LPR_final_ideal_fsmps;    %%% Noiseless ideal observer
    cacc(sign(fdistc)==-1) = cacc(sign(fdistc)==-1).*-1; % conditioning LPRs on accuracy, not choice
    cacc = 0.5+0.5.*erf(cacc);
    acc_ideal_nn(subj,1) = mean(cacc);
    cacc = sign(LPR_final_ideal_fsmps);
    acc_ideal_nn_binary(subj,1) = length(find(cacc-fdistc==0))./length(cacc);
    for b = 1:size(bnds_fsmp,1)
        c = find(fsmps_full>=bnds_fsmp(b,1) & fsmps_full<=bnds_fsmp(b,2));
        acc_fsmp_ideal_nn(subj,b) = length(find(cacc(c)-fdistc(c)==0))./length(c);
    end
            
    cacc = LLRfinal_fsmps;       %%% Last-sample only (no noise)
    cacc(sign(fdistc)==-1) = cacc(sign(fdistc)==-1).*-1; % conditioning LPRs on accuracy, not choice
    cacc = 0.5+0.5.*erf(cacc);
    acc_fsmp(subj,1) = mean(cacc);
    cacc = sign(LLRfinal_fsmps);
    acc_fsmp_binary(subj,1) = length(find(cacc-fdistc==0))./length(cacc);
    for b = 1:size(bnds_fsmp,1)
        c = find(fsmps_full>=bnds_fsmp(b,1) & fsmps_full<=bnds_fsmp(b,2));
        acc_fsmp_fsmp(subj,b) = length(find(cacc(c)-fdistc(c)==0))./length(c);
    end
    
    cacc = LPR_perfacc_fsmps;      %%% Perfect accumulation (no noise)
    cacc(sign(fdistc)==-1) = cacc(sign(fdistc)==-1).*-1; % conditioning LPRs on accuracy, not choice
    cacc = 0.5+0.5.*erf(cacc);
    acc_perfacc(subj,1) = mean(cacc);
    cacc = sign(LPR_perfacc_fsmps);
    acc_perfacc_binary(subj,1) = length(find(cacc-fdistc==0))./length(cacc);
    for b = 1:size(bnds_fsmp,1)
        c = find(fsmps_full>=bnds_fsmp(b,1) & fsmps_full<=bnds_fsmp(b,2));
        acc_fsmp_perfacc(subj,b) = length(find(cacc(c)-fdistc(c)==0))./length(c);
    end
    
    % Calculate LLR influence on choice (via binning)
    nbins = 12;
    LLRbins = linspace(-3.64,0,(nbins/2)+1); LLRbins = [LLRbins sort(abs(LLRbins(1:end-1)))]; LLRbins([1 end]) = [-inf inf];
    LLRcounts = zeros(length(choices_full),length(LLRbins)-1);
    for b = 1:length(LLRbins)-1
        for t = 1:length(choices_full)
            LLRcounts(t,b) = length(find(oLLR_full(t,:)>=LLRbins(b) & oLLR_full(t,:)<=LLRbins(b+1)));
        end
    end
    Bz = glmfit(LLRcounts,[choices_full ones(length(choices_full),1)],'binomial','constant','off');
    GA_wLLR_ssB(subj,1:size(LLRcounts,2)) = Bz;
    
    
    % Compute kernels/LLR regressions for model fits
    CP = 0.5+0.5.*erf(LPR_full_norm./GA_noise_norm(subj));  % calculate choice probabilities scaled by noise
    CP_ideal = 0.5+0.5.*erf(LPR_full_ideal./GA_noise(subj));  % calculate choice probabilities scaled by noise
    CP_perf = 0.5+0.5.*erf(LPR_full_perf./GA_noise_perf(subj));  % calculate choice probabilities scaled by noise
    
        % Normative fit
    [B,~,stats] = glmfit([nanzscore(oLLR_full,0,1) nanzscore(oLLR_full(:,2:end).*nanzscore(surprise_full_norm(:,2:end)),0,1) nanzscore(oLLR_full(:,2:end).*nanzscore(-abs(psi_full_norm(:,2:end))),0,1)],[CP ones(length(choices),1)],'binomial');
    MssB_regU(subj,1:maxsamps) = B(2:maxsamps+1);
    MssB_surpU(subj,:) = B(maxsamps+2:end-length(2:maxsamps));
    MssB_uncertU(subj,:) = B(end-length(2:maxsamps)+1:end);
    
    B = glmfit(LLRcounts,[CP ones(length(choices),1)],'binomial','constant','off');
    MwLLR_Bs(subj,1:size(LLRcounts,2)) = B;
    
        % Ideal observer w/ human noise
    [B,~,stats] = glmfit([nanzscore(oLLR_full,0,1) nanzscore(oLLR_full(:,2:end).*nanzscore(surprise_full_ideal(:,2:end)),0,1) nanzscore(oLLR_full(:,2:end).*nanzscore(-abs(psi_full_ideal(:,2:end))),0,1)],[CP_ideal ones(length(choices),1)],'binomial');
    MssB_regU_ideal(subj,1:maxsamps) = B(2:maxsamps+1);
    MssB_surpU_ideal(subj,:) = B(maxsamps+2:end-length(2:maxsamps));
    MssB_uncertU_ideal(subj,:) = B(end-length(2:maxsamps)+1:end);
    
    B = glmfit(LLRcounts,[CP_ideal ones(length(choices),1)],'binomial','constant','off');
    MwLLR_Bs_ideal(subj,1:size(LLRcounts,2)) = B;
    
        % Perfect accumulation
    [B,~,stats] = glmfit([nanzscore(oLLR_full,0,1) nanzscore(oLLR_full(:,2:end).*nanzscore(surprise_full_perf(:,2:end)),0,1) nanzscore(oLLR_full(:,2:end).*nanzscore(-abs(psi_full_perf(:,2:end))),0,1)],[CP_perf ones(length(choices),1)],'binomial');
    MssB_regU_perf(subj,1:maxsamps) = B(2:maxsamps+1);
    MssB_surpU_perf(subj,:) = B(maxsamps+2:end-length(2:maxsamps));
    MssB_uncertU_perf(subj,:) = B(end-length(2:maxsamps)+1:end);
    
    B = glmfit(LLRcounts,[CP_perf ones(length(choices),1)],'binomial','constant','off');
    MwLLR_Bs_perf(subj,1:size(LLRcounts,2)) = B;
end


% Plotting specs
fs = 12;
lw = 1.25;
axlw = 1.25;
scatsize = 40;

% Calculating ideal-observer transfer function
Lin = -10:0.1:10;
io_psi = Lin + log(((1-Htrue)/Htrue)+exp(-Lin))-log(((1-Htrue)/Htrue)+exp(Lin));

% Calculating transfer functions from Glaze model fits
for subj = 1:length(allsubj)
    Hfit_psi(:,subj) = Lin + log(((1-GA_H(subj))/GA_H(subj))+exp(-Lin))-log(((1-GA_H(subj))/GA_H(subj))+exp(Lin));
end

% Plot L->psi transfer functions & fitted hazard rates
figure, set(1,'units','centimeters','pos', [0 0 12 15],'Color',[1 1 1]),

ga_psi=[];
xext = [min(Lin) max(Lin)];   % get extent of horizontal midline

s1=subplot(2,1,1); hold on
plot(xext,[0 0],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.75)
plot([0 0],[-5 5],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.75)
for subj = 1:length(allsubj)
    plot(Lin,interp1([sort(-L_nm1) 0 L_nm1],GA_psi_n(subj,:),Lin,'spline'),'Color',[0.8 0.8 0.8],'LineWidth',0.75)
    ga_psi(:,subj) = interp1([sort(-L_nm1) 0 L_nm1],GA_psi_n(subj,:),Lin,'spline');
end
shadedErrorBar(Lin,mean(ga_psi,2),std(ga_psi,[],2)./sqrt(size(ga_psi,2)),{'Color','k','LineWidth',lw},0)
plot(Lin,io_psi,'b--','LineWidth',lw)
plot(Lin,mean(Hfit_psi,2),'r--','LineWidth',lw)
xlim(xext), ylim([-5 5])
xlabel('L_n_-_1','FontSize',19), ylabel('\psi_n','FontSize',19),
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[-10 -5 0 5 10],'YTick',[-5 0 5])

s2=subplot(2,1,2); hold on
plot([0.08 0.08],[0.5 1.5],'Color',[0 0 1],'LineWidth',2.5)
S1=scatter(sort(GA_H),ones(length(GA_H),1),scatsize); set(S1,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[1 0 0],'LineWidth',0.5)
xlim([0 0.12]), ylim([0.45 1.55]), set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','ycol','w','XTick',[0:0.04:0.12],'YTick',[])
xlabel('Fitted H','FontSize',19)

set(s1, 'Position', [0.14, 0.5, 0.75, 0.45])  % [left bottom width height
set(s2, 'Position', [0.14, 0.13, 0.75, 0.14])

% permutation tests of fitted parameters against normative values
pH = mult_comp_perm_t1(GA_H-0.08,10000);
pIU = mult_comp_perm_t1(GA_incon_norm-1,10000);



% Plot average observed and model-based kernels/surprise/UNCERTAINTY weighting functions
disp('Cluster-based permutation tests...')
nperm=10000; clustalpha=0.05; alpha=0.05;
[sig_ts1,sig_ts_uncorr1,cP1,~] = cluster_permWS_fast(cat(3,GA_ssB_regU,zeros(size(GA_ssB_regU))),nperm,clustalpha,alpha);
cbnd1=[]; s=0; con=0;
while s < length(sig_ts1), s = s+1;
    if ~con if sig_ts1(s)==1, cbnd1(end+1,1)=s; con=1; end
    else if isnan(sig_ts1(s)), cbnd1(end,2)=s-1;  con=0; elseif s==length(sig_ts1), cbnd1(end,2)=s; end
    end
end
cbnd_uncorr1=[]; s=0; con=0;
while s < length(sig_ts_uncorr1), s = s+1;
    if ~con if sig_ts_uncorr1(s)==1, cbnd_uncorr1(end+1,1)=s; con=1; end
    else if isnan(sig_ts_uncorr1(s)), cbnd_uncorr1(end,2)=s-1;  con=0; elseif s==length(sig_ts_uncorr1), cbnd_uncorr1(end,2)=s; end
    end
end

[sig_ts2,sig_ts_uncorr2,cP2,~] = cluster_permWS_fast(cat(3,GA_ssB_surpU,zeros(size(GA_ssB_surpU))),nperm,clustalpha,alpha);
cbnd2=[]; s=0; con=0;
while s < length(sig_ts2), s = s+1;
    if ~con if sig_ts2(s)==1, cbnd2(end+1,1)=s; con=1; end
    else if isnan(sig_ts2(s)), cbnd2(end,2)=s-1;  con=0; elseif s==length(sig_ts2), cbnd2(end,2)=s; end
    end
end
cbnd_uncorr2=[]; s=0; con=0;
while s < length(sig_ts_uncorr2), s = s+1;
    if ~con if sig_ts_uncorr2(s)==1, cbnd_uncorr2(end+1,1)=s; con=1; end
    else if isnan(sig_ts_uncorr2(s)), cbnd_uncorr2(end,2)=s-1;  con=0; elseif s==length(sig_ts_uncorr2), cbnd_uncorr2(end,2)=s; end
    end
end

[sig_ts3,sig_ts_uncorr3,cP3,~] = cluster_permWS_fast(cat(3,GA_ssB_uncertU,zeros(size(GA_ssB_uncertU))),nperm,clustalpha,alpha);
cbnd3=[]; s=0; con=0;
while s < length(sig_ts3), s = s+1;
    if ~con if sig_ts3(s)==1, cbnd3(end+1,1)=s; con=1; end
    else if isnan(sig_ts3(s)), cbnd3(end,2)=s-1;  con=0; elseif s==length(sig_ts3), cbnd3(end,2)=s; end
    end
end
cbnd_uncorr3=[]; s=0; con=0;
while s < length(sig_ts_uncorr3), s = s+1;
    if ~con if sig_ts_uncorr3(s)==1, cbnd_uncorr3(end+1,1)=s; con=1; end
    else if isnan(sig_ts_uncorr3(s)), cbnd_uncorr3(end,2)=s-1;  con=0; elseif s==length(sig_ts_uncorr3), cbnd_uncorr3(end,2)=s; end
    end
end

sigoffset = [0.08 0.05];

figure, set(2,'units','centimeters','pos', [0 0 36 10],'Color',[1 1 1]),
gaO = mean(GA_ssB_regU,1); seO = std(GA_ssB_regU,[],1)/sqrt(size(GA_ssB_regU,1));
gaM = mean(MssB_regU,1); seM = std(MssB_regU,[],1)/sqrt(size(MssB_regU,1));
gaMi = mean(MssB_regU_ideal,1); seMi = std(MssB_regU_ideal,[],1)/sqrt(size(MssB_regU_ideal,1));
gaMp = mean(MssB_regU_perf,1); seMp = std(MssB_regU_perf,[],1)/sqrt(size(MssB_regU_perf,1));
subplot(1,3,1), hold on
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.75)
f=fill([1:maxsamps maxsamps:-1:1],[gaMp+seMp fliplr(gaMp-seMp)],[0 0.6 0.85]); set(f,'EdgeColor','none')
f=fill([1:maxsamps maxsamps:-1:1],[gaMi+seMi fliplr(gaMi-seMi)],[0 0 1]); set(f,'EdgeColor','none')
f=fill([1:maxsamps maxsamps:-1:1],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
for p = 1:length(gaO)
    line([p p],[gaO(p)+seO(p) gaO(p)-seO(p)],'Color',[0 0 0],'LineWidth',lw)
end
plot(1:maxsamps,gaO,'Color',[0 0 0],'LineWidth',1)

for s = 1:size(cbnd_uncorr1,1), plot([cbnd_uncorr1(s,1)-0.1 cbnd_uncorr1(s,2)+0.1],[1 1].*(-0.4+diff([-0.4 3.0]).*sigoffset(1)),'Color',[0.6 0.6 0.6],'LineWidth',3), end
for s = 1:size(cbnd1,1), plot([cbnd1(s,1)-0.1 cbnd1(s,2)+0.1],[1 1].*(-0.4+diff([-0.4 3.0]).*sigoffset(2)),'Color',[0 0 0],'LineWidth',3), end

xlabel('Sample position','FontSize',19), ylabel('Weight on choice (a.u.)','FontSize',19)
xlim([0.5 maxsamps+0.5]), ylim([-0.4 3.0]) %, ylim([-0.4 2.8])
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[2:2:12],'YTick',[0 1.2 2.4])

gaO = mean(GA_ssB_surpU,1); seO = std(GA_ssB_surpU,[],1)/sqrt(size(GA_ssB_surpU,1));
gaM = mean(MssB_surpU,1); seM = std(MssB_surpU,[],1)/sqrt(size(MssB_surpU,1));
gaMi = mean(MssB_surpU_ideal,1); seMi = std(MssB_surpU_ideal,[],1)/sqrt(size(MssB_surpU_ideal,1));
gaMp = mean(MssB_surpU_perf,1); seMp = std(MssB_surpU_perf,[],1)/sqrt(size(MssB_surpU_perf,1));
subplot(1,3,2), hold on
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.75)
f=fill([2:maxsamps maxsamps:-1:2],[gaMp+seMp fliplr(gaMp-seMp)],[0 0.6 0.85]); set(f,'EdgeColor','none')
f=fill([2:maxsamps maxsamps:-1:2],[gaMi+seMi fliplr(gaMi-seMi)],[0 0 1]); set(f,'EdgeColor','none')
f=fill([2:maxsamps maxsamps:-1:2],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
for p = 1:length(gaO)
    line([p+1 p+1],[gaO(p)+seO(p) gaO(p)-seO(p)],'Color',[0 0 0],'LineWidth',lw)
end
plot(2:maxsamps,gaO,'Color',[0 0 0],'LineWidth',lw)

for s = 1:size(cbnd_uncorr2,1), plot([cbnd_uncorr2(s,1)-0.25+1 cbnd_uncorr2(s,2)+0.25+1],[1 1].*(-0.25+diff([-0.25 0.57]).*sigoffset(1)),'Color',[0.6 0.6 0.6],'LineWidth',3), end
for s = 1:size(cbnd2,1), plot([cbnd2(s,1)-0.25+1 cbnd2(s,2)+0.25+1],[1 1].*(-0.25+diff([-0.25 0.57]).*sigoffset(2)),'Color',[0 0 0],'LineWidth',3), end

xlabel('Sample position','FontSize',19), ylabel('Weight on choice (a.u.)','FontSize',19)
xlim([0.5 maxsamps+0.5]), ylim([-0.25 0.57])
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[2:2:12],'YTick',[-0.5 0 0.5 1.0 1.5])

gaO = mean(GA_ssB_uncertU,1); seO = std(GA_ssB_uncertU,[],1)/sqrt(size(GA_ssB_uncertU,1));
gaM = mean(MssB_uncertU,1); seM = std(MssB_uncertU,[],1)/sqrt(size(MssB_uncertU,1));
gaMi = mean(MssB_uncertU_ideal,1); seMi = std(MssB_uncertU_ideal,[],1)/sqrt(size(MssB_uncertU_ideal,1));
gaMp = mean(MssB_uncertU_perf,1); seMp = std(MssB_uncertU_perf,[],1)/sqrt(size(MssB_uncertU_perf,1));
subplot(1,3,3), hold on
plot([0.5 maxsamps+0.5],[0 0],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.75)
f=fill([2:maxsamps maxsamps:-1:2],[gaMp+seMp fliplr(gaMp-seMp)],[0 0.6 0.85]); set(f,'EdgeColor','none')
f=fill([2:maxsamps maxsamps:-1:2],[gaMi+seMi fliplr(gaMi-seMi)],[0 0 1]); set(f,'EdgeColor','none')
f=fill([2:maxsamps maxsamps:-1:2],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
for p = 1:length(gaO)
    line([p+1 p+1],[gaO(p)+seO(p) gaO(p)-seO(p)],'Color',[0 0 0],'LineWidth',lw)
end
plot(2:maxsamps,gaO,'Color',[0 0 0],'LineWidth',lw)

for s = 1:size(cbnd_uncorr3,1), plot([cbnd_uncorr3(s,1)-0.25+1 cbnd_uncorr3(s,2)+0.25+1],[1 1].*(-0.25+diff([-0.25 0.57]).*sigoffset(1)),'Color',[0.6 0.6 0.6],'LineWidth',3), end
for s = 1:size(cbnd3,1), plot([cbnd3(s,1)-0.25+1 cbnd3(s,2)+0.25+1],[1 1].*(-0.25+diff([-0.25 0.57]).*sigoffset(2)),'Color',[0 0 0],'LineWidth',3), end

xlabel('Sample position','FontSize',19), ylabel('Weight on choice (a.u.)','FontSize',19)
xlim([0.5 maxsamps+0.5]), ylim([-0.25 0.57])
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[2:2:12],'YTick',[-0.5 0 0.5 1.0 1.5])

% permutation test comparing surprise/uncertainty coefficients
pcomp = mult_comp_perm_t1(mean(GA_ssB_surpU,2)-mean(GA_ssB_uncertU,2),10000);


% Plot fitted stimulus polar angle -> LLR transfer functions
ga_tLLRs=[]; ga_tLLRs_norm=[];
xext = [min(xLLR)-3 max(xLLR)+3];   % get extent of horizontal midline

figure, set(4,'units','centimeters','pos', [0 0 15 11],'Color',[1 1 1]),
hold on
plot(xext,[0 0],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.75)
plot([0 0],[-17 17],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.75)
for subj = 1:length(allsubj)
    plot(xLLR,tLLRs{subj},'Color',[0.8 0.8 0.8],'LineWidth',0.75)
    ga_tLLRs(:,subj) = tLLRs{subj};
    ga_tLLRs_norm(:,subj) = tLLRs_norm{subj};
end
plot(xLLR,LLRs,'b','LineWidth',lw)
shadedErrorBar(xLLR,mean(ga_tLLRs_norm,2),std(ga_tLLRs_norm,[],2)./sqrt(size(ga_tLLRs_norm,2)),{'Color','r','LineWidth',lw},0)
shadedErrorBar(xLLR,mean(ga_tLLRs,2),std(ga_tLLRs,[],2)./sqrt(size(ga_tLLRs,2)),{'Color','k','LineWidth',lw},0)
xlim(xext), ylim([-17 17])
xlabel('Dot polar angle'), ylabel('LLR'), set(gca,'TickDir','out','box','off','XTick',[-90 -45 0 45 90],'YTick',[-20 -10 0 10 20]),
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[-90 -45 0 45 90],'YTick',[-15 0 15])


% Plotting LLR weights on choice after binning and calculating LLR counts
figure, set(5,'units','centimeters','pos', [0 0 10 9],'Color',[1 1 1]),

gaO = mean(GA_wLLR_ssB,1); seO = std(GA_wLLR_ssB,[],1)/sqrt(size(GA_wLLR_ssB,1));
gaM = mean(MwLLR_Bs,1); seM = std(MwLLR_Bs,[],1)/sqrt(size(MwLLR_Bs,1));
gaMi = mean(MwLLR_Bs_ideal,1); seMi = std(MwLLR_Bs_ideal,[],1)/sqrt(size(MwLLR_Bs_ideal,1));
binpos = [LLRbins(2:end-1)-diff(LLRbins(2:3)/2) LLRbins(end-1)+diff(LLRbins(2:3)/2)];

hold on
plot([0 0],[-1.2 1.2],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.75)
plot([binpos(1)-0.5 binpos(end)+0.5],[0 0],'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.75)
f=fill([binpos fliplr(binpos)],[gaMi+seMi fliplr(gaMi-seMi)],[0 0 1]); set(f,'EdgeColor','none')
f=fill([binpos fliplr(binpos)],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
for p = 1:length(gaO)
    line([binpos(p) binpos(p)],[gaO(p)+seO(p) gaO(p)-seO(p)],'Color',[0 0 0],'LineWidth',lw)
end
plot(binpos,gaO,'Color',[0 0 0],'LineWidth',lw)
xlabel('LLR'), ylabel('Weight on choice (a.u.)')
set(gca,'TickDir','out','box','off'), xlim([binpos(1)-0.5 binpos(end)+0.5]), ylim([-1.34 1.34])
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[-2 0 2],'YTick',[-1 0 1])



% Choice accuracy as a function of # samples from last change-point
binc = mean(bnds_fsmp,2)';

figure, hold on
gaO = mean(GA_fsmp_acc,1); seO = std(GA_fsmp_acc,[],1)/sqrt(size(GA_fsmp_acc,1));
gaM = mean(acc_fsmp_norm,1); seM = std(acc_fsmp_norm,[],1)/sqrt(size(acc_fsmp_norm,1));
gaMinn = mean(acc_fsmp_ideal_nn,1); seMinn = std(acc_fsmp_ideal_nn,[],1)/sqrt(size(acc_fsmp_ideal_nn,1));
gaMp = mean(acc_fsmp_perfacc,1); seMp = std(acc_fsmp_perfacc,[],1)/sqrt(size(acc_fsmp_perfacc,1));
gaMf = mean(acc_fsmp_fsmp,1); seMf = std(acc_fsmp_fsmp,[],1)/sqrt(size(acc_fsmp_fsmp,1));

hold on
f=fill([binc fliplr(binc)],[gaMf+seMf fliplr(gaMf-seMf)],[1 0.8 1]); set(f,'EdgeColor','none')
f=fill([binc fliplr(binc)],[gaMp+seMp fliplr(gaMp-seMp)],[0.8 1 0.8]); set(f,'EdgeColor','none')
f=fill([binc fliplr(binc)],[gaMinn+seMinn fliplr(gaMinn-seMinn)],[0.8 0.8 1]); set(f,'EdgeColor','none')
f=fill([binc fliplr(binc)],[gaM+seM fliplr(gaM-seM)],[1 0 0]); set(f,'EdgeColor','none')
for p = 1:length(binc)
    line([binc(p) binc(p)],[gaO(p)+seO(p) gaO(p)-seO(p)],'Color',[0 0 0],'LineWidth',lw)
end
plot(binc,gaO,'Color',[0 0 0],'LineWidth',lw)
xlabel('# samples from last change-point'), ylabel('Accuracy'), xlim([1 12])
set(gca,'FontName','Arial','TickDir','out','box','off','XTick',[2 4 6 8 10 12])


figpath = 'C:\Users\Peter UKE\Desktop\Submissions\Volatility_MEG_pupil\Figures\FigS2\';
save([figpath,'accXduration.mat'],'bnds_fsmp','GA_fsmp_acc',...
    'acc_fsmp_norm','acc_fsmp_ideal','acc_fsmp_incon','acc_fsmp_ideal_nn','acc_fsmp_perfacc','acc_fsmp_fsmp')


