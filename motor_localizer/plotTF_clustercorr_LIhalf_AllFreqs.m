clear, close all

% Specify inits
loadpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/TF_Motor/dB_common/';
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun

allsubj = {'DCB' 'DHB' 'EMB' 'EXF' 'HBC' 'JTB' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB'};

clustalpha = 0.01;  % % intial threhsold for forming clusters
nperms = 1000;
freqbnd = [1 100];  % frequency range to consider
timewin = [0.7 1.1];  % time range within which to run stats 

% Retrieve full gradiometer structure
load([loadpath,allsubj{1},'_output_LIhalf.mat'],'grad_full');
grad_full.chanpos = grad_full.chanpos(1:274,:);
grad_full.chanori = grad_full.chanori(1:274,:);
grad_full.chantype = grad_full.chantype(1:274);
grad_full.chanunit = grad_full.chanunit(1:274);
grad_full.label = grad_full.label(1:274);
fullgrad = grad_full;

%%%%%%%%%%%%%%%%%%
%%% TRIAL-WISE %%%
%%%%%%%%%%%%%%%%%%
all_trl_avg=[]; all_trl_avg_full=[];
allB_respF={}; gaB_respF=[];
fprintf('Loading trial-wise data...\nSubject ')

for s = 1:length(allsubj)
    fprintf('%d, ',s)
    % Load data
    load([loadpath,allsubj{s},'_output_LIhalf.mat'])
    freqsL = round(freqs);
    
    trl_avgL = trl_avg(:,freqsL>=freqbnd(1) & freqsL<=freqbnd(2),:);
    trl_avg_fullL = trl_avg_full(:,freqsL>=freqbnd(1) & freqsL<=freqbnd(2),:);
    TrespFL = TrespF(:,freqsL>=freqbnd(1) & freqsL<=freqbnd(2),:);
    
    load([loadpath,allsubj{s},'_output_LIhalf_HiFreq.mat'])
    allfreqs = [freqsL round(freqs)];
    
    trl_avg = trl_avg(:,round(freqs)>=freqbnd(1) & round(freqs)<=freqbnd(2),:);
    trl_avg_full = trl_avg_full(:,round(freqs)>=freqbnd(1) & round(freqs)<=freqbnd(2),:);
    TrespF = TrespF(:,round(freqs)>=freqbnd(1) & round(freqs)<=freqbnd(2),:);
    
    % Store data for grand-average, sensor-resolved TFR
    all_trl_avg(:,:,:,s) = cat(2,trl_avgL,trl_avg);
    all_trl_avg_full(:,:,:,s) = cat(2,trl_avg_fullL,trl_avg_full);
    gaB_respF(:,:,:,s) = cat(2,TrespFL,TrespF);
    
    all_trl_avg_tav(:,:,s) = mean(cat(2,trl_avgL(:,:,trltimes>=timewin(1) & trltimes<=timewin(2)),trl_avg(:,:,trltimes>=timewin(1) & trltimes<=timewin(2))),3);
    all_trl_avg_full_tav(:,:,s) = mean(cat(2,trl_avg_fullL(:,:,trltimes>=timewin(1) & trltimes<=timewin(2)),trl_avg_full(:,:,trltimes>=timewin(1) & trltimes<=timewin(2))),3);
    gaB_respF_tav(:,:,s) = mean(cat(2,TrespFL(:,:,trltimes>=timewin(1) & trltimes<=timewin(2)),TrespF(:,:,trltimes>=timewin(1) & trltimes<=timewin(2))),3);
    
    allfreqs = allfreqs(allfreqs>=freqbnd(1) & allfreqs<=freqbnd(2));
    
    % Create structure for this subject to mimic output from FT_FREQANALYSIS
    cstruct = struct;
    cstruct.label = grad.label;  % arbitrary label for motor lateralization 'channel'
    cstruct.fsample = 20;  % sampling rate
    cstruct.freq = allfreqs;
    cstruct.time = trltimes;
    cstruct.dimord = 'chan_freq_time';
    
    cstruct_tav = cstruct;
    cstruct_tav.time = trltimes(trltimes>=timewin(1) & trltimes<=timewin(2));
    
    % Add structures to group-level arrays
    cstruct.powspctrm = cat(2,trl_avgL,trl_avg); all_av{s} = cstruct;
    cstruct.powspctrm = cat(2,trl_avg_fullL,trl_avg_full); all_av_full{s} = cstruct;
    cstruct.powspctrm = cat(2,TrespFL,TrespF); allB_respF{s} = cstruct;
    cstruct.powspctrm = zeros(size(cstruct.powspctrm));
    allBnullF{s} = cstruct;
    
    cstruct_tav.powspctrm = cat(2,trl_avgL(:,:,trltimes>=timewin(1) & trltimes<=timewin(2)),trl_avg(:,:,trltimes>=timewin(1) & trltimes<=timewin(2))); all_av_tav{s} = cstruct_tav;
    cstruct_tav.powspctrm = cat(2,trl_avg_fullL(:,:,trltimes>=timewin(1) & trltimes<=timewin(2)),trl_avg_full(:,:,trltimes>=timewin(1) & trltimes<=timewin(2))); all_av_full_tav{s} = cstruct_tav;
    cstruct_tav.powspctrm = cat(2,TrespFL(:,:,trltimes>=timewin(1) & trltimes<=timewin(2)),TrespF(:,:,trltimes>=timewin(1) & trltimes<=timewin(2))); allB_respF_tav{s} = cstruct_tav;
    cstruct.powspctrm = zeros(size(cstruct.powspctrm));
    allBnullF_tav{s} = cstruct;
end
fprintf('done.\n')

% Trim gradiometer struct & compute neighbours
cfg=[]; cfg.grad = grad; cfg.method = 'triangulation';
neighbours = ft_prepare_neighbours(cfg, all_av{1});

% % Run cluster-based permutation tests
% cfg_stat = [];
% cfg_stat.latency     = 'all';
% cfg_stat.frequency   = 'all';
% cfg_stat.avgoverchan = 'no';
% cfg_stat.avgovertime = 'no';
% cfg_stat.avgoverfreq = 'no';
% cfg_stat.parameter   = 'powspctrm';
% cfg_stat.method      = 'montecarlo';
% cfg_stat.statistic   = 'depsamplesT';  % also possibly try 'depsamplesregrT'
% cfg_stat.alpha       = 0.05;
% cfg_stat.correctm    = 'cluster';
% % cfg_stat.correctm    = 'no';
% cfg_stat.clusteralpha = clustalpha;  % intial threhsold for forming clusters
% cfg_stat.clusterstatistic = 'maxsum';  % method for quantifying combined cluster test statistics
% cfg_stat.correcttail = 'prob';
% cfg_stat.numrandomization = nperms;
% cfg_stat.neighbours = neighbours;
% 
% Nsub = length(allsubj);
% cfg_stat.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg_stat.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% cfg_stat.ivar                = 1; % the 1st row in cfg.design contains the independent variable
% cfg_stat.uvar                = 2; % the 2nd row in cfg.design contains the subject number
% 
% stat1 = ft_freqstatistics(cfg_stat,all_av{:},allBnullF{:});
% stat2 = ft_freqstatistics(cfg_stat,allB_respF{:},allBnullF{:});
% 
% % Add null channels into grand-average matrices & stat masks
% temp_ga = zeros(length(fullgrad.chantype),size(all_trl_avg,2),size(all_trl_avg,3),size(all_trl_avg,4));
% temp_ga(ismember(fullgrad.label,grad.label),:,:,:) = all_trl_avg; all_trl_avg=temp_ga;
% temp_ga(ismember(fullgrad.label,grad.label),:,:,:) = gaB_respF; gaB_respF=temp_ga; clear temp_ga
% 
% tempmask = zeros(size(all_trl_avg,1),size(all_trl_avg,2),size(all_trl_avg,3));
% tempmask(ismember(fullgrad.label,grad.label),:,:) = stat1.mask; stat1.mask = tempmask;
% tempmask(ismember(fullgrad.label,grad.label),:,:) = stat2.mask; stat2.mask = tempmask; clear tempmask
% 
% % Plot grand-average, sensor-resolved raw TFRs
% TFR = struct;
% TFR.label = fullgrad.label;  % arbitrary label for motor lateralization 'channel'
% TFR.freq = allfreqs;
% TFR.time = stat1.time;
% TFR.dimord = 'chan_freq_time';
% TFR.grad = fullgrad;
% TFR.cfg = cfg;
% TFR.mask = stat1.mask;
% TFR.powspctrm = squeeze(nanmean(all_trl_avg,4));
% 
% cfg = [];
% cfg.colormap       = cmap;
% cfg.maskparameter  = 'mask';
% cfg.maskstyle      = 'outline';
% cfg.showlabels     = 'yes';	
% cfg.layout         = fullgrad;
% cfg.zlim           = [-1.0 1.0];
% 
% figure 
% ft_multiplotTFR(cfg, TFR);
% 
% % Plot grand-average, sensor-resolved response encoding
% TFR.mask = stat2.mask;
% TFR.powspctrm = squeeze(mean(gaB_respF,4));
% 
% cfg = [];
% cfg.colormap       = cmap;
% cfg.maskparameter  = 'mask';
% cfg.maskstyle      = 'outline';
% cfg.showlabels     = 'no';	
% cfg.layout         = grad;
% cfg.zlim           = [-1.3 1.3];
% 
% figure 
% ft_multiplotTFR(cfg, TFR);



% Run cluster-based permutation tests - ONLY CHAN*FREQ
cfg_stat = [];
cfg_stat.latency     = 'all';
cfg_stat.frequency   = 'all';
cfg_stat.avgoverchan = 'no';
cfg_stat.avgovertime = 'yes';
cfg_stat.avgoverfreq = 'no';
cfg_stat.parameter   = 'powspctrm';
cfg_stat.method      = 'montecarlo';
cfg_stat.statistic   = 'depsamplesT';  % also possibly try 'depsamplesregrT'
cfg_stat.alpha       = 0.05;
cfg_stat.correctm    = 'cluster';
% cfg_stat.correctm    = 'no';
cfg_stat.clusteralpha = clustalpha;  % intial threhsold for forming clusters
cfg_stat.clusterstatistic = 'maxsum';  % method for quantifying combined cluster test statistics
cfg_stat.correcttail = 'prob';
cfg_stat.numrandomization = nperms;
cfg_stat.neighbours = neighbours;

Nsub = length(allsubj);
cfg_stat.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg_stat.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg_stat.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg_stat.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat11 = ft_freqstatistics(cfg_stat,all_av_tav{:},allBnullF_tav{:});
stat12 = ft_freqstatistics(cfg_stat,allB_respF_tav{:},allBnullF_tav{:});


% Plot TF plot and topography weighted by contribution to cluster
w_sens=[]; TF=[];
for s = 1:size(stat12.mask,1)
    w_sens(s) = sum(stat12.stat(s,stat12.posclusterslabelmat(s,:)==1));
end
w_sens(isnan(w_sens)) = 0;
w_sens = w_sens./sum(w_sens);  % weights for averaging acrss sensors
for s = 1:size(stat12.mask,1)
    TF(s,:,:) = squeeze(mean(all_trl_avg(s,:,:,:),4)).*w_sens(s);
end
TF = squeeze(sum(TF,1));  % weighted average

cmap = colormapsmoothPM(500);

figure,
subplot(1,2,1), hold on
S = contourf(trltimes,stat12.freq,TF,100,'linestyle','none'); caxis([-max(max(abs(TF))) max(max(abs(TF)))].*0.9), colormap(cmap),
line([0 0],[min(stat12.freq) max(stat12.freq)],'LineStyle','-','Color',[0 0 0],'LineWidth',1)
line([1.3 1.3],[min(stat12.freq) max(stat12.freq)],'LineStyle','--','Color',[0 0 0],'LineWidth',1)
% cb=colorbar('North'); ylabel(cb,'\Sigma T-score','FontSize',15), set(cb,'TickLen',[0.04 0.04])
set(gca,'FontName','Arial','LineWidth',1,'FontSize',8,'TickDir','out','box','off','XTick',[0 0.4 0.8 1.3],'XTickLabel',{'Cue' '0.4' '0.8' 'Go'},'YTick',[3 7 20 50 100],'yscale','log')
xlim([-0.1 max(trltimes)]); ylim([min(stat12.freq) max(stat12.freq)])
ylabel('Frequency (Hz)','FontSize',10), xlabel('Time (s)','FontSize',10)


w_freq=[]; topo=[];
for f = 1:size(stat12.mask,2)
    w_freq(f) = sum(stat12.stat(stat12.posclusterslabelmat(:,f)==1,f));
end
w_freq(isnan(w_freq)) = 0;
w_freq = w_freq./sum(w_freq);  % weights for averaging acrss frequencies
for f = 1:size(stat12.mask,2)
    topo(:,f) = squeeze(mean(all_trl_avg_tav(:,f,:),3)).*w_freq(f);
end
topo = squeeze(sum(topo,2));  % weighted average

tempmask = zeros(length(fullgrad.chantype),1);
tempmask(ismember(fullgrad.label,grad.label)) = topo; topo = tempmask;
clear tempmask

cfgtopo                     = [];
cfgtopo.comment             = 'no';
cfgtopo.marker              = 'off';
cfgtopo.highlight           = 'off';
cfgtopo.colormap            = colormapsmoothPM(500);
cfgtopo.renderer            = 'painters';
cfgtopo.style               = 'straight';
cfgtopo.zlim                = 'maxabs';
% cfgtopo.gridscale           = 200;

topodata.dimord     = 'chan_time';
topodata.time       = 0;
topodata.label      = fullgrad.label;
topodata.grad       = fullgrad;
topodata.avg        = topo;

subplot(1,2,2), hold on
ft_topoplotER(cfgtopo, topodata);


% Plot different topos
tempmask = zeros(length(fullgrad.chantype),1);
tempmask(ismember(fullgrad.label,grad.label),:,:) = w_sens; w_sens_topo = tempmask;
tempmask(ismember(fullgrad.label,grad.label)) = squeeze(mean(mean(all_trl_avg_tav(:,allfreqs>=8 & allfreqs<=30,:),2),3)); LItopo = tempmask;
tempmask(ismember(fullgrad.label,grad.label),:,:) = squeeze(mean(mean(gaB_respF_tav(:,allfreqs>=8 & allfreqs<=30,:),2),3)); Ttopo = tempmask;
clear tempmask

figure,
subplot(1,5,1), hold on

ft_topoplotER(cfgtopo, topodata); title('Weighted LI')
subplot(1,5,2), hold on
topodata.avg        = w_sens_topo;
ft_topoplotER(cfgtopo, topodata); title('Weights')
subplot(1,5,3), hold on
topodata.avg        = LItopo;
ft_topoplotER(cfgtopo, topodata); title('LI: 8-30Hz')
subplot(1,5,4), hold on
topodata.avg        = Ttopo;
ft_topoplotER(cfgtopo, topodata); title('T-score: 8-30Hz')
subplot(1,5,5), hold on
topodata.avg        = all_trl_avg_full_tav(:,allfreqs>=8 & allfreqs<=30,:);
ft_topoplotER(cfgtopo, topodata); title('Raw: 8-30Hz')


% Plot single time course of motor localizer
w_sens_freq = stat12.stat;
w_sens_freq(stat12.posclusterslabelmat~=1) = 0;
w_sens_freq = w_sens_freq./sum(sum(w_sens_freq));

ml=[];
for subj = 1:size(all_trl_avg,4)
    for t = 1:size(all_trl_avg,3)
        ml(subj,t) = sum(sum(squeeze(all_trl_avg(:,:,t,subj)).*w_sens_freq));
    end
end

figure, hold on
shadedErrorBar(trltimes,mean(ml,1),std(ml,[],1)./sqrt(size(ml,1)),{'Color','k'},0)
plot([min(trltimes) max(trltimes)],[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k')
set(gca,'FontName','Arial','LineWidth',1,'FontSize',8,'TickDir','out','box','off','XTick',[0 0.4 0.8 1.3],'XTickLabel',{'Cue' '0.4' '0.8' 'Go'})
xlim([-0.1 max(trltimes)]);
ylabel('Lateralization (dB)','FontSize',10), xlabel('Time (s)','FontSize',10)

% Save
save([loadpath,'Clust_corrected_LIhalf_AllFreqs.mat'],'stat*','ga*','all_*','allfreqs','neighbours','fullgrad','grad','nperms','trltimes','timewin','cmap','clustalpha','w_*','ml')
save([loadpath,'Motor_loc_weights.mat'],'allfreqs','fullgrad','grad','timewin','w_*')



