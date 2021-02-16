clear, close all

% file details
allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};  % ORDER IS IMPORTANT FOR TRIAL-IDs!!!!

basetype = 'dB_common';  % options: 'dB_st', 'dB_common', 'pc_st' or 'pc_common'
modeltype = 'fitted_lin';  % switch b/w 'normative', 'fitted', 'fitted_np' & 'fitted_lin'
belief_type = 'lpr';  % switch b/w 'psi' and 'lpr'

% Paths dervied from processing options
if strcmp(modeltype,'normative'), str2 = 'output';
elseif strcmp(modeltype,'fitted'), str2 = 'output_fitted';
elseif strcmp(modeltype,'fitted_np'), str2 = 'output_fitted_np';
elseif strcmp(modeltype,'fitted_lin'), str2 = 'output_fitted_lin'; end
loadpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/TF/wMotorLoc/',basetype,filesep,str2,filesep];

% Loop through subjects/clusters
alldat=[]; allbelief=[]; allts=[];
fprintf('\nProcessing subject ')
for s = 1:length(allsubj)
    fprintf('%d,',s)
    
    load([loadpath,allsubj{s},'_CPcentered_output_appML_tonly.mat'])
    
    alldat(s,:,:) = dataCP;
    
    if strcmp(belief_type,'psi')
        allbelief(s,:,:) = psiCP;
    elseif strcmp(belief_type,'lpr')
        allbelief(s,:,:) = lprCP;
    end
    
    allts(s,:) = trlnums;
end
fprintf('done.\n')

% Specify samnple-wise granularity of plotting
bins = {[1:4],[5:8],[9:12]};

% Plot belief time-courses
ccols = [[ones(2,1) linspace(0.6,0.05,2)' linspace(0.6,0.05,2)'];...
         [0.5 0 0]];
     
figure, hold on
plot([0.5 12.5],[0 0],'LineStyle','--','Color',[0.6 0.6 0.6])
for cp = 1:length(bins)
    plot(1:12,squeeze(mean(mean(allbelief(:,bins{cp},:),2),1)),'Color',ccols(cp,:))
end
xlabel('Sample position'), ylabel('\psi'), set(gca,'TickDir','out','box','off','XTick',[2:2:12])

% Plot clusters of LI trajectories overlaid on belief trajectories
belief_scaling = 0.045;
belief_offset = 0.3;
dotsize = 20;
xlims = [0 5.4];

figure, hold on
plot([min(trltimes) max(trltimes)],[0 0],'LineStyle','--','Color',[0.6 0.6 0.6])
for cp = 1:length(bins)
    shadedErrorBar(trltimes,squeeze(mean(mean(alldat(:,bins{cp},:),2),1)),squeeze(std(mean(alldat(:,bins{cp},:),2),[],1))./sqrt(size(alldat,1)),{'Color',ccols(cp,:),'LineWidth',1.5},0)
end
for cp = 1:length(bins)
    plot(((1:12).*0.4)+belief_offset,squeeze(mean(mean(allbelief(:,bins{cp},:),2),1)).*belief_scaling,'Color',ccols(cp,:))
    s=scatter(((1:12).*0.4)+belief_offset,squeeze(mean(mean(allbelief(:,bins{cp},:),2),1)).*belief_scaling,dotsize); set(s,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',ccols(cp,:))
end
xlim(xlims), xlabel('Time (s)'), ylabel('LI'), title('Motor'), set(gca,'TickDir','out','box','off','XTick',[0:5])

