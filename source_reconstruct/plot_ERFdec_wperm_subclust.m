clear, close all

% Paths dervied from processing options
allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};  % ORDER IS IMPORTANT FOR TRIAL-IDs!!!!

addpath '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts';
addpath '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/Gen_fun';
addpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts/FMINSEARCHBND')
addpath('/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/agg/output/av')

megpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/decodeERF/SSbase_40Hz/'];

% participants to include
include = 1:17;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dorsal stream + motor, low-freq %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clusters = {{'vfcPrimary'};...
            {'vfcEarly';'vfcV3ab';'vfcIPS01';'vfcIPS23'};...
            {'JWG_IPS_PCeS';'HCPMMP1_premotor';'JWG_M1'}};  % JW/Glasser

cnames = {'V1','Extrastriate','Motor'};
tcols = [0 0.35 0; 0 0.825 0.4; 0 0 0.9];

% load data
gaLLR={};
for c = 1:length(clusters)
    gaLLR{c}=[];
    for subj = 1:length(allsubj)
        fprintf('Subj %d...\n',subj)
        cLLR=[];
        for cc = 1:length(clusters{c})
            % load decoding output
            fname = [megpath,allsubj{subj},'_',clusters{c}{cc},'_full_finegrainERF.csv'];
            fprintf('Loading %s...\n',fname)
            t = readtable(fname);
            
            % initialize maximum times vector (some samples may have different latency ranges)
            if ~exist('ts')
                ts = min(round(t.latency,3)):0.005:max(round(t.latency,3));
            end
            
            % store decoding scores
            sLLR=nan(length(ts),length(unique(t.sample)));
            for smp = unique(t.sample)'
                ctimes = round(t.latency(strcmp(t.target,'LLR') & t.sample==smp),3);
                sLLR(find(abs(ts-min(ctimes))==min(abs(ts-min(ctimes)))):find(abs(ts-(max(ctimes)))==min(abs(ts-max(ctimes)))),smp) = t.test_correlation(strcmp(t.target,'LLR') & t.sample==smp);
            end
            
            % average over sample positions
            cLLR = [cLLR; nanmean(sLLR,2)'];
        end
        gaLLR{c} = [gaLLR{c} mean(cLLR,1)'];
    end
end

latwin = [0.05 1]; % time window (s) around sample onset to search for half-max latency

% plot
axlw=0.5; fs=7;
fig_w = 8; fig_h = 10;

h =  findobj('type','figure');
cfig = length(h)+1;
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

subplot(2,1,1), hold on
plot([min(ts) max(ts)],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--')
for c = 1:length(clusters)
    if ~isempty(gaLLR{c})
        shadedErrorBar(ts,mean(gaLLR{c},2),std(gaLLR{c},[],2)./sqrt(size(gaLLR{c},2)),{'Color',tcols(c,:)},0)
    end
end
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off')
xlim([min(ts) max(ts)]), ylim([-0.05 0.45])
xlabel('Time rel. to sample onset (s)'), ylabel('T-score'), title('LLR')

subplot(2,1,2), hold on
plot([min(ts) max(ts)],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--')
for c = 1:length(clusters)
    if ~isempty(gaLLR{c})
        plot(ts,mean(gaLLR{c},2)./max(mean(gaLLR{c},2)),'Color',tcols(c,:))
    end
end
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off')
xlim([min(ts) max(ts)]), ylim([-0.15 1.05])
xlabel('Time rel. to sample onset (s)'), ylabel('Normalized T-score')


% compute latency @ half max and plot
newsmpts = min(ts):0.001:max(ts);
interp_method = 'spline';
tau_plus = 0.05;  % time (s) to go beyond peak latency when picking starting point for timescale estimation
latLLR=[]; fitLLR=[];
for c = 1:length(clusters)
    for subj = 1:size(gaLLR{c},2)
        m = interp1(ts,gaLLR{c}(:,subj)',newsmpts,interp_method);   % interpolate to ms resolution for fine-grained latency estimate
        m = m(newsmpts>=latwin(1) & newsmpts<=latwin(2)); rts=newsmpts(newsmpts>=latwin(1) & newsmpts<=latwin(2));
        m = m./max(m);
        latLLR(subj,c) = rts(find(m>=0.5,1,'first'));   % latency of half-maximum
        
        peaklatLLR(subj,c) = rts(find(m==1,1,'first'));   % peak latency
        m = interp1(ts,gaLLR{c}(:,subj)',newsmpts,interp_method);
        m = m./m(find(newsmpts<=peaklatLLR(subj,c),1,'last'));
        normLLR(subj,c,:) = m;
        
        m = interp1(ts,gaLLR{c}(:,subj)',newsmpts,interp_method);   % interpolate to ms resolution
        m = m(newsmpts>=latwin(1)); rts = newsmpts(newsmpts>=latwin(1)); % throw out everything before start of latency window
        peakts = find(rts<=latwin(2));   % only search for peak latencies within first Xs of encoding function
        rts=rts(find(m(peakts)==max(m(peakts)),1,'first'):end); m=m(find(m(peakts)==max(m(peakts)),1,'first'):end); % throw out everything before peak
        m=m(rts>=(min(rts)+tau_plus)); rts=rts(rts>=(min(rts)+tau_plus));  % throw out everything before start of timescale fitting window
        m = m./m(1);  % normalize
        data_in = [rts' m'];
        fitLLR(subj,c) = fminsearchbnd(@(params) fit_exp(params,data_in), [0.15], [0.001], [3]);
        
        avEnc(subj,c) = mean(gaLLR{c}(:,subj));  % compute time-averaged encoding strength, used for weighted averaging
    end
end

% compute roi-wise weights for weighted averages
wSS = avEnc; wSS(wSS<0)=0;  % setting negative values to zero, meaning participants who don't show av effect won't contribute
wSS = wSS./repmat(sum(wSS,1),size(wSS,1),1);

% run pair-wise comparisons
nperm=10000;
wpLLRsig=nan(length(clusters)); wpLLRsig_fit=nan(length(clusters)); w4wp=[];
e = triu(LLRsig,-0);
disp('Running pairwise weighted permutation tests...')
for c1 = 1:length(clusters)
    for c2 = 1:length(clusters)
        if e(c1,c2)==0
            w = min([avEnc(:,c1) avEnc(:,c2)],[],2); w(w<0)=0; w=w./sum(w);  % compute subject weights for this specific pairwise test
            d1 = latLLR(:,c1)-latLLR(:,c2);
            d2 = log(fitLLR(:,c1))-log(fitLLR(:,c2));
            nsubj=length(w); pdist1=zeros(1,nperm); pdist2=pdist1; 
            for i=1:nperm
                sn = (rand(nsubj,1)>.5)*2-1;  % assign random sign to each subj
                pdist1(i) = (d1.*sn)'*w;
                pdist2(i) = (d2.*sn)'*w;
            end
            wpLLRsig(c1,c2) = length(find(abs(pdist1)>abs(d1'*w)))/nperm;  % two-tailed test
            wpLLRsig_fit(c1,c2) = length(find(abs(pdist2)>abs(d2'*w)))/nperm;
        end
    end
end

% calculate weighted means/confidence intervals for normalized LLR encoding trajectories
disp('Generating confidence intervals...')
nperm = 1000;
for c = 1:size(normLLR,2)
    wmLLR(c,:) = squeeze(normLLR(:,c,:))'*wSS(:,c);
    for i = 1:nperm
        smp = randsample(size(wSS,1),size(wSS,1),true,wSS(:,c));  % weighted sample
        iLLR(i,:) = squeeze(mean(normLLR(smp,c,:),1));
    end
    wciLLR(c,:,1:2) = prctile(iLLR,[2.5 97.5],1)';  % compute 95% CIs
end

% plot latencies/timescales with weighted averages/permutation tests
fig_w = 10.5; fig_h = 9.6;

h =  findobj('type','figure');
cfig = length(h)+1;
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

subplot(2,2,1), hold on
plot(1:size(latLLR,2),mean(latLLR,1),'Color',[0.3 0.3 0.3])
for cc = 1:size(latLLR,2)
    m = latLLR(:,cc)'*wSS(:,cc);
    s1=scatter(cc,m,40); set(s1,'MarkerFaceColor',tcols(cc,:),'MarkerEdgeColor',tcols(cc,:));  % data
end
xlim([0.5 length(clusters)+0.5]), set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[1:length(clusters)],'XTickLabel',cnames,'XTickLabelRotation',45)
ylim([0.22 0.43]), 
ylabel('Half-max latency (s)'), title('LLR')

subplot(2,2,2), hold on
plot(1:size(fitLLR,2),mean(fitLLR,1),'Color',[0.3 0.3 0.3])
for cc = 1:size(fitLLR,2)
    m = fitLLR(:,cc)'*wSS(:,cc);
    s1=scatter(cc,m,40); set(s1,'MarkerFaceColor',tcols(cc,:),'MarkerEdgeColor',tcols(cc,:));  % data
end
xlim([0.5 length(clusters)+0.5]), set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[1:length(clusters)],'XTickLabel',cnames,'XTickLabelRotation',45,'yscale','log')
% ylim([0.22 0.31]), 
ylabel('Tau (s)')

cmap = [[linspace(1,1,2)' linspace(1,1,2)' linspace(1,0,2)'];
        [linspace(1,1,200)' linspace(1,0,200)' linspace(0,0,200)']];

subplot(2,2,3), hold on
imagesc(1:length(clusters),1:length(clusters),-log(wpLLRsig(1:length(clusters),1:length(clusters))),-log([0.05 0.001])), xlim([0.5 length(clusters)+0.5]), ylim([0.5 length(clusters)+0.5]), colormap(cmap)
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','ydir','normal','XTick',[1:length(clusters)],'XTickLabel',cnames,'YTick',[1:length(clusters)],'YTickLabel',cnames,'XTickLabelRotation',45),
subplot(2,2,4), hold on
imagesc(1:length(clusters),1:length(clusters),-log(wpLLRsig_fit(1:length(clusters),1:length(clusters))),-log([0.05 0.001])), xlim([0.5 length(clusters)+0.5]), ylim([0.5 length(clusters)+0.5]), colormap(cmap)
set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','ydir','normal','XTick',[1:length(clusters)],'XTickLabel',cnames,'YTick',[1:length(clusters)],'YTickLabel',cnames,'XTickLabelRotation',45),
cb=colorbar('East'); ylabel(cb,'p')
set(cb,'Position',[0.86, 0.16, 0.014, 0.18],'Ticks',-log([0.05 0.01 0.001]),'TickLabels',[0.05 0.01 0.001])



% Make final figure
cnamesT = {'V1','Extrastriate','Motor'};

fig_w2 = 8.6; % figure width
fig_h2 = 6.4; % figure height

fs = 6.5; fs_label = 8;
lw = 1; axlw = 1;

jitbnd = [-0.3 0.3];
minscatsize = 2;
maxscatsize = 10;

nsubj = length(include);

h =  findobj('type','figure');
cfig = length(h)+1;
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w2 fig_h2],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

% plot trajectories
s1=subplot(1,3,1); hold on

plot([min(ts) max(ts)],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--')
for cc = 1:size(wciLLR)
    f=fill([newsmpts fliplr(newsmpts)],[wciLLR(cc,:,1) fliplr(wciLLR(cc,:,2))],tcols(cc,:)); set(f,'FaceColor',tcols(cc,:),'EdgeColor',tcols(cc,:));
end
for cc = 1:size(wciLLR)
    plot(newsmpts,wmLLR(cc,:),'Color',tcols(cc,:))
end

set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[0 0.4 0.8 1.2],'YTick',[0 0.4 0.8])
xlim([min(ts) max(ts)]), ylim([-0.14 1.0])
xlabel('Time from sample onset (s)'), ylabel('Normalized decoding precision (a.u.)')

% plot latencies
clims=[0.1 0.205];

s2=subplot(1,3,2); hold on
rnds1 = dist_scatters(latLLR(include,1),0.1); rnds2 = dist_scatters(latLLR(include,2),0.1); rnds3 = dist_scatters(latLLR(include,3),0.1);
%rand(nsubj,1); rnds2 = rand(nsubj,1); rnds3 = rand(nsubj,1);

for s=1:size(latLLR,1)
    plot([(1+diff(jitbnd).*rnds1(s)+jitbnd(1)) (2+diff(jitbnd).*rnds2(s)+jitbnd(1))],latLLR(include(s),1:2),'LineWidth',0.4,'Color',[0.75 0.75 0.75])
    plot([(2+diff(jitbnd).*rnds2(s)+jitbnd(1)) (3+diff(jitbnd).*rnds3(s)+jitbnd(1))],latLLR(include(s),2:3),'LineWidth',0.4,'Color',[0.75 0.75 0.75])
end

sizes=wSS(:,1)./max(wSS(:,1)).*(maxscatsize-minscatsize) + minscatsize;
S=scatter(ones(nsubj,1).*1 + diff(jitbnd).*rnds1 + jitbnd(1),latLLR(include,1),sizes,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.6 0.6 0.6])
plot([1 1]+jitbnd, [latLLR(include,1)'*wSS(:,1) latLLR(include,1)'*wSS(:,1)],'LineWidth',2,'Color',tcols(1,:))

sizes=wSS(:,2)./max(wSS(:,2)).*(maxscatsize-minscatsize) + minscatsize;
S=scatter(ones(nsubj,1).*2 + diff(jitbnd).*rnds2 + jitbnd(1),latLLR(include,2),sizes,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.6 0.6 0.6])
plot([2 2]+jitbnd, [latLLR(include,2)'*wSS(:,2) latLLR(include,2)'*wSS(:,2)],'LineWidth',2,'Color',tcols(2,:))

sizes=wSS(:,3)./max(wSS(:,3)).*(maxscatsize-minscatsize) + minscatsize;
S=scatter(ones(nsubj,1).*3 + diff(jitbnd).*rnds3 + jitbnd(1),latLLR(include,3),sizes,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.6 0.6 0.6])
plot([3 3]+jitbnd, [latLLR(include,3)'*wSS(:,3) latLLR(include,3)'*wSS(:,3)],'LineWidth',2,'Color',tcols(3,:))

p_i = [wpLLRsig(2,1) wpLLRsig(3,2) wpLLRsig(3,1)];
sigx = [1.5 2.5 2]; sigy = [0.19 0.195 0.202];

if p_i(1)<0.1, plot([1 2],[sigy(1) sigy(1)],'Color',[0.6 0.6 0.6]), end
if p_i(2)<0.1, plot([2 3],[sigy(2) sigy(2)],'Color',[0.6 0.6 0.6]), end
if p_i(3)<0.1, plot([1 3],[sigy(3) sigy(3)],'Color',[0.6 0.6 0.6]), end

for sig = find(p_i<0.001), text(sigx(sig),sigy(sig),'***','HorizontalAlignment','center','Color',[0.6 0.6 0.6]), end
for sig = find(p_i<0.01 & p_i>=0.001), text(sigx(sig),sigy(sig),'**','HorizontalAlignment','center','Color',[0.6 0.6 0.6]), end
for sig = find(p_i<0.05 & p_i>=0.01), text(sigx(sig),sigy(sig),'*','HorizontalAlignment','center','Color',[0.6 0.6 0.6]), end
for sig = find(p_i<0.1 & p_i>=0.05), text(sigx(sig),sigy(sig),'~','HorizontalAlignment','center','Color',[0.6 0.6 0.6]), end

set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[],'YTick',[0.11 0.15 0.19])
%set(gca,'XTickLabel',cnamesT,'XTickLabelRotation',45)
xlim([0.3 3.5]); ylim(clims)
ylabel('Latency (s)','FontSize',fs_label)

% plot timescales
clims=[0.03 4];
s3=subplot(1,3,3); hold on
rnds1 = dist_scatters(fitLLR(include,1),0.1); rnds2 = dist_scatters(fitLLR(include,2),0.1); rnds3 = dist_scatters(fitLLR(include,3),0.1);

for s=1:size(fitLLR,1)
    plot([(1+diff(jitbnd).*rnds1(s)+jitbnd(1)) (2+diff(jitbnd).*rnds2(s)+jitbnd(1))],fitLLR(include(s),1:2),'LineWidth',0.4,'Color',[0.75 0.75 0.75])
    plot([(2+diff(jitbnd).*rnds2(s)+jitbnd(1)) (3+diff(jitbnd).*rnds3(s)+jitbnd(1))],fitLLR(include(s),2:3),'LineWidth',0.4,'Color',[0.75 0.75 0.75])
end

sizes=wSS(:,1)./max(wSS(:,1)).*(maxscatsize-minscatsize) + minscatsize;
S=scatter(ones(nsubj,1).*1 + diff(jitbnd).*rnds1 + jitbnd(1),fitLLR(include,1),sizes,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.6 0.6 0.6])
plot([1 1]+jitbnd, [fitLLR(include,1)'*wSS(:,1) fitLLR(include,1)'*wSS(:,1)],'LineWidth',2,'Color',tcols(1,:))

sizes=wSS(:,2)./max(wSS(:,2)).*(maxscatsize-minscatsize) + minscatsize;
S=scatter(ones(nsubj,1).*2 + diff(jitbnd).*rnds2 + jitbnd(1),fitLLR(include,2),sizes,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.6 0.6 0.6])
plot([2 2]+jitbnd, [fitLLR(include,2)'*wSS(:,2) fitLLR(include,2)'*wSS(:,2)],'LineWidth',2,'Color',tcols(2,:))

sizes=wSS(:,2)./max(wSS(:,3)).*(maxscatsize-minscatsize) + minscatsize;
S=scatter(ones(nsubj,1).*3 + diff(jitbnd).*rnds3 + jitbnd(1),fitLLR(include,3),sizes,'k','o');
set(S,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.6 0.6 0.6])
plot([3 3]+jitbnd, [fitLLR(include,3)'*wSS(:,3) fitLLR(include,3)'*wSS(:,3)],'LineWidth',2,'Color',tcols(3,:))

p_i = [wpLLRsig_fit(2,1) wpLLRsig_fit(3,2) wpLLRsig_fit(3,1)];
sigx = [1.5 2.5 2]; sigy = [3.25 3.5 3.8];

if p_i(1)<0.1, plot([1 2],[sigy(1) sigy(1)],'Color',[0.6 0.6 0.6]), end
if p_i(2)<0.1, plot([2 3],[sigy(2) sigy(2)],'Color',[0.6 0.6 0.6]), end
if p_i(3)<0.1, plot([1 3],[sigy(3) sigy(3)],'Color',[0.6 0.6 0.6]), end

for sig = find(p_i<0.001), text(sigx(sig),sigy(sig),'***','HorizontalAlignment','center','Color',[0.6 0.6 0.6]), end
for sig = find(p_i<0.01 & p_i>=0.001), text(sigx(sig),sigy(sig),'**','HorizontalAlignment','center','Color',[0.6 0.6 0.6]), end
for sig = find(p_i<0.05 & p_i>=0.01), text(sigx(sig),sigy(sig),'*','HorizontalAlignment','center','Color',[0.6 0.6 0.6]), end
for sig = find(p_i<0.1 & p_i>=0.05), text(sigx(sig),sigy(sig),'~','HorizontalAlignment','center','Color',[0.6 0.6 0.6]), end

set(gca,'FontName','Arial','LineWidth',axlw,'FontSize',fs,'TickDir','out','box','off','XTick',[1 2 3],'YTick',[0.1 0.5 2.4],'yscale','log')
set(gca,'XTickLabel',cnamesT,'XTickLabelRotation',45)
xlim([0.3 3.5]); ylim(clims)
ylabel('tau (s)','FontSize',fs_label)

% position panels
w0 = 0.4; w1 = 0.22; wgap = 0.16; woff = 0.12;
h0 = 0.42; h1 = 0.37; hgap = 0.06; hoff = 0.15;

set(s1, 'Position', [woff, hoff+0.2, w0, h0])   % [left bottom width height]
set(s2, 'Position', [woff+(w0+wgap), hoff+h1+hgap, w1, h1])
set(s3, 'Position', [woff+(w0+wgap), hoff, w1, h1])   % [left bottom width height]



    
