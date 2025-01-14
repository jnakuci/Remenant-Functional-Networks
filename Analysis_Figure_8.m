%% LA5c Clinical RFNs Differences

clear
close all

%% Toolboxes
addpath ~/Dropbox/Toolbox/plotSpread/

%%
grps = {'HC','SCZ','BP', 'ADHD'};
Net = {'dmri','dist','gene','recep','gene_dist','recep_dist'};
NetLbls = {'FC_S_C','FC_D_C','FC_G_C','FC_R_C','FC_G_C_d','FC_R_C_d'};

colors = lines(numel(NetLbls));
nodeId = Atlas_Info.nodeId;
Nroi = 200;

mk = {'.','.','s','s','d','d'};
mks = [20, 20, 8, 8, 8, 8];


%% Colors 
colors = lines(1);
colors(2,:) = [0.9290 0.6940 0.1250];
colors(3,:) = [ 0 0.5 0];
colors(4,:) = [ 0.65 0 0];
colors(5,:) = [ 0 0.75 0];
colors(6,:) = [ 0.9 0 0];


%% Load Data
load Data.mat

%% Create RFN

dn = 0.16;

% HC
for s = 1:size(LA5c_data.HC,3)
    hc_full(s,:) = unwrapMat(LA5c_data.HC(:,:,s));
    Masked.HC.sc_masked(s,:) = unwrapMat(create_RFN(LA5c_data.HC(:,:,s), sc, dn));
    Masked.HC.ec_masked(s,:) = unwrapMat(create_RFN(LA5c_data.HC(:,:,s), ec, dn));
    Masked.HC.gc_masked(s,:) = unwrapMat(create_RFN(LA5c_data.HC(:,:,s), gc, dn));
    Masked.HC.gcdist_masked(s,:) = unwrapMat(create_RFN(LA5c_data.HC(:,:,s), gc_dist, dn));
    Masked.HC.rc_masked(s,:) = unwrapMat(create_RFN(LA5c_data.HC(:,:,s), rc, dn));
    Masked.HC.rcdist_masked(s,:) = unwrapMat(create_RFN(LA5c_data.HC(:,:,s), rc_dist, dn));
end

% SCZ 
for s = 1:size(LA5c_data.SCZ,3)
    scz_full(s,:) = unwrapMat(LA5c_data.SCZ(:,:,s));
    Masked.SCZ.sc_masked(s,:)  = unwrapMat(create_RFN(LA5c_data.SCZ(:,:,s), sc, dn));
    Masked.SCZ.ec_masked(s,:)  = unwrapMat(create_RFN(LA5c_data.SCZ(:,:,s), ec, dn));
    Masked.SCZ.gc_masked(s,:)  = unwrapMat(create_RFN(LA5c_data.SCZ(:,:,s), gc, dn));
    Masked.SCZ.gcdist_masked(s,:)  = unwrapMat(create_RFN(LA5c_data.SCZ(:,:,s), gc_dist, dn));
    Masked.SCZ.rc_masked(s,:)  = unwrapMat(create_RFN(LA5c_data.SCZ(:,:,s), rc, dn));
    Masked.SCZ.rcdist_masked(s,:)  = unwrapMat(create_RFN(LA5c_data.SCZ(:,:,s), rc_dist, dn));
end

% BP 
for s = 1:size(LA5c_data.BP,3)
    bp_full(s,:) = unwrapMat(LA5c_data.BP(:,:,s));
    Masked.BP.sc_masked(s,:) = unwrapMat(create_RFN(LA5c_data.BP(:,:,s), sc, dn));
    Masked.BP.ec_masked(s,:) = unwrapMat(create_RFN(LA5c_data.BP(:,:,s), ec, dn));
    Masked.BP.gc_masked(s,:) = unwrapMat(create_RFN(LA5c_data.BP(:,:,s), gc, dn));
    Masked.BP.gcdist_masked(s,:) = unwrapMat(create_RFN(LA5c_data.BP(:,:,s), gc_dist, dn));
    Masked.BP.rc_masked(s,:) = unwrapMat(create_RFN(LA5c_data.BP(:,:,s), rc, dn));
    Masked.BP.rcdist_masked(s,:) = unwrapMat(create_RFN(LA5c_data.BP(:,:,s), rc_dist, dn));
end

% ADHD
for s = 1:size(LA5c_data.ADHD,3)
    adhd_full(s,:) = unwrapMat(LA5c_data.ADHD(:,:,s));
    Masked.ADHD.sc_masked(s,:) = unwrapMat(create_RFN(LA5c_data.ADHD(:,:,s), sc, dn));
    Masked.ADHD.ec_masked(s,:) = unwrapMat(create_RFN(LA5c_data.ADHD(:,:,s), ec, dn));
    Masked.ADHD.gc_masked(s,:) = unwrapMat(create_RFN(LA5c_data.ADHD(:,:,s), gc, dn));
    Masked.ADHD.gcdist_masked(s,:) = unwrapMat(create_RFN(LA5c_data.ADHD(:,:,s), gc_dist, dn));
    Masked.ADHD.rc_masked(s,:) = unwrapMat(create_RFN(LA5c_data.ADHD(:,:,s), rc, dn));
    Masked.ADHD.rcdist_masked(s,:) = unwrapMat(create_RFN(LA5c_data.ADHD(:,:,s), rc_dist, dn));
end


%% RFNs Estimated Differences (T-value) 

% SCZ
[~,~,~,stats] = ttest2(scz_full,hc_full);
tval.scz.full = stats.tstat;
scz_tfull = nanmean(stats.tstat);

[~,~,~,stats] = ttest2(Masked.SCZ.sc_masked, Masked.HC.sc_masked);
tval.scz.dmri = stats.tstat;
Tvals.dmri.SCZ = 100.*(nanmean(stats.tstat) - scz_tfull )./scz_tfull;

[~,~,~,stats] = ttest2(Masked.SCZ.ec_masked,Masked.HC.ec_masked);
tval.scz.dist = stats.tstat;
Tvals.dist.SCZ = 100.*(nanmean(stats.tstat) - scz_tfull )./scz_tfull;

[~,~,~,stats] = ttest2(Masked.SCZ.gc_masked, Masked.HC.gc_masked);
tval.scz.gene = stats.tstat;
Tvals.gene.SCZ = 100.*(nanmean(stats.tstat) - scz_tfull )./scz_tfull;

[~,~,~,stats] = ttest2(Masked.SCZ.gcdist_masked, Masked.HC.gcdist_masked);
tval.scz.gene_dist = stats.tstat;
Tvals.gene_dist.SCZ = 100.*(nanmean(stats.tstat) - scz_tfull )./scz_tfull;

[~,~,~,stats] = ttest2(Masked.SCZ.rc_masked,Masked.HC.rc_masked);
tval.scz.recep = stats.tstat;
Tvals.recep.SCZ = 100.*(nanmean(stats.tstat) - scz_tfull )./scz_tfull;

[~,~,~,stats] = ttest2(Masked.SCZ.rcdist_masked,Masked.HC.rcdist_masked);
tval.scz.recep_dist = stats.tstat;
Tvals.recep_dist.SCZ = 100.*(nanmean(stats.tstat) - scz_tfull )./scz_tfull;

clear stats 

% BP
[~,~,~,stats] = ttest2(bp_full,hc_full);
tval.bp.full = stats.tstat;
bp_tfull = nanmean(stats.tstat);

[~,~,~,stats] = ttest2(Masked.BP.sc_masked,Masked.HC.sc_masked);
tval.bp.dmri = stats.tstat;
Tvals.dmri.BP = 100.*(nanmean(stats.tstat) - bp_tfull)./bp_tfull;

[~,~,~,stats] = ttest2(Masked.BP.ec_masked,Masked.HC.ec_masked);
tval.bp.dist = stats.tstat;
Tvals.dist.BP = 100.*(nanmean(stats.tstat) - bp_tfull)./bp_tfull;

[~,~,~,stats] = ttest2(Masked.BP.gc_masked,Masked.HC.gc_masked);
tval.bp.gene = stats.tstat;
Tvals.gene.BP = 100.*(nanmean(stats.tstat) - bp_tfull)./bp_tfull;

[~,~,~,stats] = ttest2(Masked.BP.gcdist_masked,Masked.HC.gcdist_masked);
tval.bp.gene_dist = stats.tstat;
Tvals.gene_dist.BP = 100.*(nanmean(stats.tstat) - bp_tfull)./bp_tfull;

[~,~,~,stats] = ttest2(Masked.BP.rc_masked,Masked.HC.rc_masked);
tval.bp.recep = stats.tstat;
Tvals.recep.BP = 100.*(nanmean(stats.tstat) - bp_tfull)./bp_tfull;

[~,~,~,stats] = ttest2(Masked.BP.rcdist_masked,Masked.HC.rcdist_masked);
tval.bp.recep_dist = stats.tstat;
Tvals.recep_dist.BP = 100.*(nanmean(stats.tstat) - bp_tfull)./bp_tfull;

clear stats

% ADHD
[~,~,~,stats] = ttest2(adhd_full,hc_full);
tval.adhd.full = stats.tstat;
adhd_tfull = nanmean(stats.tstat);

[~,~,~,stats] = ttest2(Masked.ADHD.sc_masked,Masked.HC.sc_masked);
tval.adhd.dmri = stats.tstat;
Tvals.dmri.ADHD = 100.*(nanmean(stats.tstat) - adhd_tfull)./adhd_tfull;

[~,~,~,stats] = ttest2(Masked.ADHD.ec_masked,Masked.HC.ec_masked);
tval.adhd.dist = stats.tstat;
Tvals.dist.ADHD = 100.*(nanmean(stats.tstat) - adhd_tfull)./adhd_tfull;

[~,~,~,stats] = ttest2(Masked.ADHD.gc_masked,Masked.HC.gc_masked);
tval.adhd.gene = stats.tstat;
Tvals.gene.ADHD = 100.*(nanmean(stats.tstat) - adhd_tfull)./adhd_tfull;

[~,~,~,stats] = ttest2(Masked.ADHD.gcdist_masked,Masked.HC.gcdist_masked);
tval.adhd.gene_dist = stats.tstat;
Tvals.gene_dist.ADHD = 100.*(nanmean(stats.tstat) - adhd_tfull)./adhd_tfull;

[~,~,~,stats] = ttest2(Masked.ADHD.rc_masked,Masked.HC.rc_masked);
tval.adhd.recep = stats.tstat;
Tvals.recep.ADHD = 100.*(nanmean(stats.tstat) - adhd_tfull)./adhd_tfull;

[~,~,~,stats] = ttest2(Masked.ADHD.rcdist_masked,Masked.HC.rcdist_masked);
tval.adhd.recep_dist = stats.tstat;
Tvals.recep_dist.ADHD = 100.*(nanmean(stats.tstat) - adhd_tfull)./adhd_tfull;

clear stats

%% Average T-value

% SCZ
nanmean(tval.scz.full);
nanstd(tval.scz.full)./sqrt(numel(tval.scz.full))

% BP
nanmean(tval.bp.full)
nanstd(tval.bp.full)./sqrt(numel(tval.bp.full))

% ADHD
nanmean(tval.adhd.full)
nanstd(tval.adhd.full)./sqrt(numel(tval.adhd.full))


%% Random Network

for n = 1:1000

    n

    % Create Random Network
    nc = rand(Nroi); 
    nc = triu(nc,1) + triu(nc,1)';

    % Apply Rand Mask
    for s = 1:size(LA5c_data.HC,3)
        HC_rand(s,:) = unwrapMat(create_RFN(LA5c_data.HC(:,:,s),nc,dn));
    end
    
    for s = 1:size(LA5c_data.SCZ,3)
        SCZ_rand(s,:) = unwrapMat(create_RFN(LA5c_data.SCZ(:,:,s),nc,dn));
    end

    for s = 1:size(LA5c_data.BP,3)
        BP_rand(s,:) = unwrapMat(create_RFN(LA5c_data.BP(:,:,s),nc,dn));
    end
    
    for s = 1:size(LA5c_data.ADHD,3)
        ADHD_rand(s,:) = unwrapMat(create_RFN(LA5c_data.ADHD(:,:,s),nc,dn));
    end

    % T-value
    [~,~,~,stats] = ttest2(SCZ_rand,HC_rand);
    Tvals.rand.SCZ(n) = 100.*(nanmean(stats.tstat) - scz_tfull)./scz_tfull;

     [~,~,~,stats] = ttest2(BP_rand,HC_rand);
    Tvals.rand.BP(n) = 100.*(nanmean(stats.tstat) - bp_tfull)./bp_tfull;

    [~,~,~,stats] = ttest2(ADHD_rand,HC_rand);
    Tvals.rand.ADHD(n) = 100.*(nanmean(stats.tstat) - adhd_tfull)./adhd_tfull;
    
end


%% Figures 

% SCZ
gr = repmat(Tvals.rand.SCZ',1,6);

close all
figure('Units','centimeters','Position',[10 10 7 4])
hold on
plotSpread(gr)

for i = 1:numel(Net)
    Net{i}
   plot(i,Tvals.(Net{i}).SCZ,'marker',mk(i),'Markersize',mks(i),'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:))
end

ylim([-18 18])
xlim([0.5 numel(Net)+0.5])

set(gca,'XTick',1:numel(NetLbls),'XtickLabel',NetLbls, ...
    'FontSize',8,'XColor','k','YColor','k')
ylabel('\Delta (%)','FontSize',10)
% title('SCZ - HC','FontSize',10)

set(gcf,'Renderer','painters')
% f2s = ['LA5c_Figures/SCZ-HCP_indirect_tval_diff.eps'];
% saveas(gcf,f2s,'epsc')


% BP
gr = repmat(Tvals.rand.BP',1,6);
close all
figure('Units','centimeters','Position',[10 10 7 4])
hold on
plotSpread(gr)

for i = 1:numel(Net)
    Net{i}
   plot(i,Tvals.(Net{i}).BP,'marker',mk(i),'Markersize',mks(i),'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:))
end

ylim([-7 7])

xlim([0.5 numel(Net)+0.5])

set(gca,'XTick',1:numel(NetLbls),'XtickLabel',NetLbls, ...
    'FontSize',8,'XColor','k','YColor','k')
ylabel('\Delta (%)','FontSize',10)
% title('BP - HC','FontSize',10)

set(gcf,'Renderer','painters')
% f2s = ['LA5c_Figures/BP-HCP_indirect_tval_diff.eps'];
% saveas(gcf,f2s,'epsc')

% ADHD 
gr = repmat(Tvals.rand.ADHD',1,6);
close all
figure('Units','centimeters','Position',[10 10 7 4])
hold on
plotSpread(gr)

for i = 1:numel(Net)
    Net{i}
   plot(i,Tvals.(Net{i}).ADHD,'marker',mk(i),'Markersize',mks(i),'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:))
end

ylim([-7 7])
xlim([0.5 numel(Net)+0.5])

set(gca,'XTick',1:numel(NetLbls),'XtickLabel',NetLbls, ...
    'FontSize',8,'XColor','k','YColor','k')
ylabel('\Delta (%)','FontSize',10)
% title('ADHD - HC','FontSize',10)

set(gcf,'Renderer','painters')
% f2s = ['LA5c_Figures/ADHD-HCP_indirect_tval_diff.eps'];
% saveas(gcf,f2s,'epsc')

