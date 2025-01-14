

clear
close all

%% Load Data
load Data.mat

fc_abs = abs(fc); 

%% Connections to Preserve
NetLbls = {'FC_S_C','FC_D_C','FC_G_C','FC_R_C','FC_G_C_d','FC_R_C_d','All'};

dn1 = [0.5:-.05:0.05];

for i = 1:numel(dn1)
    
%     i = 1
    dn = dn1(i)

    [perc_overlap_all_bp(i,:),bp_overlap(:,i)] = percent_overlap(fc,sc,ec,gc,rc,gc_dist,rc_dist,dn)

end

%% Figures
close all
figure('Units','centimeters','Position',[10 10 10 5])
hold on
bar(perc_overlap_all_bp)
ylim([0 100])
% xlim([0.25 1.55])
set(gca,'XTick',1:numel(dn1),'XtickLabel',100*dn1, ...
            'FontSize',8,'XColor','k','YColor','k')

xlabel('% Connections Preserved','FontSize',10)
ylabel('% Shared Connections')
f2s = ['HCP_Corr_Figures/HCP_absFC_Shared_Conn_BioPhys.eps'];
saveas(gcf,f2s,'epsc')


% close all
figure('Units','centimeters','Position',[10 10 10 5])
imagesc(bp_overlap)
colorbar
colormap(flipud(hot))
caxis([0 65])
set(gca,'XTick',1:numel(dn1),'XtickLabel',100*dn1, ...
            'FontSize',8,'XColor','k','YColor','k')
set(gca,'YTick',1:numel(NetLbls),'YtickLabel',NetLbls, ...
            'FontSize',8,'XColor','k','YColor','k')
xlabel('% Connections Preserved','FontSize',10)
% f2s = ['HCP_Corr_Figures/HCP_absFC_RFN_Overlap.eps'];
% saveas(gcf,f2s,'epsc')
