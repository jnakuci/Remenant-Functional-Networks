%% Analysis: Gradients
clear
close all

%% Toolboxes 
addpath(genpath('Code/')) 
addpath ~/Dropbox/Toolbox/BCT/2019_03_03_BCT/
addpath(genpath('~/Dropbox/Toolbox/BrainSpace/'))  
addpath(genpath('~/Dropbox/Toolbox/matlab_GIfTI/'))  
addpath ~/Dropbox/Toolbox/Code/
addpath ../Atlas_Stuff/
addpath ~/Dropbox/Toolbox/plotSpread/

%% Load Data
load Data.mat

%%
load Node_Community_Info.mat

%% 
Net = {'dmri','dist','gene','recep','gene_dist','recep_dist'};
Net1 = {'full','dmri','dist','gene','recep','gene_dist','recep_dist'};
NetLbls = {'SC','DC','GC','RC','GC_E_D','RC_E_D'};
NetLblsLeg = {'FC_F_u_l_l','FC_S_C','FC_D_C','FC_G_C','FC_R_C','FC_G_C_d','FC_R_C_d'};
NetLblsLeg1 = {'FC_S_C','FC_D_C','FC_G_C','FC_R_C','FC_G_C_d','FC_R_C_d'};

%% Colors 
colors = lines(1);
colors(2,:) = [0.9290 0.6940 0.1250];
colors(3,:) = [ 0 0.5 0];
colors(4,:) = [ 0.6 0 0];
colors(5,:) = [ 0 0.75 0];
colors(6,:) = [ 1 0 0];

mk = {'.','.','s','s','d','d'};
mks = [20, 20, 8, 8, 8, 8];


%% Create RFN
dn = 0.16; % density
 
fc_abs = abs(fc);
RFN.fc_sc_rfn = create_RFN(fc_abs, sc, dn);
RFN.fc_ec_rfn = create_RFN(fc_abs, ec,dn);
RFN.fc_gc_rfn = create_RFN(fc_abs, gc,dn);
RFN.fc_gcdist_rfn = create_RFN(fc_abs, gc_dist, dn);
RFN.fc_rc_rfn = create_RFN(fc_abs, rc, dn);
RFN.fc_rcdist_rfn = create_RFN(fc_abs, rc_dist, dn);


%% Create Random Masked Networks
Nrand = 1000;

for i = 1:Nrand
     RFN.rand_masked(:,:,i)= my_mask_fc_rand(fc_abs,dn);
end

%% Gradients
clear g
g{1} = fc;
g{2} = RFN.fc_sc_rfn;
g{3} = RFN.fc_ec_rfn;
g{4} = RFN.fc_gc_rfn;
g{5} = RFN.fc_gcdist_rfn;
g{6} = RFN.fc_rc_rfn;
g{7} = RFN.fc_rcdist_rfn;

for i = 1:Nrand
    g{end+1} = RFN.rand_masked(:,:,i);
end 

% Fit Gradienst
gm = GradientMaps('kernel','cs','approach','dm','alignment','pa','n_components',2);
gm = gm.fit(g);


%% First two Gradients
close all
for i = 1:numel(Net1)

    v = gm.gradients{i};
    figure('Units','centimeters','Position',[10 10 5 5])
    hold on
    for j = 1:200
        plot(v(j,1),v(j,2),'.','MarkerSize', 7,'Color',Atlas_Info.nodeColor_Orig(j,:))
    end
  
    set(gca,'XColor','k','YColor','k','FontSize',8)
    title(NetLblsLeg{i},'FontSize',10)

    axis square
    set(gcf, 'renderer', 'painters')
%     f2s = ['HCP_Corr_Figures/HCP_' NetLblsLeg{i} '_GrpLevel_Gradients_DM.eps'];
%     saveas(gcf,f2s,'epsc')

end 

%% Gradient Dissimilarity Figures
close all

for c = 1:2

    v = [];
    for i = 1:Nrand+7
        v = [v gm.gradients{i}(:,c)];
    end
   
    sim = 1-abs(corr(v(:,1),v(:,2:end)));

    r = repmat(sim(1,8:end)',1,6); 

    figure('Units','centimeters','Position',[10 10 10 7])
    hold on
    plotSpread(r)

    for i = 1:6
        plot(i,sim(1,i),'marker',mk(i),'Markersize',mks(i),'Color',colors(i,:),'MarkerFaceColor',colors(i,:))
    end

    xlim([0.5 numel(Net)+0.5])
    set(gca,'XTick',1:numel(NetLbls),'XTickLabel',NetLblsLeg1, ...
        'FontSize',10,'XColor','k','YColor','k')

    ylabel('\delta','FontSize',12)
    set(gcf,'Renderer','painters')
%     f2s = ['HCP_Corr_Figures/HCP_Gradient_Similarity_GrpLevel_Dim' num2str(c) '_neg_corr.eps'];
%     saveas(gcf,f2s,'epsc')

    clear r sim v

end


%% Variance Explained

close all
figure('Units','centimeters','Position',[10 10 7 5])
hold on
for i = 2:7
    lambdas = gm.lambda{1,i};
    if i == 1
         plot(100.*lambdas ./ sum(lambdas),'.-','Color','k','MarkerSize',15);
    else
        plot(100.*lambdas ./ sum(lambdas),'.-','Color',colors(i-1,:),'MarkerSize',15);
    end
end

% ylim([0 35])

xlim([0.5 7])
% legend(NetLblsLeg,'FontSize',6)
set(gca,'XTick',1:7)
xlabel('Gradient','FontSize',10);

ylabel('Variance Explained (%)','FontSize',10);

f2s = ['HCP_Corr_Figures/Gradients_GrpLevel_' conn_type{1} '_Gradient_Var_Explained.eps'];
saveas(gcf,f2s,'epsc')




