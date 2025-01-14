

clear
close all

%% Toolboxes 

addpath(genpath('Code/'))  
addpath ~/Dropbox/Toolbox/Code/  
addpath ~/Dropbox/Toolbox/BCT/2019_03_03_BCT/
addpath ~/Dropbox/Toolbox/SWP/
addpath ~/Dropbox/Toolbox/plotSpread/
addpath Code/Util/

%% Load Data
load Data.mat

%% Create RFN
dn = 0.16; % density

fc_abs = abs(fc); 

RFN.fc_sc_rfn = create_RFN(fc_abs, sc, dn);
RFN.fc_ec_rfn = create_RFN(fc_abs, ec,dn);
RFN.fc_gc_rfn = create_RFN(fc_abs, gc,dn);
RFN.fc_gcdist_rfn = create_RFN(fc_abs, gc_dist, dn);
RFN.fc_rc_rfn = create_RFN(fc_abs, rc, dn);
RFN.fc_rcdist_rfn = create_RFN(fc_abs, rc_dist, dn);
 

%% Calculate Network Properties
nitr = 1000;
 
% Full Network
Net_Prop.GrpLevel.full = my_calc_netprop(fc);

% Masked Network
Net_Prop.GrpLevel.sc = my_calc_netprop(RFN.fc_sc_rfn);
Net_Prop.GrpLevel.ed = my_calc_netprop(RFN.fc_ec_rfn);
Net_Prop.GrpLevel.gc = my_calc_netprop(RFN.fc_gc_rfn);
Net_Prop.GrpLevel.gc_dist = my_calc_netprop(RFN.fc_gcdist_rfn);
Net_Prop.GrpLevel.rc = my_calc_netprop(RFN.fc_rc_rfn);
Net_Prop.GrpLevel.rc_dist = my_calc_netprop(RFN.fc_rcdist_rfn);

% Random Edges 
[Net_Prop.GrpLevel.rand,Net_Prop.GrpLevel.fc_rand_rfn] = my_calc_netprop_rand(fc,dn,nitr);

save('HCP_Corr_Based_Network_Prop.mat','Net_Prop')  


%% If all ready calculated Network Properties

load HCP_Corr_Based_Network_Prop.mat

%% 

Meas = {'Degree','Q','ClustCoef','SpecRad','PL','Synch','EigCent',}; 
Net = {'sc','ec','gc','rc','gc_dist','rc_dist'};
NetLbls = {'FC_S_C','FC_D_C','FC_G_C','FC_R_C','FC_G_C_d','FC_R_C_d'};

mk = {'.','.','s','s','d','d'};
mks = [20, 20, 8, 8, 8, 8];


%% Colors 
colors = lines(1);
colors(2,:) = [0.9290 0.6940 0.1250];
colors(3,:) = [ 0 0.5 0];
colors(4,:) = [ 0.6 0 0];
colors(5,:) = [ 0 0.75 0];
colors(6,:) = [ 1 0 0];

%% Figure 4
for m = 1:numel(Meas)
    
    if strcmp(Meas{m},'Degree') || strcmp(Meas{m},'EigCent') || strcmp(Meas{m},'ClustCoef')
       
        for n = 1:numel(Net)
            m1 = mean(Net_Prop.GrpLevel.(Net{n}).(Meas{m}));
            m2 = mean(Net_Prop.GrpLevel.full.(Meas{m}));
            delta(n,m) = my_delta(m1,m2);
            clear m1 m2
   
        end
        
    else
   
        for n = 1:numel(Net)
            m1 = Net_Prop.GrpLevel.(Net{n}).(Meas{m});
            m2 = Net_Prop.GrpLevel.full.(Meas{m});
            delta(n,m) = my_delta(m1,m2);
            clear m1 m2
        end

    end

end

delta_abs = abs(delta);

% Figures: Average Graph Properties
[~,idx] = sort(mean(delta),'descend');

delta_sort = delta_abs(:,idx);

close all
figure('Units','centimeters','Position',[10 10 10 7])
hold on
for i = 1:numel(Meas)
    errorbar(i,mean(delta_sort(:,i)),std(delta_sort(:,i))./sqrt(numel(Meas)),'ok','MarkerSize',10,'LineWidth',1)
end

ylim([0 35])
xlim([0.5 7.5])
set(gca,'XTick',1:numel(Meas),'XtickLabel',Meas(idx), ...
        'FontSize',10,'XColor','k','YColor','k')
ylabel('|\Delta| (%)','FontSize',12)

% f2s = ['Poster_Figures/HCP_Network_Measures_Devation.eps'];
% saveas(gcf,f2s,'epsc')


% Figures: Average RFN
close all
figure('Units','centimeters','Position',[10 10 10 7])
hold on
for i = 1:numel(NetLbls)
    errorbar(i,mean(delta_abs(i,:)),std(delta_abs(i,:))./sqrt(numel(NetLbls)),...
        'marker',mk(i),'Markersize',mks(i),'Color',colors(i,:),'MarkerFaceColor',colors(i,:))
end

ylim([0 40])
xlim([0.5 6.5])
set(gca,'XTick',1:numel(NetLbls),'XtickLabel',NetLbls, ...
        'FontSize',10,'XColor','k','YColor','k')
ylabel('|\Delta| (%)','FontSize',12)

% f2s = ['HCP_Corr_Figures/HCP_RFN_Devation.eps'];
% saveas(gcf,f2s,'epsc')



%% Figure 5

close all

for m = 1:numel(Meas)
    
    if strcmp(Meas{m},'Degree') || strcmp(Meas{m},'EigCent') || strcmp(Meas{m},'ClustCoef')
    
        r = repmat(Net_Prop.GrpLevel.rand.(Meas{m}),1,numel(Net)); 
        gr = my_delta(r,mean(Net_Prop.GrpLevel.full.(Meas{m})));
       
        for n = 1:numel(Net)
            m1 = mean(Net_Prop.GrpLevel.(Net{n}).(Meas{m}));
            m2 = mean(Net_Prop.GrpLevel.full.(Meas{m}));
            delta(n,1) = my_delta(m1,m2);
            clear m1 m2
   
        end
        
    else
   
        r = repmat(Net_Prop.GrpLevel.rand.(Meas{m}),1,numel(Net)); 
        gr = my_delta(r,Net_Prop.GrpLevel.full.(Meas{m}));
    
        for n = 1:numel(Net)
            m1 = Net_Prop.GrpLevel.(Net{n}).(Meas{m});
            m2 = Net_Prop.GrpLevel.full.(Meas{m});
            delta(n,1) = my_delta(m1,m2);
            clear m1 m2
        end

    end


    % Figures
    figure('Units','centimeters','Position',[10 10 9 6])
    hold on
    plotSpread(gr)
    
    for i = 1:numel(NetLbls)
       plot(i,delta(i),'marker',mk(i),'Markersize',mks(i),'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:))
    end
    
%     ylim([ylm(m,1) ylm(m,2)])
    xlim([0.5 numel(Net)+0.5])
    
    set(gca,'XTick',[], ...
        'FontSize',10,'XColor','k','YColor','k')
    
    ylabel('\Delta (%)','FontSize',10)
    title(Meas{m},'FontSize',12,'FontWeight','normal')
    
    set(gcf,'Renderer','painters')
%     f2s = ['HCP_Corr_Figures/HCP_GrpLevel_Corr_Connections_HCP-dMRI_' Meas{m} '_delat_percent_thresh-0.esp'];
%     saveas(gcf,f2s,'epsc')
    
    clear r delta zm zs

end


%% Figure 6D

% Degree
m1 = Net_Prop.GrpLevel.recep.Degree;
m2 = Net_Prop.GrpLevel.full.Degree;
delta(:,1) = abs(my_delta(m1,m2));

close all
figure('Units','centimeters','Position',[10 10 10 7])
hold on
plot(m2, delta,'.', ...
    'Markersize',10,'MarkerFaceColor',[0.5 0.5 0.5])

set(gca,'FontSize',10,'XColor','k','YColor','k')
ylabel('|\Delta| (%)','FontSize',12)
xlabel('Degree','FontSize',12)


% Clustering Coefficient
m1 = Net_Prop.GrpLevel.recep.ClustCoef;
m2 = Net_Prop.GrpLevel.full.ClustCoef;
delta(:,1) = abs(my_delta(m1,m2));

close all
figure('Units','centimeters','Position',[10 10 10 7])
hold on
plot(m2, delta,'.', ...
    'Markersize',10,'MarkerFaceColor',[0.5 0.5 0.5])

set(gca,'FontSize',10,'XColor','k','YColor','k')
ylabel('|\Delta| (%)','FontSize',12)
xlabel('ClustCoef','FontSize',12)

% Eigenvector centrality
m1 = Net_Prop.GrpLevel.recep.EigCent;
m2 = Net_Prop.GrpLevel.full.EigCent;
delta(:,1) = abs(my_delta(m1,m2));

close all
figure('Units','centimeters','Position',[10 10 10 7])
hold on
plot(Net_Prop.GrpLevel.recep.Degree, delta,'.', ...
    'Markersize',10,'MarkerFaceColor',[0.5 0.5 0.5])

set(gca,'FontSize',10,'XColor','k','YColor','k')
ylabel('|\Delta| (%)','FontSize',12)
xlabel('EigCent','FontSize',12)



