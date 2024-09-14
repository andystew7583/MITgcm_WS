% %%%
% %%% paper3_plotHeatBudget.m
% %%%
% %%% Plots heat/salt budget in quasi-latitude coordinates
% %%%
% 
% %%% Load experiment
% expdir = '../experiments';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% startyear = 2008; %%% Year in which simulation starts
% loadexp;
% 
% %%% Reference surface freezing temperature
% theta0 = -1.9;
% 
% %%% Reference salinity
% salt0 = 34.6;
% % salt0 = 34.72;
% 
% %%% Cross-shelf locations of ice front and shelf break for heat budget calculation
% eta_icefront = -1.1;
% eta_shelfbreak = 3.5;  
% 
% %%% Depth of ice front for heat budget calculation
% zidx_icefront = 25;    
% 
% %%% Set true to deform coordinates in the cavity
% deform_cavity = false;
% 
% %%% Set true to use barotropic streamfunction as the coordinate system
% use_PsiBT = false;
% 
% %%% Set true to use depth-averaged temperature as the coordinate system
% use_meanT = false;
% 
% %%% Set true to decompose eddy fluxes
% calc_eddy_decomp = false;
% 
% %%% Physical parameters
% rho0 = 1000;
% Cp = 4000;
% 
% %%% Define coordinate system for integrating to compute heatfunction
% ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
% eta = -9:.1:11;
% Neta = length(eta);
% 
% %%% Bounds and horizontal mask for heat budget volume
% eidx_icefront = find(abs(eta-eta_icefront)==min(abs(eta-eta_icefront)));
% eidx_shelfbreak = find(abs(eta-eta_shelfbreak)==min(abs(eta-eta_shelfbreak)));
% msk_wholeshelf = (ETA > eta_icefront) & (ETA < eta_shelfbreak);
% msk_slope = (ETA > eta_shelfbreak & ETA < 5); 
% msk_outershelf = (ETA > 2 & ETA < eta_shelfbreak); 
% 
% %%% Load HeatFunction data file
% outfname = [expname,'_HeatFunction'];
% if (use_PsiBT)
%   outfname = [outfname,'_PsiBT'];
% else 
%   if (use_meanT)
%     outfname = [outfname,'_meanT'];
%   else 
%     if (deform_cavity)
%       outfname = [outfname,'_deform'];
%     end
%   end
% end
% load(fullfile('products',outfname));
% Neta = length(eta);
% 
% %%% Load positive and negative heatfunction components
% outfname = [expname,'_PosNegHeatFunction'];
% if (use_PsiBT)
%   outfname = [outfname,'_PsiBT'];
% else
%   if (use_meanT)
%     outfname = [outfname,'_meanT'];
%   else 
%     if (deform_cavity)
%       outfname = [outfname,'_deform'];
%     end
%   end
% end
% outfname = [outfname,'.mat'];
% load(fullfile('products',outfname));
% 
% %%% Load eddy-induced components of heatfunction/transport
% outfname = [expname,'_HeatFunctionEddyDecomp'];
% if (use_PsiBT)
%   outfname = [outfname,'_PsiBT'];
% else
%   if (deform_cavity)
%     outfname = [outfname,'_deform'];
%   end
% end
% outfname = [outfname,'.mat'];
% load(fullfile('products',outfname));
%  
% %%% Load shelf heat budget diagnostics
% outfname = [expname,'_ShelfHeatBudget.mat'];
% load(fullfile('products',outfname));
% 
%%% Compute volume averages
EKE_avg = zeros(1,length(tt));
EKE_slope = zeros(1,length(tt));
EKE_outershelf = zeros(1,length(tt));
PEtoEKE_avg = zeros(1,length(tt));
PEtoEKE_slope = zeros(1,length(tt));
PEtoEKE_outershelf = zeros(1,length(tt));
MKEtoEKE_avg = zeros(1,length(tt));
salt_upper_avg = zeros(1,length(tt));
salt_lower_avg = zeros(1,length(tt));
theta_upper_avg = zeros(1,length(tt));
theta_lower_avg = zeros(1,length(tt));
theta_pos_lower_avg = zeros(1,length(tt));
theta_neg_lower_avg = zeros(1,length(tt));
volW_lower = sum(sum(sum(RAW.*msk_wholeshelf.*DRF(:,:,zidx_icefront:end).*hFacW(:,:,zidx_icefront:end))));
volS_lower = sum(sum(sum(RAS.*msk_wholeshelf.*DRF(:,:,zidx_icefront:end).*hFacS(:,:,zidx_icefront:end))));
volC_lower = sum(sum(sum(RAC.*msk_wholeshelf.*DRF(:,:,zidx_icefront:end).*hFacC(:,:,zidx_icefront:end))));
volW_upper = sum(sum(sum(RAW.*msk_wholeshelf.*DRF(:,:,1:zidx_icefront-1).*hFacW(:,:,1:zidx_icefront-1))));
volS_upper = sum(sum(sum(RAS.*msk_wholeshelf.*DRF(:,:,1:zidx_icefront-1).*hFacS(:,:,1:zidx_icefront-1))));
volC_upper = sum(sum(sum(RAC.*msk_wholeshelf.*DRF(:,:,1:zidx_icefront-1).*hFacC(:,:,1:zidx_icefront-1))));
volC_slope = sum(sum(sum(RAC.*msk_slope.*DRF.*hFacC)));
volW_slope = sum(sum(sum(RAW.*msk_slope.*DRF.*hFacW)));
volS_slope = sum(sum(sum(RAS.*msk_slope.*DRF.*hFacS)));
volC_outershelf = sum(sum(sum(RAC.*msk_outershelf.*DRF.*hFacC)));
volW_outershelf = sum(sum(sum(RAW.*msk_outershelf.*DRF.*hFacW)));
volS_outershelf = sum(sum(sum(RAS.*msk_outershelf.*DRF.*hFacS)));
for n=1:length(tt)
  EKE_avg(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk_wholeshelf)) / (volW_upper+volW_lower)  ...
             + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk_wholeshelf)) / (volS_upper+volW_lower);
  EKE_slope(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk_slope)) / (volW_slope)  ...
             + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk_slope)) / (volS_slope);  
  EKE_outershelf(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk_outershelf)) / (volW_outershelf)  ...
             + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk_outershelf)) / (volS_outershelf);  
  PEtoEKE_avg(n) = sum(sum(PEtoEKE_int(:,:,n).*RAC.*msk_wholeshelf)) / (volC_upper+volC_lower);
  PEtoEKE_slope(n) = sum(sum(PEtoEKE_int(:,:,n).*RAC.*msk_slope)) / (volC_slope);
  PEtoEKE_outershelf(n) = sum(sum(PEtoEKE_int(:,:,n).*RAC.*msk_outershelf)) / (volC_outershelf);
  MKEtoEKE_avg(n) = sum(sum(MKEtoEKE_int(:,:,n).*RAC.*msk_wholeshelf)) / (volC_upper+volC_lower);             
  salt_lower_avg(n) = sum(sum(salt_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
  salt_upper_avg(n) = sum(sum(salt_int_upper(:,:,n).*RAC.*msk_wholeshelf)) / volC_upper; 
  theta_lower_avg(n) = sum(sum(theta_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
  theta_pos_lower_avg(n) = sum(sum(theta_pos_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
  theta_neg_lower_avg(n) = sum(sum(theta_neg_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
  theta_upper_avg(n) = sum(sum(theta_int_upper(:,:,n).*RAC.*msk_wholeshelf)) / volC_upper; 
end

%%% Area-integrated vertical heat flux
wt_tot_mod = wt_tot_flux - w_flux.*theta0;
wt_mean_mod = wt_mean_flux - w_flux.*theta0;

%%% Heat function
psiT_mod = psiT_tot-psi_tot*theta0;
psiT_mean_mod = psiT_mean-psi_tot*theta0;
% psiT_eddy_adv_mod = psiT_eddy_adv - psi_eddy*theta0;
psiT_eddy_adv_mod = psiT_eddy_adv; %%% TODO FIX THIS!
psiT_eddy_stir = psiT_eddy - psiT_eddy_adv_mod;

%%% Decomposition of vertical heat flux
hflux_vert = squeeze(sum(sum(wt_tot_mod.*RAC.*msk_wholeshelf,1),2)*rho0*Cp/1e12);
hflux_mean_vert = squeeze(sum(sum(wt_mean_mod.*RAC.*msk_wholeshelf,1),2)*rho0*Cp/1e12);
hflux_eddy_vert = squeeze(sum(sum(wt_eddy_flux.*RAC.*msk_wholeshelf,1),2)*rho0*Cp/1e12);
hflux_eddy_adv_vert = squeeze(sum(sum(w_eddy_flux.*(theta_bnd-theta0).*RAC.*msk_wholeshelf,1),2)*rho0*Cp/1e12);
hflux_eddy_stir_vert = hflux_eddy_vert - hflux_eddy_adv_vert;

%%% Decomposition of cross-shelf break heat flux
hflux_shelfbreak = squeeze(-psiT_mod(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;
hflux_mean_shelfbreak = squeeze(-psiT_mean_mod(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;
hflux_eddy_shelfbreak = squeeze(-psiT_eddy(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;
hflux_eddy_adv_shelfbreak = squeeze(-psiT_eddy_adv_mod(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;
hflux_eddy_stir_shelfbreak = squeeze(-psiT_eddy_stir(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

%%% Decomposition of cross-ice front heat flux
hflux_icefront = squeeze(-psiT_mod(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_mean_icefront = squeeze(-psiT_mean_mod(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_eddy_icefront = squeeze(-psiT_eddy(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_eddy_adv_icefront = squeeze(-psiT_eddy_adv_mod(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_eddy_stir_icefront = squeeze(-psiT_eddy_stir(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;

%%% "Positive" and "negative" components of heat flux
hflux_pos_icefront = squeeze(-psiT_pos_mean(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_pos_shelfbreak = squeeze(-psiT_pos_mean(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;
hflux_neg_icefront = squeeze(-psiT_neg_mean(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_neg_shelfbreak = squeeze(-psiT_neg_mean(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

%%% Estimate heat tendency
theta_fit = polyfit(tt,theta_lower_avg,1);
theta_lin = theta_fit(1)*tt + theta_fit(2); %%% Linear best fit to monthly mean time series
theta_res = theta_lower_avg - theta_lin; %%% Residual from linear fit
Nt = length(tt);
M = zeros(Nt,Nt);
M(1,1) = 3/4;
M(1,2) = 1/8;
M(1,Nt) = 1/8;
for n=2:Nt-1
  M(n,n-1) = 1/8;
  M(n,n) = 3/4;
  M(n,n+1) = 1/8;
end
M(Nt,Nt-1) = 1/8;
M(Nt,Nt) = 3/4;
M(Nt,1) = 1/8;
% for n=1:Nt-1
%   M(n,n) = 0.5;
%   M(n,n+1) = 0.5;
% end
% M(Nt,Nt) = 1;
% % M(Nt,1) = 0.5;
theta_lower_inst = M \ theta_res';
htend = (theta_lower_inst([2:end 1]) - theta_lower_inst([end 1:end-1]))/(2*(times(2)-times(1)))*volC_lower/1e12*rho0*Cp;
% htend = (theta_lower_avg([2:end 1]) - theta_lower_avg([end 1:end-1]))/(2*(times(2)-times(1)))*volC_lower/1e12*rho0*Cp;
htend = htend + theta_fit(1)*volC_lower/1e12*rho0*Cp; %%% Add back in linear trend


%%% Mask for ice/land
psiT_tot_mean = mean(psiT_tot-psi_tot*theta0,3)* rho0*Cp/1e12;
msk = ones(size(psiT_tot_mean));
msk_ice = NaN*msk;
for j=1:Neta  
  idx = find(psiT_tot_mean(j,:)==psiT_tot_mean(j,1));
  idx(end) = [];
  msk(j,idx) = NaN;
  if (~isempty(idx))
    msk_ice(j,1:idx(end)) = 1;
  end
  idx = find(abs(psiT_tot_mean(j,:))<1e-12,1,'first');
  msk(j,idx+1:end) = NaN;
end

htot = - hflux_vert - hflux_icefront + hflux_shelfbreak - htend;

%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(5,4);
axpos(1,:) = [0.07 0.75 .4 .23];
axpos(2,:) = [0.55 0.7 .4 .28];
axpos(3,:) = [0.07 0.5 .25 .2];
axpos(4,:) = [0.37 0.5 .56 .2];
axpos(5,:) = [0.07 0.25 .86 .2];
axpos(6,:) = [0.07 0.06 .86 .2];
cbpos1 = [0.5 0.75 0.01 .23];
cbpos2 = [0.96 0.75 0.01 .23];
axlabels = {'(a)','(b)','(c)','(d)','(e)'};
rho0 = 1027;
Cp = 4000;
arrowcolor = [0.301960784313725 0.55098039215686 0.93333333333333];
colororder = get(gca,'ColorOrder');
linewidth = 1.5;
ylim = [0 2000];
xlim = [-8.6 4];
icecolor = [186 242 239]/255;
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
psimax = 4;
psistep = 0.2;
psisteps = [0:psistep:-psistep psistep:psistep:psimax];

%%% Plotting range for salinity figure
latMin_b = min(min(YC));
latMax_b = YC(1,end-spongethickness);
lonMin_b = min(min(XC));
lonMax_b = XC(end-spongethickness,1);




%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(204)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34  1000  1200]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HEAT BUDGET SCHEMATIC %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eidx_icefront = find(abs(eta-eta_icefront)<1e-8);
zminidx_icefront = find(isnan(msk(eidx_icefront,:)),1);
eidx_shelfbreak = find(abs(eta-eta_shelfbreak)<1e-8);
zminidx_shelfbreak = find(isnan(msk(eidx_shelfbreak,:)),1);



%%% TODO add time series of vertical eddy heat flux and eddy energy
%%% production, pos/neg heat flux into cavity

%%% TODO shade region used to compute averages of vertical buoyancy fluxes
%%% Schematic 

psiT_tot_plot = -mean(psiT_mod,3) * rho0*Cp/1e12 .* msk;
axes('Position',axpos(1,:));
pcolor(EE,ZZ,psiT_tot_plot);
shading flat;
hold on;
plot([eta_icefront,eta_shelfbreak],[-RC(zidx_icefront),-RC(zidx_icefront)],'k--','LineWidth',2);
plot([eta_icefront,eta_icefront],[-RC(zidx_icefront),-RC(zminidx_icefront)],'k--','LineWidth',2);
plot([eta_shelfbreak,eta_shelfbreak],[-RC(zidx_icefront),-RC(zminidx_shelfbreak)],'k--','LineWidth',2);
hold off;
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(gca,cmocean('amp',length(psisteps)));
cbhandle = colorbar;
set(cbhandle,'Position',cbpos1);
title(cbhandle,'TW','FontSize',fontsize);
xlabel('Cross-shelf coordinate, \eta');
ylabel('Depth (m)');
title('Total heat function');
set(gca,'FontSize',fontsize);
caxis([0 psimax]);

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
pcolor(EE,ZZ,msk_ice);
shading flat;
colormap(ax2,icecolor);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
set(ax2,'YLim',ylim);
set(ax2,'YDir','reverse');
set(ax2,'XLim',xlim);
set(ax2,'Color','None')








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SURFACE HEAT FLUX MAP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plotting options
clim = [-100 100];
cmap = cmocean('balance',50);
% cmap = cmap(8:23,:);
% cmap = cmocean('thermal',16);
% cmap = cmocean('amp',16);
% cmap = pmkmp(16,'swtth');
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];


%%% Set up map plot
subplot('Position',axpos(2,:)+[.1 0 0 0]);
axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_b latMax_b], ...
  'MapLonLimit',[lonMin_b lonMax_b], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', [-70:10:-30],...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');  
axis off;
setm(gca,'MLabelParallel',-20)

%%% Plot surface heat flux
vhf_plot = mean(wt_eddy_flux,3)*rho0*Cp;
vhf_plot(sum(hFacC,3)==0) = NaN;
pcolorm(YC,XC,vhf_plot);
shading interp
caxis(clim);
colormap(gca,cmap);


%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',cbpos2)
tightmap;
title(h,'W/m$^2$','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.05 .05 .05],'BackgroundColor','none','Edgecolor','none')       
% Nlon = 101;
% dLon = (lonMax_c-lonMin_c)/(Nlon-1);
% plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'r--','LineWidth',2);
hold off;

%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos(2,:));
hold off

% Create textbox
annotation(gcf,'textbox',...
  [0.163 0.949863930885531 0.2685 0.0242980561555075],...
  'String',{'Vertical heat flux at 255m depth'},'EdgeColor','None','FontSize',fontsize+2,'interpreter','latex');











%%%%%%%%%%%%%%%%%%
%%% BAR CHARTS %%%
%%%%%%%%%%%%%%%%%%


axes('Position',axpos(3,:));


fluxlist = {'Q_s_h_e_l_f','Q_v_e_r_t','Q_c_a_v_i_t_y','Q_t_e_n_d','Residual'};
thelabels = categorical(fluxlist);
thelabels = reordercats(thelabels,fluxlist);
thebars = [mean(hflux_shelfbreak),mean(hflux_vert),mean(hflux_icefront),mean(htend),mean(htot)];
bar(thelabels,thebars)        
set(gca,'FontSize',fontsize);




thecolororder = get(gca,'ColorOrder');




axes('Position',axpos(4,:));



fluxlist = {'Q_v_e_r_t','Q_s_h_e_l_f','Q_c_a_v_i_t_y'};
thelabels = categorical(fluxlist);
thelabels = reordercats(thelabels,fluxlist);
thebars = [mean(hflux_vert),mean(hflux_mean_vert),mean(hflux_eddy_vert),mean(hflux_eddy_adv_vert),mean(hflux_eddy_stir_vert);
          mean(hflux_shelfbreak),mean(hflux_mean_shelfbreak),mean(hflux_eddy_shelfbreak),mean(hflux_eddy_adv_shelfbreak),mean(hflux_eddy_stir_shelfbreak);
          mean(hflux_icefront),mean(hflux_mean_icefront),mean(hflux_eddy_icefront),mean(hflux_eddy_adv_icefront),mean(hflux_eddy_stir_icefront)];
bar(thelabels,thebars)        
set(gca,'FontSize',fontsize);




thecolororder = get(gca,'ColorOrder');





%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.04 axpos(4,2)-0.05 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(5,1)-0.04 axpos(5,2)-0.05 0.03 0.03],'String',axlabels{5},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');






%%%%%%%%%%%%%%%%%%%
%%% TIME SERIES %%%
%%%%%%%%%%%%%%%%%%%



figure(16);
ax = plotyy(startyear+times/t1year,PEtoEKE_outershelf*rho0*Cp,startyear+times/t1year,hflux_eddy_vert);
xlabel('Year');
ylabel(ax(1),'Baroclinic EKE production (W/m^3)');
ylabel(ax(2),'Vertical eddy heat flux (TW)');
set(gca,'Position',[0.1 0.1 0.8 0.85]);
set(gca,'FontSize',fontsize);
set(ax(2),'FontSize',fontsize);












%%% ANNOTATIONS %%%

figure1 = gcf;

% Create arrow
annotation(figure1,'arrow',[0.675 0.675],...
  [0.900647948164147 0.970842332613394],...
  'Color',arrowcolor,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create arrow
annotation(figure1,'arrow',[0.826000000000001 0.756000000000002],...
  [0.886609071274299 0.887688984881212],...
  'Color',arrowcolor,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create arrow
annotation(figure1,'arrow',[0.594 0.52],...
  [0.83585313174946 0.835853131749461],...
  'Color',arrowcolor,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

