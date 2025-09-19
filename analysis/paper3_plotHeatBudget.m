%%%
%%% paper3_plotHeatBudget.m
%%%
%%% Plots heat/salt budget in quasi-latitude coordinates
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
startyear = 2008; %%% Year in which simulation starts
loadexp;

%%% Reference surface freezing temperature
theta0 = -1.9;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to use depth-averaged temperature as the coordinate system
use_meanT = false;

%%% set true to use grounding line coordinate
gl_coord = true;

%%% Set true to decompose eddy fluxes
calc_eddy_decomp = false;

%%% Cross-shelf locations of ice front and shelf break for heat budget calculation
if (gl_coord)
  eta_icefront = 0;
else
  eta_icefront = -1.1;
end
eta_shelfbreak = 3.5; 

%%% Index of the upper grid cell face dividing the upper and lower portions
%%% of the water column
if (gl_coord)
  zidx_icefront = 15;
else
  zidx_icefront = 25;
end

%%% Physical parameters
rho0 = 1000;
Cp = 4000;

%%% Define coordinate system for integrating to compute heatfunction
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity,gl_coord);
eta = -9:.05:11;
Neta = length(eta);

%%% Bounds and horizontal mask for heat budget volume
eidx_icefront = find(abs(eta-eta_icefront)==min(abs(eta-eta_icefront)));
eidx_shelfbreak = find(abs(eta-eta_shelfbreak)==min(abs(eta-eta_shelfbreak)));
msk_wholeshelf = (ETA > eta_icefront) & (ETA < eta_shelfbreak);
msk_slope = (ETA > eta_shelfbreak & ETA < 5); 
msk_outershelf = (ETA > 1.5 & ETA < eta_shelfbreak); 
msk_innershelf = (ETA > 0 & ETA < 1.5); 

%%% Load HeatFunction data file
outfname = [expname,'_HeatFunction'];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else 
  if (use_meanT)
    outfname = [outfname,'_meanT'];
  else 
    if (deform_cavity)
      outfname = [outfname,'_deform'];
    elseif (gl_coord)
      outfname = [outfname,'_GLcoord'];
    end
  end
end
load(fullfile('products',outfname));
Neta = length(eta);

%%% Load positive and negative heatfunction components
outfname = [expname,'_PosNegHeatFunction'];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else
  if (use_meanT)
    outfname = [outfname,'_meanT'];
  else 
    if (deform_cavity)
      outfname = [outfname,'_deform'];
    elseif (gl_coord)
      outfname = [outfname,'_GLcoord'];
    end
  end
end
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

%%% Load eddy-induced components of heatfunction/transport
outfname = [expname,'_HeatFunctionEddyDecomp'];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else
  if (deform_cavity)
    outfname = [outfname,'_deform'];
  elseif (gl_coord)
    outfname = [outfname,'_GLcoord'];
  end
end
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

%%% Load shelf heat budget diagnostics
outfname = [expname,'_ShelfHeatBudget'];
if (gl_coord)
  outfname = [outfname,'_GLcoord'];
end
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

%%% Compute volume averages
EKE_avg = zeros(1,length(tt));
EKE_slope = zeros(1,length(tt));
EKE_outershelf = zeros(1,length(tt));
PEtoEKE_avg = zeros(1,length(tt));
PEtoEKE_slope = zeros(1,length(tt));
PEtoEKE_outershelf = zeros(1,length(tt));
MKEtoEKE_avg = zeros(1,length(tt));
MKEtoEKE_outershelf = zeros(1,length(tt));
salt_upper_avg = zeros(1,length(tt));
salt_lower_avg = zeros(1,length(tt));
salt_outershelf_avg = zeros(1,length(tt));
salt_innershelf_avg = zeros(1,length(tt));
theta_upper_avg = zeros(1,length(tt));
theta_lower_avg = zeros(1,length(tt));
theta_upper_outer_avg = zeros(1,length(tt));
theta_lower_outer_avg = zeros(1,length(tt));
theta_pos_lower_avg = zeros(1,length(tt));
theta_neg_lower_avg = zeros(1,length(tt));
volW_lower = sum(sum(sum(RAW.*msk_wholeshelf.*DRF(:,:,zidx_icefront:end).*hFacW(:,:,zidx_icefront:end))));
volS_lower = sum(sum(sum(RAS.*msk_wholeshelf.*DRF(:,:,zidx_icefront:end).*hFacS(:,:,zidx_icefront:end))));
volC_lower = sum(sum(sum(RAC.*msk_wholeshelf.*DRF(:,:,zidx_icefront:end).*hFacC(:,:,zidx_icefront:end))));
volW_upper = sum(sum(sum(RAW.*msk_wholeshelf.*DRF(:,:,1:zidx_icefront-1).*hFacW(:,:,1:zidx_icefront-1))));
volS_upper = sum(sum(sum(RAS.*msk_wholeshelf.*DRF(:,:,1:zidx_icefront-1).*hFacS(:,:,1:zidx_icefront-1))));
volC_upper = sum(sum(sum(RAC.*msk_wholeshelf.*DRF(:,:,1:zidx_icefront-1).*hFacC(:,:,1:zidx_icefront-1))));
volC_upper_outer = sum(sum(sum(RAC.*msk_outershelf.*DRF(:,:,1:zidx_icefront-1).*hFacC(:,:,1:zidx_icefront-1))));
volC_lower_outer = sum(sum(sum(RAC.*msk_outershelf.*DRF(:,:,zidx_icefront:end).*hFacC(:,:,zidx_icefront:end))));
volC_slope = sum(sum(sum(RAC.*msk_slope.*DRF.*hFacC)));
volW_slope = sum(sum(sum(RAW.*msk_slope.*DRF.*hFacW)));
volS_slope = sum(sum(sum(RAS.*msk_slope.*DRF.*hFacS)));
volC_innershelf = sum(sum(sum(RAC.*msk_innershelf.*DRF.*hFacC)));
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
  MKEtoEKE_outershelf(n) = sum(sum(MKEtoEKE_int(:,:,n).*RAC.*msk_outershelf)) / (volC_outershelf);
  MKEtoEKE_avg(n) = sum(sum(MKEtoEKE_int(:,:,n).*RAC.*msk_wholeshelf)) / (volC_upper+volC_lower);             
  salt_lower_avg(n) = sum(sum(salt_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
  salt_upper_avg(n) = sum(sum(salt_int_upper(:,:,n).*RAC.*msk_wholeshelf)) / volC_upper; 
  salt_outershelf_avg(n) = sum(sum((salt_int_lower(:,:,n)+salt_int_upper(:,:,n)).*RAC.*msk_outershelf)) / volC_outershelf; 
  salt_innershelf_avg(n) = sum(sum((salt_int_lower(:,:,n)+salt_int_upper(:,:,n)).*RAC.*msk_innershelf)) / volC_innershelf; 
  theta_lower_avg(n) = sum(sum(theta_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
  theta_lower_outer_avg(n) = sum(sum(theta_int_lower(:,:,n).*RAC.*msk_outershelf)) / volC_lower_outer; 
  theta_pos_lower_avg(n) = sum(sum(theta_pos_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
  theta_neg_lower_avg(n) = sum(sum(theta_neg_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
  theta_upper_avg(n) = sum(sum(theta_int_upper(:,:,n).*RAC.*msk_wholeshelf)) / volC_upper; 
  theta_upper_outer_avg(n) = sum(sum(theta_int_upper(:,:,n).*RAC.*msk_outershelf)) / volC_upper_outer; 
end

%%% Area-integrated vertical heat flux
wt_tot_mod = wt_tot_flux - w_flux.*theta0;
wt_mean_mod = wt_mean_flux - w_flux.*theta0;

%%% Heat function
psiT_mod = psiT_tot-psi_tot*theta0;
psiT_mean_mod = psiT_mean-psi_tot*theta0;
psiT_eddy_adv_mod = psiT_eddy_adv - psi_eddy*theta0;
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
msk_ice_land = ones(size(psiT_tot_mean));
msk_ice = NaN*msk_ice_land;
for j=1:Neta  
  idx = find(psiT_tot_mean(j,:)==psiT_tot_mean(j,1));
  idx(end) = [];
  msk_ice_land(j,idx) = NaN;
  if (~isempty(idx))
    msk_ice(j,1:idx(end)) = 1;
  end
  idx = find(abs(psiT_tot_mean(j,:))<1e-12,1,'first');
  msk_ice_land(j,idx+1:end) = NaN;
end

%%% Residual heat flux
htot = - hflux_vert - hflux_icefront + hflux_shelfbreak - htend;



%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(5,4);
axpos(1,:) = [0.07 0.75 .4 .21];
axpos(2,:) = [0.52 0.69 .4 .28];
axpos(3,:) = [0.07 0.5 .25 .2];
axpos(4,:) = [0.37 0.5 .56 .2];
axpos(5,:) = [0.07 0.27 .86 .18];
axpos(6,:) = [0.07 0.04 .86 .18];
cbpos1 = [0.5 0.75 0.01 .23];
cbpos2 = [0.94 0.75 0.01 .23];
axlabels = {'\textbf{A}','\textbf{B}','\textbf{C}','\textbf{D}','\textbf{E}','\textbf{F}'};
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

qshelfcolor = colororder(4,:);
qvertcolor = colororder(2,:);
qcavitycolor = colororder(1,:);
qtendcolor = colororder(3,:);
qrescolor = colororder(7,:);
saltcolor = colororder(5,:);
bcprodcolor = colororder(6,:);
btprodcolor = [.7 .7 .7];

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
zminidx_icefront = find(~isnan(msk_ice_land(eidx_icefront,:)),1,'last');
eidx_shelfbreak = find(abs(eta-eta_shelfbreak)<1e-8);
zminidx_shelfbreak = find(isnan(msk_ice_land(eidx_shelfbreak,:)),1);



%%% TODO add time series of vertical eddy heat flux and eddy energy
%%% production, pos/neg heat flux into cavity

%%% TODO shade region used to compute averages of vertical buoyancy fluxes
%%% Schematic 

psiT_tot_plot = -mean(psiT_mod,3) * rho0*Cp/1e12 .* msk_ice_land;
axes('Position',axpos(1,:));
pcolor(EE,ZZ,psiT_tot_plot*0);
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
% cbhandle = colorbar;
% set(cbhandle,'Position',cbpos1);
% title(cbhandle,'TW','FontSize',fontsize);
xlabel('Cross-shelf coordinate, \eta');
ylabel('Depth (m)');
% title('Total heat function');
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

%%% Add reference latitudes
ax3 = axes('Position',get(ax1,'Position'));
set(ax3,'Color','None');
set(ax3,'XAxisLocation','Top');
set(ax3,'YAxisLocation','Right');
set(ax3,'YLim',get(ax1,'YLim'));
set(ax3,'YTick',[]);
box off;
set(ax3,'XLim',xlim-77);
set(ax3,'FontSize',fontsize);
set(get(ax3,'XLabel'),'String','Reference latitude');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VERTICAL HEAT FLUX MAP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plotting options
if (gl_coord)
  clim = [-50 50];
else
  clim = [-100 100];
end
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










%%%%%%%%%%%%%%%%%%
%%% BAR CHARTS %%%
%%%%%%%%%%%%%%%%%%


axes('Position',axpos(3,:));


fluxlist = {'Q_s_h_e_l_f','Q_v_e_r_t','Q_c_a_v_i_t_y','Q_t_e_n_d','Q_r_e_s'};
thelabels = categorical(fluxlist);
thelabels = reordercats(thelabels,fluxlist);
thebars = [mean(hflux_shelfbreak),mean(hflux_vert),mean(hflux_icefront),mean(htend),mean(htot)];
thestd = [std(hflux_shelfbreak),std(hflux_vert),std(hflux_icefront),std(htend),std(htot)];
b = bar(thelabels(1),thebars(1));        
set(b,'FaceColor',qshelfcolor);
hold on;
b = bar(thelabels(2),thebars(2));        
set(b,'FaceColor',qvertcolor);
b = bar(thelabels(3),thebars(3));        
set(b,'FaceColor',qcavitycolor);
b = bar(thelabels(4),thebars(4));        
set(b,'FaceColor',qtendcolor);
b = bar(thelabels(5),thebars(5));        
set(b,'FaceColor',qrescolor);
set(gca,'FontSize',fontsize);
hold on
er = errorbar(thelabels,thebars,thestd,thestd);    
hold off;
er.Color = [0 0 0];
er.LineStyle = 'None';
ylabel('Heat flux (TW)');



axes('Position',axpos(4,:));




fluxlist = {'Q_s_h_e_l_f','Q_v_e_r_t','Q_c_a_v_i_t_y'};
thelabels = categorical(fluxlist);
thelabels = reordercats(thelabels,fluxlist);
thebars = [mean(hflux_mean_shelfbreak),mean(hflux_eddy_adv_shelfbreak),mean(hflux_eddy_stir_shelfbreak);
           mean(hflux_mean_vert),mean(hflux_eddy_adv_vert),mean(hflux_eddy_stir_vert);           
           mean(hflux_mean_icefront),mean(hflux_eddy_adv_icefront),mean(hflux_eddy_stir_icefront)];
thestd = [ std(hflux_mean_shelfbreak),std(hflux_eddy_adv_shelfbreak),std(hflux_eddy_stir_shelfbreak);
           std(hflux_mean_vert),std(hflux_eddy_adv_vert),std(hflux_eddy_stir_vert);
           std(hflux_mean_icefront),std(hflux_eddy_adv_icefront),std(hflux_eddy_stir_icefront)];
tmp = [1 2 3];
b2 = bar(tmp,thebars);      
b2(1).LineStyle = '--';
b2(2).LineStyle = ':';
b2(3).LineStyle = '-.';
b2(1).LineWidth = 1.5;
b2(2).LineWidth = 1.5;
b2(3).LineWidth = 1.5;
theedgecolor = [.3 .3 .3];
barcolors1 = 1-(1-[qshelfcolor;qvertcolor;qcavitycolor])*.75;
barcolors2 = 1-(1-[qshelfcolor;qvertcolor;qcavitycolor])*.5;
barcolors3 = 1-(1-[qshelfcolor;qvertcolor;qcavitycolor])*.25;
% b2(1).EdgeColor = theedgecolor;
% b2(2).EdgeColor = theedgecolor;
% b2(3).EdgeColor = theedgecolor;
b2(1).CData = barcolors1;
b2(1).FaceColor = 'flat';
b2(2).CData = barcolors2;
b2(2).FaceColor = 'flat';
b2(3).CData = barcolors3;
b2(3).FaceColor = 'flat';
set(gca,'FontSize',fontsize);
hold on
ngroups = size(thebars, 1);
nbars = size(thebars, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
  x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
  er = errorbar(x, thebars(:,i), thestd(:,i), '.');
  er.Color = [0 0 0];
  er.LineStyle = 'None';
end
hold off
set(gca,'XTickLabel',fluxlist);
leghandle = legend('Mean flow','Eddy advection','Eddy Stirring');
set(gca,'YLim',[-1 4]);





%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1) axpos(2,2)+0.01 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.04 axpos(4,2)-0.05 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(5,1)-0.05 axpos(5,2)-0.05 0.03 0.03],'String',axlabels{5},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(6,1)-0.05 axpos(6,2)-0.04 0.03 0.03],'String',axlabels{6},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');






%%%%%%%%%%%%%%%%%%%
%%% TIME SERIES %%%
%%%%%%%%%%%%%%%%%%%



%%% Panel 1: time series of shelf-averaged salinity or surface buoyancy
%%% loss, accompanied by time series of terms in heat budget
axes('Position',axpos(5,:));
tyears = startyear+times/t1year;
plot(tyears,hflux_icefront,'LineWidth',1.5,'Color',qcavitycolor);
hold on;
plot(tyears,hflux_shelfbreak,'LineWidth',1.5,'Color',qshelfcolor);
plot(tyears,hflux_eddy_icefront,'--','LineWidth',1.5,'Color',qcavitycolor);
plot(tyears,hflux_eddy_shelfbreak,'--','LineWidth',1.5,'Color',qshelfcolor);
clear ax;
ax(1) = gca;
ax(2) = axes('Position',get(ax(1),'Position'));
h2 = plot(ax(2),tyears,salt_lower_avg,'Color',saltcolor);
set(ax(1),'YAxisLocation','Left');
set(ax(2),'YAxisLocation','Right');
set(ax(2),'XAxisLocation','Top');
set(ax(2),'Color','None');
set(h2,'LineWidth',1.5);
hold on;
area(ax(2),[2011 2012],[34.45,34.62;34.45,34.62],'FaceColor','y','FaceAlpha',0.25);
hold off;
legend(ax(1),'Q_c_a_v_i_t_y','Q_s_h_e_l_f','Location','NorthWest');
set(ax(1),'XLim',[2009 2015])
set(ax(2),'XLim',[2009 2015])
if (gl_coord)
  set(ax(2),'YLim',[34.55 34.71]);
else
  set(ax(2),'YLim',[34.64 34.80]);
end
% set(ax(2),'YTick',[34.46:0.02:34.62]);
% xlabel(ax(1),'Year');
set(ax(2),'XTick',[]);
set(ax,'Box','off');
set(ax,'FontSize',fontsize);
ylabel(ax(1),'Heat flux (TW)');
ylabel(ax(2),'Shelf-averaged salinity (g/kg)');
% text(ax(2),2008,34.46,'Spin-up period','FontSize',fontsize);
% text(ax(2),2011,34.46,'12-hourly output','FontSize',fontsize);
ax(1).YColor = qshelfcolor;
ax(2).YColor = saltcolor;


%%% Panel 2: Time series of EKE and EKE production terms in outer shelf
%%% Shade period of high coastal polynya productivity
axes('Position',axpos(6,:));
tyears = startyear+times/t1year;
plot(tyears,hflux_vert,'LineWidth',1.5,'Color',qvertcolor);
hold on;
area([2011 2012],[4;4],'FaceColor','y','FaceAlpha',0.25);
area([2011 2012],[-1;-1],'FaceColor','y','FaceAlpha',0.25);
hold off;
clear ax;
ax(1) = gca;
ax(2) = axes('Position',get(ax(1),'Position'));
plot(ax(2),tyears,1000*PEtoEKE_avg*1e6,'Color',bcprodcolor,'LineWidth',1.5);
hold on;
plot(ax(2),tyears,1000*MKEtoEKE_avg*1e6,'Color',btprodcolor,'LineWidth',1.5);
set(ax(1),'YAxisLocation','Left');
set(ax(2),'YAxisLocation','Right');
set(ax(2),'XAxisLocation','Top');
set(ax(2),'Color','None');
hold off;
set(ax(1),'XLim',[2009 2015])
set(ax(2),'XLim',[2009 2015])
if (gl_coord)
  set(ax(2),'YLim',[-8 32]);
end
% set(ax(2),'YTick',[34.46:0.02:34.62]);
xlabel(ax(1),'Year');
set(ax(2),'XTick',[]);
set(ax,'Box','off');
set(ax,'FontSize',fontsize);
ylabel(ax(1),'Vertical heat flux (TW)');
ylabel(ax(2),'Eddy energy production (10^-^6 W/m^3)');
% text(ax(2),2008,34.46,'Spin-up period','FontSize',fontsize);
% text(ax(2),2011,34.46,'12-hourly output','FontSize',fontsize);
ax(2).YColor = bcprodcolor;
ax(1).YColor = qvertcolor;
legend(ax(2),'Baroclinic production','Barotropic production','Location','NorthEast');

%%% Compute and remove mean seasonal cycles
hflux_eddy_vert_mss = zeros(1,12);
PEtoEKE_avg_mss = zeros(1,12);
dtheta_avg_mss = zeros(1,12);
dtheta_avg = theta_lower_outer_avg-theta_upper_outer_avg;
% dtheta_avg = theta_lower_avg-theta_upper_avg;
hflux_eddy_param = ((salt_innershelf_avg-salt_outershelf_avg).*EKE_avg.*dtheta_avg);
hflux_eddy_param_mss = zeros(1,12);
iter_cntr = zeros(1,12);
for n=1:length(hflux_eddy_vert)
  hflux_eddy_vert_mss(mod(n-1,12)+1) = hflux_eddy_vert_mss(mod(n-1,12)+1) + hflux_eddy_vert(n);
  PEtoEKE_avg_mss(mod(n-1,12)+1) = PEtoEKE_avg_mss(mod(n-1,12)+1) + PEtoEKE_avg(n);
  dtheta_avg_mss(mod(n-1,12)+1) = dtheta_avg_mss(mod(n-1,12)+1) + dtheta_avg(n);
  hflux_eddy_param_mss(mod(n-1,12)+1) = hflux_eddy_param_mss(mod(n-1,12)+1) + hflux_eddy_param(n);
  iter_cntr(mod(n-1,12)+1) = iter_cntr(mod(n-1,12)+1) + 1;
end
hflux_eddy_vert_mss = hflux_eddy_vert_mss./iter_cntr;
PEtoEKE_avg_mss = PEtoEKE_avg_mss./iter_cntr;
dtheta_avg_mss = dtheta_avg_mss./iter_cntr;
hflux_eddy_param_mss = hflux_eddy_param_mss./iter_cntr;
hflux_eddy_vert_noss= hflux_eddy_vert' - repmat(hflux_eddy_vert_mss,[1 length(hflux_eddy_vert)/12]);
PEtoEKE_avg_noss= PEtoEKE_avg - repmat(PEtoEKE_avg_mss,[1 length(PEtoEKE_avg)/12]);
dtheta_avg_noss= dtheta_avg - repmat(dtheta_avg_mss,[1 length(dtheta_avg)/12]);
hflux_eddy_param_noss= hflux_eddy_param - repmat(hflux_eddy_param_mss,[1 length(hflux_eddy_param)/12]);

[r,p] = corr(PEtoEKE_avg',hflux_eddy_vert)
[r,p] = corr(PEtoEKE_avg_noss',hflux_eddy_vert_noss')
[r,p] = corr(EKE_avg',hflux_eddy_vert)
[r,p] = corr(dtheta_avg',hflux_eddy_vert)
[r,p] = corr(dtheta_avg'.*EKE_avg',hflux_eddy_vert)
[r,p] = corr(dtheta_avg_noss',hflux_eddy_vert_noss')
[r,p] = corr(hflux_eddy_param',hflux_eddy_vert)
[r,p] = corr(hflux_eddy_param_noss',hflux_eddy_vert_noss')

%%% ANNOTATIONS %%%

figure1 = gcf;





% Create textbox
annotation(figure1,'textbox',...
  [0.629 0.973197264218865 0.3085 0.0242980561555076],...
  'String',{['Vertical eddy heat flux at ',num2str(round(-RF(zidx_icefront))),'m depth']},...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none');

if (gl_coord)

  % Create arrow
annotation(figure1,'arrow',[0.471 0.408000000000002],...
  [0.92 0.920188984881213],'Color',[0.494 0.184 0.556],'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create arrow
annotation(figure1,'arrow',[0.373 0.299],...
  [0.880853131749461 0.880853131749462],'Color',[0 0.447 0.741],...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create arrow
annotation(figure1,'arrow',[0.38 0.381],[0.9175 0.96167566594673],...
  'Color',[0.85 0.325 0.098],...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create textbox
annotation(figure1,'textbox',...
  [0.34 0.897363930885536 0.0499999999999997 0.0242980561555076],...
  'Color',[0.85 0.325 0.098],...
  'String','$Q_\mathrm{vert}$',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
  [0.4 0.875697264218866 0.0499999999999997 0.0242980561555076],...
  'Color',[0.494 0.184 0.556],...
  'String','$Q_\mathrm{shelf}$',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
  [0.269 0.847363930885537 0.0499999999999998 0.0242980561555076],...
  'Color',[0 0.447 0.741],...
  'String','$Q_\mathrm{cavity}$',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none');


else

% Create arrow
annotation(figure1,'arrow',[0.471 0.408000000000002],...
  [0.92 0.920188984881213],'Color',1-(1-qshelfcolor)*1,'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create arrow
annotation(figure1,'arrow',[0.373 0.373],...
  [0.932500000000001 0.983342332613395],'Color',1-(1-qvertcolor)*1,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create textbox
annotation(figure1,'textbox',...
  [0.377 0.975697264218866 0.0499999999999997 0.0242980561555076],...
  'String','$Q_\mathrm{vert}$',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
   'Color',qvertcolor);

% Create textbox
annotation(figure1,'textbox',...
  [0.38 0.884863930885532 0.0499999999999997 0.0242980561555076],...
  'String','$Q_\mathrm{shelf}$',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
   'Color',qshelfcolor);

% Create textbox
annotation(figure1,'textbox',...
  [0.24 0.848197264218867 0.0499999999999998 0.0242980561555076],...
  'String','$Q_\mathrm{cavity}$',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
   'Color',qcavitycolor);

% Create arrow
annotation(figure1,'arrow',[0.343 0.269],...
  [0.879186465082794 0.879186465082795],'Color',1-(1-qcavitycolor)*1,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

end

