% %%%
% %%% paper3_plotFreshwaterBudget.m
% %%%
% %%% Plots freshwater budget in quasi-latitude coordinates
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
% %%% Set true to use grounding line coordinate
% gl_coord = false;
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
% rho0 = rhoConst;
% Cp = 4000;
% rhoFresh = 1000;
% 
% %%% Define coordinate system for integrating to compute heatfunction
% ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity,gl_coord);
% eta = -9:.1:11;
% Neta = length(eta);
% 
% %%% Bounds and horizontal mask for heat budget volume
% eidx_icefront = find(abs(eta-eta_icefront)==min(abs(eta-eta_icefront)));
% eidx_shelfbreak = find(abs(eta-eta_shelfbreak)==min(abs(eta-eta_shelfbreak)));
% msk_wholeshelf = (ETA > eta_icefront) & (ETA < eta_shelfbreak);
% msk_slope = (ETA > eta_shelfbreak & ETA < 5); 
% msk_outershelf = (ETA > 1.5 & ETA < eta_shelfbreak); 
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
% %%% Compute volume averages
% EKE_avg = zeros(1,length(tt));
% EKE_slope = zeros(1,length(tt));
% EKE_outershelf = zeros(1,length(tt));
% PEtoEKE_avg = zeros(1,length(tt));
% PEtoEKE_slope = zeros(1,length(tt));
% PEtoEKE_outershelf = zeros(1,length(tt));
% MKEtoEKE_avg = zeros(1,length(tt));
% MKEtoEKE_outershelf = zeros(1,length(tt));
% salt_upper_avg = zeros(1,length(tt));
% salt_lower_avg = zeros(1,length(tt));
% theta_upper_avg = zeros(1,length(tt));
% theta_lower_avg = zeros(1,length(tt));
% theta_pos_lower_avg = zeros(1,length(tt));
% theta_neg_lower_avg = zeros(1,length(tt));
% volW_lower = sum(sum(sum(RAW.*msk_wholeshelf.*DRF(:,:,zidx_icefront:end).*hFacW(:,:,zidx_icefront:end))));
% volS_lower = sum(sum(sum(RAS.*msk_wholeshelf.*DRF(:,:,zidx_icefront:end).*hFacS(:,:,zidx_icefront:end))));
% volC_lower = sum(sum(sum(RAC.*msk_wholeshelf.*DRF(:,:,zidx_icefront:end).*hFacC(:,:,zidx_icefront:end))));
% volW_upper = sum(sum(sum(RAW.*msk_wholeshelf.*DRF(:,:,1:zidx_icefront-1).*hFacW(:,:,1:zidx_icefront-1))));
% volS_upper = sum(sum(sum(RAS.*msk_wholeshelf.*DRF(:,:,1:zidx_icefront-1).*hFacS(:,:,1:zidx_icefront-1))));
% volC_upper = sum(sum(sum(RAC.*msk_wholeshelf.*DRF(:,:,1:zidx_icefront-1).*hFacC(:,:,1:zidx_icefront-1))));
% volC_slope = sum(sum(sum(RAC.*msk_slope.*DRF.*hFacC)));
% volW_slope = sum(sum(sum(RAW.*msk_slope.*DRF.*hFacW)));
% volS_slope = sum(sum(sum(RAS.*msk_slope.*DRF.*hFacS)));
% volC_outershelf = sum(sum(sum(RAC.*msk_outershelf.*DRF.*hFacC)));
% volW_outershelf = sum(sum(sum(RAW.*msk_outershelf.*DRF.*hFacW)));
% volS_outershelf = sum(sum(sum(RAS.*msk_outershelf.*DRF.*hFacS)));
% for n=1:length(tt)
%   EKE_avg(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk_wholeshelf)) / (volW_upper+volW_lower)  ...
%              + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk_wholeshelf)) / (volS_upper+volW_lower);
%   EKE_slope(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk_slope)) / (volW_slope)  ...
%              + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk_slope)) / (volS_slope);  
%   EKE_outershelf(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk_outershelf)) / (volW_outershelf)  ...
%              + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk_outershelf)) / (volS_outershelf);  
%   PEtoEKE_avg(n) = sum(sum(PEtoEKE_int(:,:,n).*RAC.*msk_wholeshelf)) / (volC_upper+volC_lower);
%   PEtoEKE_slope(n) = sum(sum(PEtoEKE_int(:,:,n).*RAC.*msk_slope)) / (volC_slope);
%   PEtoEKE_outershelf(n) = sum(sum(PEtoEKE_int(:,:,n).*RAC.*msk_outershelf)) / (volC_outershelf);
%   MKEtoEKE_outershelf(n) = sum(sum(MKEtoEKE_int(:,:,n).*RAC.*msk_outershelf)) / (volC_outershelf);
%   MKEtoEKE_avg(n) = sum(sum(MKEtoEKE_int(:,:,n).*RAC.*msk_wholeshelf)) / (volC_upper+volC_lower);             
%   salt_lower_avg(n) = sum(sum(salt_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
%   salt_upper_avg(n) = sum(sum(salt_int_upper(:,:,n).*RAC.*msk_wholeshelf)) / volC_upper; 
%   theta_lower_avg(n) = sum(sum(theta_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
%   theta_pos_lower_avg(n) = sum(sum(theta_pos_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
%   theta_neg_lower_avg(n) = sum(sum(theta_neg_int_lower(:,:,n).*RAC.*msk_wholeshelf)) / volC_lower; 
%   theta_upper_avg(n) = sum(sum(theta_int_upper(:,:,n).*RAC.*msk_wholeshelf)) / volC_upper; 
% end
% 
% %%% Area-integrated vertical heat flux
% ws_tot_mod = ws_tot_flux - w_flux.*salt0;
% ws_mean_mod = ws_mean_flux - w_flux.*salt0;
% 
% %%% Heat function
% psiS_mod = psiS_tot-psi_tot*salt0;
% psiS_mean_mod = psiS_mean-psi_tot*salt0;
% psiS_eddy_adv_mod = psiS_eddy_adv - psi_eddy*salt0;
% psiS_eddy_stir = psiS_eddy - psiS_eddy_adv_mod;
% 
% %%% Decomposition of vertical salt flux
% sflux_vert = squeeze(sum(sum(ws_tot_mod.*RAC.*msk_wholeshelf,1),2) / -salt0 * rhoFresh*t1year / 1e12); %%% Convert to Gt/yr);
% sflux_mean_vert = squeeze(sum(sum(ws_mean_mod.*RAC.*msk_wholeshelf,1),2)/ -salt0 * rhoFresh*t1year / 1e12);
% sflux_eddy_vert = squeeze(sum(sum(ws_eddy_flux.*RAC.*msk_wholeshelf,1),2)/ -salt0 * rhoFresh*t1year / 1e12);
% sflux_eddy_adv_vert = squeeze(sum(sum(w_eddy_flux.*(salt_bnd-salt0).*RAC.*msk_wholeshelf,1),2)/ -salt0 * rhoFresh*t1year / 1e12);
% sflux_eddy_stir_vert = sflux_eddy_vert - sflux_eddy_adv_vert;
% 
% %%% Decomposition of cross-shelf break salt flux
% sflux_shelfbreak = squeeze(-psiS_mod(eidx_shelfbreak,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% sflux_mean_shelfbreak = squeeze(-psiS_mean_mod(eidx_shelfbreak,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% sflux_eddy_shelfbreak = squeeze(-psiS_eddy(eidx_shelfbreak,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% sflux_eddy_adv_shelfbreak = squeeze(-psiS_eddy_adv_mod(eidx_shelfbreak,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% sflux_eddy_stir_shelfbreak = squeeze(-psiS_eddy_stir(eidx_shelfbreak,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% 
% %%% Decomposition of cross-ice front salt flux
% sflux_icefront = squeeze(psiS_mod(eidx_icefront,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% sflux_mean_icefront = squeeze(psiS_mean_mod(eidx_icefront,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% sflux_eddy_icefront = squeeze(psiS_eddy(eidx_icefront,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% sflux_eddy_adv_icefront = squeeze(psiS_eddy_adv_mod(eidx_icefront,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% sflux_eddy_stir_icefront = squeeze(psiS_eddy_stir(eidx_icefront,zidx_icefront,:))/ -salt0 * rhoFresh*t1year / 1e12;
% 
% %%% Estimate heat tendency
% salt_fit = polyfit(tt,salt_lower_avg,1);
% salt_lin = salt_fit(1)*tt + salt_fit(2); %%% Linear best fit to monthly mean time series
% salt_res = salt_lower_avg - salt_lin; %%% Residual from linear fit
% Nt = length(tt);
% M = zeros(Nt,Nt);
% M(1,1) = 3/4;
% M(1,2) = 1/8;
% M(1,Nt) = 1/8;
% for n=2:Nt-1
%   M(n,n-1) = 1/8;
%   M(n,n) = 3/4;
%   M(n,n+1) = 1/8;
% end
% M(Nt,Nt-1) = 1/8;
% M(Nt,Nt) = 3/4;
% M(Nt,1) = 1/8;
% % for n=1:Nt-1
% %   M(n,n) = 0.5;
% %   M(n,n+1) = 0.5;
% % end
% % M(Nt,Nt) = 1;
% % % M(Nt,1) = 0.5;
% salt_lower_inst = M \ salt_res';
% stend = (salt_lower_inst([2:end 1]) - salt_lower_inst([end 1:end-1]))/(2*(times(2)-times(1)))*volC_lower/ -salt0 * rhoFresh*t1year / 1e12;
% stend = stend + salt_fit(1)*volC_lower/ -salt0 * rhoFresh*t1year / 1e12; %%% Add back in linear trend
% 
% 
% %%% Mask for ice/land
% psiS_tot_mean = mean(psiS_tot-psi_tot*salt0,3)* rho0/1e9;
% msk = ones(size(psiS_tot_mean));
% msk_ice = NaN*msk;
% for j=1:Neta  
%   idx = find(psiS_tot_mean(j,:)==psiS_tot_mean(j,1));
%   idx(end) = [];
%   msk(j,idx) = NaN;
%   if (~isempty(idx))
%     msk_ice(j,1:idx(end)) = 1;
%   end
%   idx = find(abs(psiS_tot_mean(j,:))<1e-12,1,'first');
%   msk(j,idx+1:end) = NaN;
% end
% 
% stot = - sflux_vert + sflux_icefront + sflux_shelfbreak - stend;
% 
% 

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

fshelfcolor = colororder(4,:);
fvertcolor = colororder(2,:);
fcavitycolor = colororder(1,:);
ftendcolor = colororder(3,:);
frescolor = colororder(7,:);
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
figure(205)
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

psiS_tot_plot = -mean(psiS_mod,3) * rho0/1e9 .* msk;
axes('Position',axpos(1,:));
pcolor(EE,ZZ,psiS_tot_plot*0);
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
% set(ax3,'XTick',eticklocs);
% set(ax3,'XTickLabel',eticklabels);
box off;
set(ax3,'XLim',xlim-77);
set(ax3,'FontSize',fontsize);
set(get(ax3,'XLabel'),'String','Reference latitude');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SURFACE HEAT FLUX MAP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plotting options
clim = [-1 1]*2e1;
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
vhf_plot = -mean(ws_eddy_flux,3)/salt0*t1year;
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
title(h,'m/yr','Fontsize',fontsize,'interpreter','latex');

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
fluxlist = {'F_s_h_e_l_f','F_v_e_r_t','F_c_a_v_i_t_y','F_t_e_n_d','Residual'};
thelabels = categorical(fluxlist);
thelabels = reordercats(thelabels,fluxlist);
thebars = [mean(sflux_shelfbreak),mean(sflux_vert),mean(sflux_icefront),mean(stend),mean(stot)];
thestd = [std(sflux_shelfbreak),std(sflux_vert),std(sflux_icefront),std(stend),std(stot)];
b = bar(thelabels(1),thebars(1));        
set(b,'FaceColor',fshelfcolor);
hold on;
b = bar(thelabels(2),thebars(2));        
set(b,'FaceColor',fvertcolor);
b = bar(thelabels(3),thebars(3));        
set(b,'FaceColor',fcavitycolor);
b = bar(thelabels(4),thebars(4));        
set(b,'FaceColor',ftendcolor);
b = bar(thelabels(5),thebars(5));        
set(b,'FaceColor',frescolor);
set(gca,'FontSize',fontsize);
hold on
er = errorbar(thelabels,thebars,thestd,thestd);    
hold off;
er.Color = [0 0 0];
er.LineStyle = 'None';
ylabel('Freshwater flux (Gt/yr)');
set(gca,'YLim',[-300 1000]);


axes('Position',axpos(4,:));
fluxlist = {'F_s_h_e_l_f','F_v_e_r_t','F_c_a_v_i_t_y'};
thelabels = categorical(fluxlist);
thelabels = reordercats(thelabels,fluxlist);
thebars = [mean(sflux_mean_shelfbreak),mean(sflux_eddy_adv_shelfbreak),mean(sflux_eddy_stir_shelfbreak);
           mean(sflux_mean_vert),mean(sflux_eddy_adv_vert),mean(sflux_eddy_stir_vert);           
           mean(sflux_mean_icefront),mean(sflux_eddy_adv_icefront),mean(sflux_eddy_stir_icefront)];
thestd = [ std(sflux_mean_shelfbreak),std(sflux_eddy_adv_shelfbreak),std(sflux_eddy_stir_shelfbreak);
           std(sflux_mean_vert),std(sflux_eddy_adv_vert),std(sflux_eddy_stir_vert);
           std(sflux_mean_icefront),std(sflux_eddy_adv_icefront),std(sflux_eddy_stir_icefront)];
tmp = [1 2 3];
b2 = bar(tmp,thebars);      
b2(1).LineStyle = '--';
b2(2).LineStyle = ':';
b2(3).LineStyle = '-.';
b2(1).LineWidth = 1.5;
b2(2).LineWidth = 1.5;
b2(3).LineWidth = 1.5;
theedgecolor = [.3 .3 .3];
barcolors1 = 1-(1-[fshelfcolor;fvertcolor;fcavitycolor])*.75;
barcolors2 = 1-(1-[fshelfcolor;fvertcolor;fcavitycolor])*.5;
barcolors3 = 1-(1-[fshelfcolor;fvertcolor;fcavitycolor])*.25;
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
set(gca,'YLim',[-300 1000]);





%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None','FontWeight','bold');
annotation('textbox',[axpos(2,1) axpos(2,2)+0.01 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None','FontWeight','bold');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None','FontWeight','bold');
annotation('textbox',[axpos(4,1)-0.04 axpos(4,2)-0.05 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None','FontWeight','bold');
annotation('textbox',[axpos(5,1)-0.05 axpos(5,2)-0.05 0.03 0.03],'String',axlabels{5},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None','FontWeight','bold');
annotation('textbox',[axpos(6,1)-0.05 axpos(6,2)-0.04 0.03 0.03],'String',axlabels{6},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None','FontWeight','bold');






%%%%%%%%%%%%%%%%%%%
%%% TIME SERIES %%%
%%%%%%%%%%%%%%%%%%%



%%% Panel 1: time series of shelf-averaged salinity or surface buoyancy
%%% loss, accompanied by time series of terms in heat budget
axes('Position',axpos(5,:));
tyears = startyear+times/t1year;
plot(tyears,sflux_icefront,'LineWidth',1.5,'Color',fcavitycolor);
hold on;
plot(tyears,sflux_shelfbreak,'LineWidth',1.5,'Color',fshelfcolor);
clear ax;
ax(1) = gca;
ax(2) = axes('Position',get(ax(1),'Position'));
h2 = plot(ax(2),tyears,(1-salt_lower_avg/salt0)*1e3,'Color',saltcolor);
set(ax(1),'YAxisLocation','Left');
set(ax(2),'YAxisLocation','Right');
set(ax(2),'XAxisLocation','Top');
set(ax(2),'Color','None');
set(h2,'LineWidth',1.5);
hold on;
area(ax(2),[2011 2012],flip(1-[34.64,34.80;34.64,34.80]/salt0)*1e3,'FaceColor','y','FaceAlpha',0.25);
hold off;
legend(ax(1),'F_c_a_v_i_t_y','F_s_h_e_l_f','Location','NorthWest');
set(ax(1),'XLim',[2009 2015])
set(ax(2),'XLim',[2009 2015])
set(ax(2),'YLim',flip(1-[34.64 34.80]/salt0)*1e3);
% set(ax(2),'YTick',[34.46:0.02:34.62]);
% xlabel(ax(1),'Year');
set(ax(2),'XTick',[]);
set(ax,'Box','off');
set(ax,'FontSize',fontsize);
ylabel(ax(1),'Freshwater flux (Gt/yr)');
ylabel(ax(2),'Shelf-averaged freshwater fraction (x 10^-^3)');
% text(ax(2),2008,34.46,'Spin-up period','FontSize',fontsize);
% text(ax(2),2011,34.46,'12-hourly output','FontSize',fontsize);
ax(1).YColor = fshelfcolor;
ax(2).YColor = saltcolor;


%%% Panel 2: Time series of EKE and EKE production terms in outer shelf
%%% Shade period of high coastal polynya productivity
axes('Position',axpos(6,:));
tyears = startyear+times/t1year;
plot(tyears,sflux_vert,'LineWidth',1.5,'Color',fvertcolor);
hold on;
area([2011 2012],[0,2000;0,2000],'FaceColor','y','FaceAlpha',0.25);
hold off;
clear ax;
ax(1) = gca;
ax(2) = axes('Position',get(ax(1),'Position'));
% plot(ax(2),tyears,1000*PEtoEKE_outershelf*1e6,'Color',bcprodcolor,'LineWidth',1.5);
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
set(ax(2),'YLim',[-1 25]);
% set(ax(2),'YTick',[34.46:0.02:34.62]);
xlabel(ax(1),'Year');
set(ax(2),'XTick',[]);
set(ax,'Box','off');
set(ax,'FontSize',fontsize);
ylabel(ax(1),'Vertical freshwater flux (Gt/yr)');
ylabel(ax(2),'Eddy energy production (10^-^6 W/m^3)');
% text(ax(2),2008,34.46,'Spin-up period','FontSize',fontsize);
% text(ax(2),2011,34.46,'12-hourly output','FontSize',fontsize);
ax(2).YColor = bcprodcolor;
ax(1).YColor = fvertcolor;
legend(ax(2),'Baroclinic production','Barotropic production','Location','NorthEast');




%%% ANNOTATIONS %%%

figure1 = gcf;




% Create arrow
annotation(figure1,'arrow',[0.274 0.348],...
  [0.879186465082794 0.879186465082795],'Color',1-(1-fcavitycolor)*1,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create textbox
annotation(figure1,'textbox',...
  [0.59 0.973197264218865 0.3385 0.0242980561555076],...
  'String',{'Vertical eddy freshwater flux at 246m depth'},...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none');

% Create arrow
annotation(figure1,'arrow',[0.471 0.408000000000002],...
  [0.91 0.910188984881213],'Color',1-(1-fshelfcolor)*1,'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create arrow
annotation(figure1,'arrow',[0.373 0.373],...
  [0.912500000000001 0.963342332613395],'Color',1-(1-fvertcolor)*1,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create textbox
annotation(figure1,'textbox',...
  [0.397 0.935697264218866 0.0499999999999997 0.0242980561555076],...
  'String','$F_\mathrm{vert}$',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
   'Color',fvertcolor);

% Create textbox
annotation(figure1,'textbox',...
  [0.38 0.874863930885532 0.0499999999999997 0.0242980561555076],...
  'String','$F_\mathrm{shelf}$',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
   'Color',fshelfcolor);

% Create textbox
annotation(figure1,'textbox',...
  [0.24 0.848197264218867 0.0499999999999998 0.0242980561555076],...
  'String','$F_\mathrm{cavity}$',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
   'Color',fcavitycolor);


