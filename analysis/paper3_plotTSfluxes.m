%%%
%%% paper3_plotTSfluxes.m
%%%
%%% Plots heat/salt fluxes in quasi-latitude coordinates
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% loadexp;

%%% Options (see calcTSfluxes)
deform_cavity = false;
gl_coord = true;
outfname = ['','_TSfluxes'];
if (deform_cavity)
  outfname = [outfname,'_deform'];
elseif (gl_coord)
  outfname = [outfname,'_GLcoord'];
end
outfname = [outfname,'.mat'];
outfname = [expname,outfname];

%%% Latent heat of freezing
Lf = 3.34e5;

%%% Load pre-computed fluxes
load(fullfile('products',outfname),'thflux_stand','thflux_fluc','thflux_eddy','tflux_tavg','SHIfwFlx_tavg','eta','ETA');
thflux_mean_plot = thflux_stand;
thflux_fluc_plot = thflux_fluc;
thflux_eddy_plot = mean(thflux_eddy,2);
thflux_tot_plot = thflux_stand+thflux_fluc+mean(thflux_eddy,2);

surfQflux = tflux_tavg;
surfQflux(hFacC(:,:,1)==0) = SHIfwFlx_tavg(hFacC(:,:,1)==0) * Lf;
surfQflux(sum(hFacC,3)==0) = 0;

%%% Load positive and negative heatfunction components
outfname = [expname,'_PosNegHeatFunction'];
outfname = [outfname,'.mat'];
load(fullfile('products',outfname),'psiT_pos_mean','psiT_neg_mean');

%%% Depth-integrated heat fluxes due to positive and negative heat function
%%% components
thflux_pos = mean(psiT_pos_mean(:,1,:),3);
thflux_neg = mean(psiT_neg_mean(:,1,:),3);



%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0 0.45 .5 .5];
axpos(2,:) = [0.5 0.45 .5 .5];
axpos(3,:) = [0.07 0.06 .41 .39];
axpos(4,:) = [0.56 0.06 .41 .39];
cbpos = zeros(2,4);
cbpos(1,:) = [0.42 0.56 0.01 .1];
cbpos(2,:) = [0.92 0.56 0.01 .1];
axlabels = {'(a)','(b)','(c)','(d)'};
rho0 = 1027;
Cp = 4000;
colororder = get(gca,'ColorOrder');
linewidth = 1.5;

%%% Plotting range for salinity figure
latMin_b = min(min(YC));
latMax_b = YC(1,end-spongethickness);
lonMin_b = min(min(XC));
lonMax_b = XC(end-spongethickness,1);







%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(202)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34  1000  926]);







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
subplot('Position',axpos(1,:)+[.1 0 0 0]);
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
surfQflux_plot = surfQflux;
surfQflux_plot(sum(hFacC,3)==0) = NaN;
pcolorm(YC,XC,surfQflux_plot);
shading interp
caxis(clim);
colormap(gca,cmap);


%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',cbpos(1,:))
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
set(gca,'Position',axpos(1,:));
hold off

% Create textbox
annotation(gcf,'textbox',...
  [0.163 0.949863930885531 0.2685 0.0242980561555075],...
  'String',{'Downward surface heat flux'},'EdgeColor','None','FontSize',fontsize+2,'interpreter','latex');





%%%%%%%%%%%%%%%%%%%%%%
%%% ETA COORDINATE %%%
%%%%%%%%%%%%%%%%%%%%%%



%%% Plotting options
clim = [-9 11];
% cmap = cmocean('balance',50);
cmap = flip(haxby(20),1);
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

%%% Plot MOC coordinate
ETA_plot = ETA;
ETA_plot(sum(hFacC,3)==0) = NaN;
pcolorm(YC,XC,ETA_plot);
shading interp
caxis(clim);
colormap(gca,cmap);


%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',cbpos(2,:))
tightmap;

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

% Create title
annotation(gcf,'textbox',...
  [0.678 0.952023758099355 0.3005 0.0242980561555075],...
  'String',{'Pseudo-latitudinal coordinate, $\eta$'},'EdgeColor','None','FontSize',fontsize+2,'interpreter','latex');






for j=1:length(eta)
  surfQint(j) = sum(sum(surfQflux(ETA<eta(j)).*RAC(ETA<eta(j))));
end




%%% Heat flux
axes('Position',axpos(3,:));
plot(eta,-rho0*Cp*thflux_mean_plot/1e12,'Color',colororder(6,:),'LineWidth',linewidth);
hold on;
plot(eta,-rho0*Cp*thflux_fluc_plot/1e12,'Color',colororder(3,:),'LineWidth',linewidth);
plot(eta,-rho0*Cp*thflux_eddy_plot/1e12,'Color',colororder(5,:),'LineWidth',linewidth);
plot(eta,-rho0*Cp*thflux_tot_plot/1e12,'Color',colororder(4,:),'LineWidth',linewidth);
% plot(eta,-surfQint/1e12,'Color',colororder(4,:),'LineWidth',linewidth);
% plot(eta,0*eta,'k--');
plot([4 4],[-5 15],'--','Color',[.3 .3 .3],'LineWidth',2);
plot([0 0],[-5 15],'--','Color',[.3 .3 .3],'LineWidth',2);
hold off;
ylabel('Shoreward heat flux (TW)');
xlabel('\eta');
axis([-9 5 -3.5 7]);
set(gca,'FontSize',fontsize);
leghandle = legend('Mean','Seasonal/interannual','Eddy','Total','Location','NorthWest');
set(leghandle,'FontSize',fontsize);
text(-4.5,-3,'FRIS','FontSize',fontsize);
text(1.3,-3,'Shelf','FontSize',fontsize);
grid on;
box off;

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Top');
set(ax2,'YAxisLocation','Right');
set(ax2,'YLim',get(ax1,'YLim'));
set(ax2,'YTick',[]);
box off;
set(ax2,'XLim',get(ax1,'XLim')-78.16);
set(ax2,'FontSize',fontsize);
set(get(ax2,'XLabel'),'String','Approx. latitude');

%%% Positive/negative heat flux
axes('Position',axpos(4,:));
plot(eta,-rho0*Cp*thflux_tot_plot/1e12,'Color',colororder(4,:),'LineWidth',linewidth);
hold on;
plot(eta,-rho0*Cp*thflux_pos/1e12,'Color',colororder(2,:),'LineWidth',linewidth);
plot(eta,-rho0*Cp*thflux_neg/1e12,'Color',colororder(1,:),'LineWidth',linewidth);
% plot(eta,-rho0*Cp*(thflux_pos+thflux_neg)/1e12,'--','Color',colororder(4,:),'LineWidth',linewidth);
% plot(eta,0*eta,'k--')
plot([4 4],[-5 15],'--','Color',[.3 .3 .3],'LineWidth',2);
plot([0 0],[-5 15],'--','Color',[.3 .3 .3],'LineWidth',2);
hold off;
% ylabel('Heat flux (TW)');
xlabel('\eta');
axis([-9 4 -3.5 7]);
set(gca,'FontSize',fontsize);
leghandle = legend('Total','``Warm'''' component','``Cold'''' component','Location','NorthWest');
set(leghandle,'FontSize',fontsize);
text(-4.5,-3,'FRIS','FontSize',fontsize);
text(1.3,-3,'Shelf','FontSize',fontsize);
grid on;
box off;

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Top');
set(ax2,'YAxisLocation','Right');
set(ax2,'YLim',get(ax1,'YLim'));
set(ax2,'YTick',[]);
box off;
set(ax2,'XLim',get(ax1,'XLim')-78.16);
set(ax2,'FontSize',fontsize);
set(get(ax2,'XLabel'),'String','Approx. latitude');



%%% Add panel labels
annotation('textbox',[axpos(1,1)+0.01 axpos(1,2)+0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)+0.01 axpos(2,2)+0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.06 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.06 axpos(4,2)-0.05 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
