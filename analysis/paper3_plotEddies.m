
%%%
%%% paper3_plotEddies.m
%%%
%%% Plots the model domain showing an instantaneous barotropic vorticity snapshot and surface vorticity/divergence plots.
%%%

% %%% Options
% expdir = '../experiments';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% loadexp;
% 
% %%% Iteration number of model output to plot
% iter = 1949760;
% 
% %%% Coriolis parameter
% Omega = 2*pi*366/365/86400;
% ff = 2*Omega*sind(YG);
% 
% 
% %%% Load velocity snapshot
% U = rdmdsWrapper(fullfile(exppath,'results','UVEL_12hourly'),iter);
% V = rdmdsWrapper(fullfile(exppath,'results','VVEL_12hourly'),iter);
% Us = U(:,:,1);
% Vs = V(:,:,1); %%% Just need surface velocity
% zlev = 1;
% Us(hFacW(:,:,zlev)==0) = NaN;
% Vs(hFacS(:,:,zlev)==0) = NaN;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD BATHYMETRY DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Load RTOPO data
datadir = '../data/RTOPO';
RTOPO_latitude = ncread(fullfile(datadir,'RTOPO2.nc'),'lat');
RTOPO_longitude = ncread(fullfile(datadir,'RTOPO2.nc'),'lon');
RTOPO_bathymetry = ncread(fullfile(datadir,'RTOPO2.nc'),'bedrock_topography');
RTOPO_ice = ncread(fullfile(datadir,'RTOPO2.nc'),'ice_base_topography');
RTOPO_elev = ncread(fullfile(datadir,'RTOPO2.nc'),'surface_elevation');
subfac = 10;
RTOPO_latitude = RTOPO_latitude(1:subfac:end);
RTOPO_longitude = RTOPO_longitude(1:subfac:end);
RTOPO_bathymetry = RTOPO_bathymetry(1:subfac:end,1:subfac:end);
RTOPO_ice = RTOPO_ice(1:subfac:end,1:subfac:end);
RTOPO_elev = RTOPO_elev(1:subfac:end,1:subfac:end);

%%% Plotting range for salinity figure
latMin_o = -83.5;
latMax_o = -64;
lonMin_o = -83;
lonMax_o = 20;

%%% Plotting range for salinity figure
latMin_b = min(min(YC));
latMax_b = YC(1,end-spongethickness);
lonMin_b = min(min(XC));
lonMax_b = XC(end-spongethickness,1);

%%% Plotting ranges for vorticity/divergence figures
latMin_c = -78.5;
latMax_c = -74;
lonMin_c = -55;
lonMax_c = -30;

%%% Color for ice shelves in plots
% icecolor = [7*16+3 9*16+11 14*16]/255;
icecolor = [186 242 239]/255;

%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(202)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[ 308         150        1261        1212]);
% fontsize = 18;
fontsize = 16;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEPTH-AVERAGED VORTICITY PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Barotropic vorticity
vort = zeros(Nx,Ny);
% U(hFacW==0) = NaN;
% V(hFacS==0) = NaN;

ubt = sum(U.*DRF.*hFacW,3) ./ sum(DRF.*hFacW,3);
vbt = sum(V.*DRF.*hFacS,3) ./ sum(DRF.*hFacS,3);
vort(:,2:Ny) = - (ubt(:,2:Ny)-ubt(:,1:Ny-1))./DYC(:,2:Ny);
vort = vort + (vbt([2:Nx 1],:)-vbt(:,:))./DXC; 

%%% Plotting options
clim = [-.5 .5];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
axpos = [0.05 0.39 .85 .575];

%%% Set up map plot
subplot('Position',axpos+[.1 0 0 0]);
axesm('eqaconicstd',...
  'fontsize',fontsize,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_b latMax_b], ...
  'MapLonLimit',[lonMin_b lonMax_b], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');  
axis off;
setm(gca,'MLabelParallel',-20)

%%% Plot salinity in color contours
pcolorm(YC,XC,vort./abs(ff));
shading interp
colormap(gca,cmocean('balance',40));
caxis(clim);

%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',[0.7 0.45 0.01 .15])
tightmap;
title(h,'$\zeta_{\mathrm{BT}}/|f|$','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',12,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')       
Nlon = 101;
dLon = (lonMax_c-lonMin_c)/(Nlon-1);
plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'k--','LineWidth',2);
hold off;

%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos);
hold off






%%%%%%%%%%%%%%%%%%%%%%%
%%% VORTICITY PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Compute vorticity
vort = zeros(Nx,Ny);
zlev = 1;
vort(:,2:Ny) = - (Us(:,2:Ny,zlev)-Us(:,1:Ny-1,zlev))./DYC(:,2:Ny);
vort(2:Nx,:) = vort(2:Nx,:) + (Vs(2:Nx,:,zlev)-Vs(1:Nx-1,:,zlev))./DXC(2:Nx,:);

%%% Plotting options
clim = [-1 1];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
axpos = [0.05 0 .4 .39];

%%% Set up map plot
subplot('Position',axpos);
axesm('eqaconicstd',...
  'fontsize',fontsize,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_c latMax_c], ...
  'MapLonLimit',[lonMin_c lonMax_c], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)

%%% Plot vorticity in color contours
pcolorm(YG,XG,vort./abs(ff));
shading interp
colormap(gca,cmocean('balance',40));
caxis(clim);

%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',[0.45 0.05 0.01 .15]);
tightmap;
title(h,'$\zeta_{\mathrm{surf}}/|f|$','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
thick_plot = SHELFICEtopo-bathy;
% thick_plot((SHELFICEtopo<0) | (bathy==SHELFICEtopo)) = NaN;
[cs,C] = contourm(YC,XC,thick_plot,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',12,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')       
        
%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos);
hold off

subplot('Position',[0 0 0.01 0.01])
ax2 = axesm('eqaconicstd',...
  'fontsize',fontsize,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_c latMax_c], ...
  'MapLonLimit',[lonMin_c lonMax_c], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 0, ...
  'MLineLocation', 0,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(ax2,'MLabelParallel',-20)


ice_plot = SHELFICEtopo;
ice_plot((ice_plot == 0) | (SHELFICEtopo==bathy)) = NaN;
pcolorm(YG,XG,ice_plot);
shading interp
colormap(ax2,icecolor);
drawnow;
set(ax2,'Position',axpos+[0 0.002 0 0]);


%%%%%%%%%%%%%%%%%%%%%%%%
%%% DIVERGENCE PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Compute horizontal divergence
div = zeros(Nx,Ny);
zlev = 1;
div(:,1:Ny-1) = (Vs(:,2:Ny,zlev)-Vs(:,1:Ny-1,zlev))./DYG(:,1:Ny-1);
div(1:Nx-1,:) = div(2:Nx,:) + (Us(2:Nx,:,zlev)-Us(1:Nx-1,:,zlev))./DXG(1:Nx-1,:);

%%% Plotting options
latMin = -78.5;
latMax = -74;
lonMin = -55;
lonMax = -30;
clim = [-.5 .5];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
axpos = [0.55 0 0.4 .39];

%%% Set up map plot
subplot('Position',axpos)
axesm('eqaconicstd',...
  'fontsize',fontsize,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_c latMax_c], ...
  'MapLonLimit',[lonMin_c lonMax_c], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)

%%% Plot vorticity in color contours
pcolorm(YC,XC,div./abs(ff));
shading interp
colormap(gca,cmocean('curl',40));
caxis(clim);

%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',[0.95 0.05 0.01 .15]);
tightmap;
title(h,'$\delta_{\mathrm{surf}}/|f|$','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
thick_plot = SHELFICEtopo-bathy;
% thick_plot((SHELFICEtopo<0) | (bathy==SHELFICEtopo)) = NaN;
[cs,C] = contourm(YC,XC,thick_plot,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',12,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')               
hold off
set(gca,'Position',axpos);
ax1 = gca;

subplot('Position',[0 0 0.01 0.01])
ax2 = axesm('eqaconicstd',...
  'fontsize',fontsize,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_c latMax_c], ...
  'MapLonLimit',[lonMin_c lonMax_c], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 0, ...
  'MLineLocation', 0,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(ax2,'MLabelParallel',-20)
set(ax2,'Position',axpos);

ice_plot = SHELFICEtopo;
ice_plot((ice_plot == 0) | (SHELFICEtopo==bathy)) = NaN;
pcolorm(YG,XG,ice_plot);
shading interp
colormap(ax2,icecolor);
drawnow;
set(ax2,'Position',axpos+[0 0.002 0 0]);



% Create textbox
annotation(gcf,'textbox',...
  [0.252 0.42149028077754 0.05 0.0500000000000001],'String','\textbf{A}',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',fontsize,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.0440000000000003 0.0413606911447106 0.05 0.0500000000000001],...
  'String','\textbf{B}',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',fontsize,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.539 0.0478401727861793 0.05 0.0500000000000001],'String','\textbf{C}',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',fontsize,...
  'FitBoxToText','off');