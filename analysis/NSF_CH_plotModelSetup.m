
%%%
%%% NSF_plotModelSetup
%%%
%%% Plots the model domain showing an instantaneous salinity snapshot and surface vorticity/divergence plots.
%%%

%%% Options
gendir = '/Volumes/Stewart-RAID1-B/UCLA/Julia/MITgcm_WS';
basedir = fullfile(gendir,'experiments');
expdir = basedir;
expname = 'hires_nest_onethirtysecond_notides_RTOPO2';
loadexp;



%%% Coriolis parameter
Omega = 2*pi*366/365/86400;
ff = 2*Omega*sind(YG);
rho_i = 920; %%% Density of ice

zlev = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load velocity snapshot %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Iteration number of model output to plot
iter = 67*t1day/deltaT;

U = rdmdsWrapper(fullfile(exppath,'results','UVEL_inst'),iter);
V = rdmdsWrapper(fullfile(exppath,'results','VVEL_inst'),iter);
U = U(:,:,zlev);
V = V(:,:,zlev); %%% Just need surface velocity
U(hFacW(:,:,zlev)==0) = NaN;
V(hFacS(:,:,zlev)==0) = NaN;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD VELOCITY AND MELT SLICES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter1 = 60*t1day/deltaT;
iter2 = 66*t1day/deltaT;
iter1 = 85*t1day/deltaT;
iter2 = 87*t1day/deltaT;

U1 = rdmdsWrapper(fullfile(exppath,'results','UVEL_inst'),iter1);
V1 = rdmdsWrapper(fullfile(exppath,'results','VVEL_inst'),iter1);
melt1 = rdmdsWrapper(fullfile(exppath,'results','SHIfwFlx_inst'),iter1);
U1 = U1(:,:,zlev);
V1 = V1(:,:,zlev); %%% Just need surface velocity
U1(hFacW(:,:,zlev)==0) = NaN;
V1(hFacS(:,:,zlev)==0) = NaN;

U2 = rdmdsWrapper(fullfile(exppath,'results','UVEL_inst'),iter2);
V2 = rdmdsWrapper(fullfile(exppath,'results','VVEL_inst'),iter2);
melt2 = rdmdsWrapper(fullfile(exppath,'results','SHIfwFlx_inst'),iter2);
U2 = U2(:,:,zlev);
V2 = V2(:,:,zlev); %%% Just need surface velocity
U2(hFacW(:,:,zlev)==0) = NaN;
V2(hFacS(:,:,zlev)==0) = NaN;

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


%%% Plotting range for speed figure
latMin_o = -83.5;
latMax_o = -64;
lonMin_o = -83;
lonMax_o = 20;

% %%% Plotting range for salinity figure
% latMin_o = -75;
% latMax_o = -69;
% lonMin_o = -32;
% lonMax_o = -12;

%%% Plotting range for salinity figure
latMin_b = min(min(YC));
latMax_b = YC(1,end);
lonMin_b = min(min(XC));
lonMax_b = XC(end,1);

%%% Plotting ranges for vorticity/divergence figures
latMin_c = YC(1,1+spongethickness);
latMax_c = -72;
lonMin_c = -28;
lonMax_c = -18;

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
set(gcf,'Position',[417    34  1000  1000]);
% fontsize = 18;
fontsize = 12;






%%%%%%%%%%%%%%%%%
%%% MAP PANEL %%%
%%%%%%%%%%%%%%%%%

%%% Plotting options
clim = [-5000 0];
axpos = [0.1 0.6 .8 .4];

latMin = -90;
latMax = -60;
lonMin = -180;
lonMax = 180;

%%% Set up map plot
subplot('Position',axpos);
axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_o-0 latMax_o+0], ...
  'MapLonLimit',[lonMin_o-0 lonMax_o+0], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on'); 
axis off;
% setm(gca,'MLabelParallel',-20)


RTOPO_bathy_plot = RTOPO_bathymetry;
RTOPO_bathy_plot(RTOPO_ice <= RTOPO_bathymetry) = NaN;
RTOPO_elev_plot = RTOPO_elev;
RTOPO_elev_plot(RTOPO_ice <= RTOPO_bathymetry) = NaN;
RTOPO_elev_plot(RTOPO_elev_plot==0) = NaN;

%%% Plot bathymetry in color contours
[LA,LO] = meshgrid(RTOPO_latitude,RTOPO_longitude);
pcolorm(double(LA),double(LO),double(RTOPO_bathy_plot));
shading interp
colormap(gca,haxby(50));
caxis(clim);
set(gca,'Position',axpos);
% hold on;
ax1 = gca;
axpos = get(ax1,'Position');

h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.8 0.65 0.01 .1]);
title(h,'m','Fontsize',14,'interpreter','latex');




%%% Absolute surface horizontal flow speed for plotting
uabs_plot = sqrt(U.^2 + V.^2);

subplot('Position',axpos+[0 0.5 0 0]);
ax3 = axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_o-0 latMax_o+0], ...
  'MapLonLimit',[lonMin_o-0 lonMax_o+0], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');  
axis off;
set(ax3,'Position',axpos);
pcolorm(YC,XC,uabs_plot);
shading flat;
colormap(ax3,cmocean('amp'));
caxis([0 1]);

h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.6 0.7 0.01 .1]);
title(h,'m/s','Fontsize',14,'interpreter','latex');








subplot('Position',axpos+[0 0.5 0 0]);
ax2 = axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_o-0 latMax_o+0], ...
  'MapLonLimit',[lonMin_o-0 lonMax_o+0], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');  
axis off;
set(ax2,'Position',axpos);
pcolorm(double(LA),double(LO),double(RTOPO_elev_plot));
hold on;
Nlon = 101;
dLon = (lonMax_b-lonMin_b)/(Nlon-1);
plotm([latMin_b*ones(1,Nlon) latMax_b*ones(1,Nlon) latMin_b],[lonMin_b:dLon:lonMax_b lonMax_b:-dLon:lonMin_b lonMin_b],'k--','LineWidth',1);
Nlon = 101;
dLon = (lonMax_o-lonMin_o)/(Nlon-1);
plotm([latMin_o*ones(1,Nlon) latMax_o*ones(1,Nlon) latMin_o],[lonMin_o:dLon:lonMax_o lonMax_o:-dLon:lonMin_o lonMin_o],'k-','LineWidth',2);
hold off;
shading interp;
% colormap(ax2,[220,243,255]/255);
colormap(ax2,icecolor);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRST VORTICITY PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Compute vorticity
vort1 = zeros(Nx,Ny);
vort1(:,2:Ny) = - (U1(:,2:Ny)-U1(:,1:Ny-1))./DYC(:,2:Ny);
vort1(2:Nx,:) = vort1(2:Nx,:) + (V1(2:Nx,:)-V1(1:Nx-1,:))./DXC(2:Nx,:);

%%% Plotting options
clim = [-1 1];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
axpos = [0.05 0.25 .4 .3];

%%% Set up map plot
subplot('Position',axpos);
axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_c latMax_c], ...
  'MapLonLimit',[lonMin_c lonMax_c], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 1, ...
  'MLineLocation', 2,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)

%%% Plot vorticity in color contours
pcolorm(YG,XG,vort1./abs(ff));
shading interp
colormap(gca,cmocean('delta',40));
caxis(clim);

%%% Add colorbar and title
h = colorbar;
set(h,'FontSize',12);
set(h,'Position',[0.5 0.43 0.01 .12]);
tightmap;
title(h,'$\zeta/|f|$','Fontsize',14,'interpreter','latex');

%%% Add bathymetry contours
hold on;
thick_plot = SHELFICEtopo-bathy;
[cs,C] = contourm(YC,XC,thick_plot,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')     
hold off;

hold on;
Nlon = 101;
dLon = (lonMax_c-lonMin_c)/(Nlon-1);
plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'k-','LineWidth',1);
hold off
        
%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos);
ax = gca;

%%% Add plot of ice shelf melt rate
subplot('Position',[0 0 0.01 0.01])
ax2 = axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_c latMax_c], ...
  'MapLonLimit',[lonMin_c lonMax_c], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 0, ...
  'MLineLocation', 0,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');   
axis off;
setm(ax2,'MLabelParallel',-20)

melt1(melt1==0) = NaN;
pcolorm(YC,XC,-melt1/rho_i*t1year);
shading flat
colormap(ax2,cmocean('balance'));
hold on;
Nlon = 101;
dLon = (lonMax_c-lonMin_c)/(Nlon-1);
plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'k-','LineWidth',1);
hold off
drawnow;
set(ax2,'Position',axpos+[-0.04 -0.013 0.08 0.067]);
caxis([-10 10]);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SECOND VORTICITY PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Compute vorticity
vort2 = zeros(Nx,Ny);
vort2(:,2:Ny) = - (U2(:,2:Ny)-U2(:,1:Ny-1))./DYC(:,2:Ny);
vort2(2:Nx,:) = vort2(2:Nx,:) + (V2(2:Nx,:)-V2(1:Nx-1,:))./DXC(2:Nx,:);

%%% Plotting options
clim = [-1 1];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
axpos = [0.55 0.25 .4 .3];

%%% Set up map plot
subplot('Position',axpos);
axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_c latMax_c], ...
  'MapLonLimit',[lonMin_c lonMax_c], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 1, ...
  'MLineLocation', 2,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)

%%% Plot vorticity in color contours
pcolorm(YG,XG,vort2./abs(ff));
shading interp
colormap(gca,cmocean('delta',40));
caxis(clim);

%%% Add bathymetry contours
hold on;
thick_plot = SHELFICEtopo-bathy;
[cs,C] = contourm(YC,XC,thick_plot,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')     
hold off;

hold on;
Nlon = 101;
dLon = (lonMax_c-lonMin_c)/(Nlon-1);
plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'k-','LineWidth',1);
hold off

%%% Add colorbar and title (UNUSED)
h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[1.1 0.4 0.01 .15]);
tightmap;
title(h,'$\zeta/|f|$','Fontsize',14,'interpreter','latex');
        
%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos);
ax = gca;

%%% Add plot of ice shelf melt rate
subplot('Position',[0 0 0.01 0.01])
ax2 = axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_c latMax_c], ...
  'MapLonLimit',[lonMin_c lonMax_c], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 0, ...
  'MLineLocation', 0,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');   
axis off;
setm(ax2,'MLabelParallel',-20)

melt2(melt2==0) = NaN;
pcolorm(YC,XC,-(melt2)/rho_i*t1year);
shading flat
colormap(ax2,cmocean('balance'));
hold on;
Nlon = 101;
dLon = (lonMax_c-lonMin_c)/(Nlon-1);
plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'k-','LineWidth',1);
hold off


%%% Add colorbar and title
h = colorbar;
set(h,'FontSize',12);
set(h,'Position',[0.5 0.25 0.01 .12]);
tightmap;
title(h,'m/yr','Fontsize',14,'interpreter','latex');



% hold on;
% [cs,C] = contourm(YC,XC,SHELFICEtopo,[-100 -100],'EdgeColor','k'); 
% hold off;


drawnow;
set(ax2,'Position',axpos);%+[-0.04 -0.013 0.08 0.067]);
caxis([-10 10]);


load('EddyDetections/April 1 centers.mat');
hold on;
plotm(centers4.y1,centers4.x1,'x','Color','g');
hold off;




%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time series plot %%%
%%%%%%%%%%%%%%%%%%%%%%%%

load('./products/SWIT_timeseries.mat');

axpos = [0.1 0.05 .8 .15];

%%% Set up map plot
subplot('Position',axpos);
ax = plotyy(datevec,SWIT_melt/1e12*t1year,datevec,rhoConst*SWIT_KE/1e12);
datetick('x');
set(ax(1),'FontSize',13);
set(ax(2),'FontSize',13);
xlabel(ax(1),'2008');
ylabel(ax(1),'SWIT melt (Gt/yr)');
ylabel(ax(2),'SWIT KE (TJ)');
set(ax(1),'ylim',[0 30]);
set(ax(1),'YTick',[0 10 20 30]);
set(ax(2),'YTick',[0 25 50]);



% Create textbox
annotation(gcf,'textbox',[0.18 0.589524838012959 0.15 0.03],...
  'String','Mar 25, 2008',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',[0.7 0.589524838012959 0.15 0.03],...
  'String','Mar 27, 2008',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');


% Create textbox
annotation(gcf,'textbox',[0.2 0.689524838012959 0.03 0.03],...
  'String','(a)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.04 0.229524838012959 0.03 0.03],'String','(b)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.55 0.229524838012959 0.03 0.03],...
  'String','(c)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.04 0.03 0.03 0.01],'String','(d)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');