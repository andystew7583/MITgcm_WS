
%%%
%%% NSF_plotModelSetup
%%%
%%% Plots the model domain showing an instantaneous salinity snapshot and surface vorticity/divergence plots.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;

%%% Iteration number of model output to plot
iter = 1949760;

%%% Coriolis parameter
Omega = 2*pi*366/365/86400;
ff = 2*Omega*sind(YG);


%%% Load velocity snapshot
U = rdmdsWrapper(fullfile(exppath,'results','UVEL_12hourly'),iter);
V = rdmdsWrapper(fullfile(exppath,'results','VVEL_12hourly'),iter);
U = U(:,:,1);
V = V(:,:,1); %%% Just need surface velocity
zlev = 1;
U(hFacW(:,:,zlev)==0) = NaN;
V(hFacS(:,:,zlev)==0) = NaN;




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD SALINITY SLICE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load salinity snapshot
S = rdmdsWrapper(fullfile(exppath,'results','SALT_12hourly'),iter);

%%% For interpolating to mid-depth
kmax = ones(Nx,Ny);
kmin = ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    idx = find(squeeze(hFacC(i,j,:))>0);
    if (~isempty(idx))
      kmin(i,j) = min(idx);
      kmax(i,j) = max(idx);
    end
  end
end
kn = ones(Nx,Ny);
kp= ones(Nx,Ny);
wn = 0.5*ones(Nx,Ny);
wp = 0.5*ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    if (sum(hFacC(i,j,:),3)==0)
      continue;
    end
    zmid = 0.5 * (SHELFICEtopo(i,j) + bathy(i,j));
    kmid = max(find(squeeze(zz)>zmid));
    if (isempty(kmid))
      continue;
    end
    kp(i,j) = kmid;
    kn(i,j) = kp(i,j) + 1;
    wp(i,j) = (zmid-zz(kn(i,j))) / (zz(kp(i,j))-zz(kn(i,j)));
    wn(i,j) = 1 - wp(i,j);
  end
end

%%% Interpolate to mid-depth
FF = zeros(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    FF(i,j) = wp(i,j)*S(i,j,kp(i,j)) + wn(i,j)*S(i,j,kn(i,j));
  end
end
FF(sum(hFacC,3)==0) = NaN;

%%% Free up memory
clear('S');

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
set(gcf,'Position',[417    34  1000  926]);
% fontsize = 18;
fontsize = 12;






%%%%%%%%%%%%%%%%%
%%% MAP PANEL %%%
%%%%%%%%%%%%%%%%%

%%% Plotting options
clim = [-5000 0];
axpos = [0 0.5 .3 .3];

latMin = -90;
latMax = -60;
lonMin = -180;
lonMax = 180;

%%% Set up map plot
subplot('Position',axpos);
axesm('stereo',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ...   
  'PLineLocation', 10, ...
  'MLineLocation', 60,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');  
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
set(h,'Position',[0.28 0.5 0.01 .1]);
title(h,'m','Fontsize',14,'interpreter','latex');

subplot('Position',axpos+[0 0.5 0 0]);
ax2 = axesm('stereo',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ...   
  'PLineLocation', 10, ...
  'MLineLocation', 60,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');  
axis off;
set(ax2,'Position',axpos);
pcolorm(double(LA),double(LO),double(RTOPO_elev_plot));
hold on;
Nlon = 101;
dLon = (lonMax_b-lonMin_b)/(Nlon-1);
plotm([latMin_b*ones(1,Nlon) latMax_b*ones(1,Nlon) latMin_b],[lonMin_b:dLon:lonMax_b lonMax_b:-dLon:lonMin_b lonMin_b],'k-','LineWidth',3);
Nlon = 101;
dLon = (lonMax_o-lonMin_o)/(Nlon-1);
plotm([latMin_o*ones(1,Nlon) latMax_o*ones(1,Nlon) latMin_o],[lonMin_o:dLon:lonMax_o lonMax_o:-dLon:lonMin_o lonMin_o],'k--','LineWidth',2);
hold off;
shading interp;
% colormap(ax2,[220,243,255]/255);
colormap(ax2,icecolor);




%%%%%%%%%%%%%%%%%%%%%%
%%% SALINITY PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%



%%% Plotting options
clim = [34.2 35];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
axpos = [0.2 0.4 .85 .55];

%%% Set up map plot
subplot('Position',axpos+[.1 0 0 0]);
axesm('eqaconicstd',...
  'fontsize',13,...
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
pcolorm(YC,XC,FF);
shading interp
colormap(gca,cmocean('haline',16));
caxis(clim);

%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.9 0.45 0.01 .15])
tightmap;
title(h,'g/kg','Fontsize',14,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')       
Nlon = 101;
dLon = (lonMax_c-lonMin_c)/(Nlon-1);
plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'r--','LineWidth',2);
hold off;

%%% Plot transect location for VHF figure
startLat = -82.3;
endLat = -72;
startLon = -62;
endLon = -40;
Nsec = 2001;
dLat = (endLat-startLat)/(Nsec-1);
dLon = (endLon-startLon)/(Nsec-1);
secLats = startLat:dLat:endLat;
secLons = startLon:dLon:endLon;
hold on;
plotm(secLats,secLons,'k:','LineWidth',2);
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
vort(:,2:Ny) = - (U(:,2:Ny,zlev)-U(:,1:Ny-1,zlev))./DYC(:,2:Ny);
vort(2:Nx,:) = vort(2:Nx,:) + (V(2:Nx,:,zlev)-V(1:Nx-1,:,zlev))./DXC(2:Nx,:);

%%% Plotting options
clim = [-1 1];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
axpos = [0.05 0 .4 .4];

%%% Set up map plot
subplot('Position',axpos);
axesm('eqaconicstd',...
  'fontsize',13,...
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
set(gca,'FontSize',10);
set(h,'Position',[0.45 0.05 0.01 .15]);
tightmap;
title(h,'$\zeta/|f|$','Fontsize',14,'interpreter','latex');

%%% Add bathymetry contours
hold on;
thick_plot = SHELFICEtopo-bathy;
% thick_plot((SHELFICEtopo<0) | (bathy==SHELFICEtopo)) = NaN;
[cs,C] = contourm(YC,XC,thick_plot,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')       
        
%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos);
hold off

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
div(:,1:Ny-1) = (V(:,2:Ny,zlev)-V(:,1:Ny-1,zlev))./DYG(:,1:Ny-1);
div(1:Nx-1,:) = div(2:Nx,:) + (U(2:Nx,:,zlev)-U(1:Nx-1,:,zlev))./DXG(1:Nx-1,:);

%%% Plotting options
latMin = -78.5;
latMax = -74;
lonMin = -55;
lonMax = -30;
clim = [-.5 .5];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
axpos = [0.55 0 0.4 .4];

%%% Set up map plot
subplot('Position',axpos)
axesm('eqaconicstd',...
  'fontsize',13,...
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
set(gca,'FontSize',10);
set(h,'Position',[0.95 0.05 0.01 .15]);
tightmap;
title(h,'$\delta/|f|$','Fontsize',14,'interpreter','latex');

%%% Add bathymetry contours
hold on;
thick_plot = SHELFICEtopo-bathy;
% thick_plot((SHELFICEtopo<0) | (bathy==SHELFICEtopo)) = NaN;
[cs,C] = contourm(YC,XC,thick_plot,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')               
hold off
set(gca,'Position',axpos);
ax1 = gca;

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
annotation(gcf,'textbox',[0.017 0.489524838012959 0.05 0.05],...
  'String','(a)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.392 0.42149028077754 0.05 0.0500000000000001],'String','(b)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.0440000000000003 0.0413606911447106 0.05 0.0500000000000001],...
  'String','(c)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.539 0.0478401727861793 0.05 0.0500000000000001],'String','(d)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');