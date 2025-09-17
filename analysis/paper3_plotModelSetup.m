
%%%
%%% NSF_plotModelSetup
%%%
%%% Plots the model domain showing an instantaneous salinity snapshot and surface vorticity/divergence plots.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% loadexp;

%%% Iteration number of model output to plot
% iter = 1949760;
iter = 1949760 - 50*86400/60;

%%% Define grid to extract data sections
startLat = -82.3;
endLat = -72;
startLon = -62;
endLon = -40;
Nsec = 2001;
dLat = (endLat-startLat)/(Nsec-1);
dLon = (endLon-startLon)/(Nsec-1);
secLats = startLat:dLat:endLat;
secLons = startLon:dLon:endLon;
% plotm(secLats,secLons,'k-');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD TEMPERATURE SLICE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% Load temperature snapshot
% T = rdmdsWrapper(fullfile(exppath,'results','THETA_12hourly'),iter);
% 
% %%% For interpolating to mid-depth
% kmax = ones(Nx,Ny);
% kmin = ones(Nx,Ny);
% for i=1:Nx
%   for j=1:Ny
%     idx = find(squeeze(hFacC(i,j,:))>0);
%     if (~isempty(idx))
%       kmin(i,j) = min(idx);
%       kmax(i,j) = max(idx);
%     end
%   end
% end
% kn = ones(Nx,Ny);
% kp= ones(Nx,Ny);
% wn = 0.5*ones(Nx,Ny);
% wp = 0.5*ones(Nx,Ny);
% for i=1:Nx
%   for j=1:Ny
%     if (sum(hFacC(i,j,:),3)==0)
%       continue;
%     end
%     zmid = 0.5 * (SHELFICEtopo(i,j) + bathy(i,j));
%     kmid = max(find(squeeze(zz)>zmid));
%     if (isempty(kmid))
%       continue;
%     end
%     kp(i,j) = kmid;
%     kn(i,j) = kp(i,j) + 1;
%     wp(i,j) = (zmid-zz(kn(i,j))) / (zz(kp(i,j))-zz(kn(i,j)));
%     wn(i,j) = 1 - wp(i,j);
%   end
% end
% 
% %%% Interpolate to mid-depth
% FF = zeros(Nx,Ny);
% for i=1:Nx
%   for j=1:Ny
%     FF(i,j) = wp(i,j)*T(i,j,kp(i,j)) + wn(i,j)*T(i,j,kn(i,j));
%   end
% end
% FF(sum(hFacC,3)==0) = NaN;
% 
% %%% Free up memory
% clear('T');

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
figure(201)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34  1000  926]);
% fontsize = 18;
fontsize = 12;
axlabels = {'\textbf{A}','\textbf{B}','\textbf{C}'};
axpos = zeros(3,4);




%%%%%%%%%%%%%%%%%
%%% MAP PANEL %%%
%%%%%%%%%%%%%%%%%

%%% Plotting options
clim = [-5000 0];
axpos(1,:) = [0 0.5 .3 .3];

latMin = -90;
latMax = -60;
lonMin = -180;
lonMax = 180;

%%% Set up map plot
subplot('Position',axpos(1,:));
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
set(gca,'Position',axpos(1,:));
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








%%%%%%%%%%%%%%%%%%%%%%%
%%% TEMPERATURE MAP %%%
%%%%%%%%%%%%%%%%%%%%%%%



%%% Plotting options
clim = [-2.6 -1];
cmap = cmocean('balance',32);
cmap = cmap(8:23,:);
% cmap = cmocean('thermal',16);
% cmap = cmocean('amp',16);
% cmap = pmkmp(16,'swtth');
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
axpos(2,:) = [0.2 0.4 .85 .55];


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

%%% Plot temperature in color contours
pcolorm(YC,XC,FF);
shading interp
caxis(clim);
colormap(gca,cmap);


%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.9 0.45 0.01 .15])
tightmap;
title(h,'$^\circ{}C$','Fontsize',14,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.05 .05 .05],'BackgroundColor','none','Edgecolor','none')       
% Nlon = 101;
% dLon = (lonMax_c-lonMin_c)/(Nlon-1);
% plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'r--','LineWidth',2);
hold off;

%%% Plot transect location for VHF figure
hold on;
plotm(secLats,secLons,'k:','LineWidth',2);
hold off;

%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos(2,:));
hold off




% Create textbox
annotation(gcf,'textbox',...
  [0.769174442190669 0.857932603573533 0.0588235294117647 0.025568181818182],...
  'String',{'Transect'},...
  'FitBoxToText','off',...
  'EdgeColor','none','Color','k');

















%%% Iteration number of model output to plot
% iter = 1949760 - 316800;
% iter = 1949760 - 100*86400/60;
iter = 1949760 - 50*86400/60;
% iter = 1949760 - 250*86400/60;

% %%% Extract data along defined sections
% S = rdmdsWrapper(fullfile(exppath,'results','SALT_12hourly'),iter);
% T = rdmdsWrapper(fullfile(exppath,'results','THETA_12hourly'),iter);
% S(hFacC==0) = NaN;
% T(hFacC==0) = NaN;
% secS = zeros(Nsec,Nr);
% secT = zeros(Nsec,Nr);
% secB = zeros(Nsec,1);
% secI = zeros(Nsec,1);
% secH = zeros(Nsec,Nr);
% for n=1:Nsec
%   jm = find(secLats(n)>=YC(1,:),1,'last');
%   im = find(secLons(n)>=XC(:,1),1,'last');
%   jp = jm+1;
%   ip = im+1;
%   wp_y = (secLats(n)-YC(1,jm))/(YC(1,jp)-YC(1,jm));
%   wm_y = 1-wp_y;
%   wp_x = (secLons(n)-XC(im,1))/(XC(ip,1)-XC(im,1));
%   wm_x = 1-wp_x;
%   if (wp_x > wm_x)
%     xidx = ip;
%   else
%     xidx = im;
%   end
%   if (wp_y > wm_y)
%     yidx = jp;
%   else
%     yidx = jm;
%   end
%   secS(n,:) = squeeze(S(xidx,yidx,:));
%   secT(n,:) = squeeze(T(xidx,yidx,:));  
%   secH(n,:) = squeeze(hFacC(xidx,yidx,:));
%   secB(n) = bathy(xidx,yidx);
%   secI(n) = SHELFICEtopo(xidx,yidx);
% end











                
%%% Plotting options                
fontsize = 14;
labelspacing = 200;
axpos(3,:) = [0.06 0.07 0.87 0.3];
cb_pos = [0.95 0.07 0.015 0.3];
axlabel = '(c)';
icecolor = [186 242 239]/255;



%%%%%%%%%%%%%%%%%%
%%% TEMP SLICE %%%
%%%%%%%%%%%%%%%%%%

[ZZ,LA] = meshgrid(RC,secLats);
for n=1:Nsec
  kmin = find(secH(n,:)>0,1,'first');
  kmax = find(secH(n,:)>0,1,'last');
  if (isempty(kmin) || isempty(kmax))
    continue;
  end
  if (kmin == kmax)
    secT(n,kmin) = NaN;
    secS(n,kmin) = NaN;
    continue;
  end
  ZZ(n,kmin) = RF(kmin);% - (1-sechFac(n,kmin))*DRF(kmin);
  ZZ(n,kmax) = RF(kmax+1);% - (1-sechFac(n,kmax))*DRF(kmax);
  
end


axes('Position',axpos(3,:));
pcolor(LA,-ZZ,secT);
set(gca,'YDir','reverse');
shading interp;
hold on
[C,h] = contour(LA,-ZZ,secS,[34.1:.1:35],'EdgeColor','k');
clabel(C,h,'FontSize',12);
plot(secLats,-secB,'k-','LineWidth',3);
Iidx = find(secI==0,1,'first');
plot(secLats(1:Iidx),-secI(1:Iidx),'k-','LineWidth',3);
hold off;
set(gca,'YLim',[0 1400]);
set(gca,'XLim',[startLat (-73)]);
caxis(clim);
colormap(gca,cmap);
cbhandle = colorbar;
set(cbhandle,'Position',cb_pos);
title(cbhandle,'$^\circ$C','Fontsize',14,'interpreter','latex');
set(gca,'FontSize',fontsize);
xlabel('Latitude')
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);
% title(['Instantaneous potential temperature, salinity, ',datestr(datenum('01-Jan-2008')+iter*deltaT/86400)]);

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
msk = ones(Nsec,Nr);
msk(ZZ<repmat(secI,[1 Nr])) = NaN;
pcolor(LA,-ZZ,msk);
shading interp;
colormap(ax2,icecolor);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
set(ax2,'YLim',[0 1400]);
set(ax2,'YDir','reverse');
set(ax2,'XLim',[startLat (-73)]);
set(ax2,'Color','None')





%%% Add panel labels
annotation('textbox',[axpos(1,1)+0.01 axpos(1,2)-0.02 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)+0.15 axpos(2,2)+0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.06 axpos(3,2)-0.06 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');



%%% Add annotations
figure1 = gcf;

% Create line
annotation(figure1,'line',[0.916 0.868],...
  [0.117710583153348 0.130669546436285]);

% Create line
annotation(figure1,'line',[0.188 0.174],...
  [0.190064794816415 0.215421166306696]);

% Create line
annotation(figure1,'line',[0.555 0.586],...
  [0.312095032397408 0.269978401727862]);

% Create textbox
annotation(figure1,'textbox',...
  [0.755000000000001 0.126969762419007 0.13 0.0242980561555076],...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'Circumpolar Deep Water'});

% Create textbox
annotation(figure1,'textbox',...
  [0.546000000000001 0.243600431965443 0.1295 0.0242980561555076],...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'High Salinity Shelf Water'});

% Create textbox
annotation(figure1,'textbox',...
  [0.138 0.236041036717063 0.087 0.0242980561555076],...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'Ice Shelf Water'});