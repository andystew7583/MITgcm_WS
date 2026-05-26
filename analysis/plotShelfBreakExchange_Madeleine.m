%%%
%%% paper3_plotDepthIntegratedHeatFluxMaps.m
%%%
%%% Plots maps of depth-integrate heat fluxes, divided into positive and
%%% negative components.
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_RTOPO2';
loadexp;

compute_diags = false;

%%% Define reference density
k_pot_dens = 1;
PP = -RC(k_pot_dens)*gravity*rhonil/1e4*ones(Nx,Ny,Nr);
% salt_ref = 34.67;
% theta_ref = 0;
salt_ref = 34.37;
theta_ref = -1.8;
dens_ref = densjmd95(salt_ref,theta_ref,PP(1,1,1));

%%% Define grid to extract data sections 
stat_lat = [-71.0, -72.8149, -73.0962, -74.9461, -75.1044, -75.1718];
stat_lon = [-60, -52.4542, -48.8007, -39.8791, -31.9771, -27.0490];
Nsec = 401;
dLon = (stat_lon(end)-stat_lon(1))/(Nsec-1);
secLons = stat_lon(1):dLon:stat_lon(end);
secLats = interp1(stat_lon,stat_lat,secLons,'spline');
startLat = secLats(1);
endLat = secLats(end);
startLon = secLons(1);
endLon = secLons(end);

if (compute_diags)
  
  %%% Compute bottom density
  load(fullfile('products',[expname,'_MOC_PD0_layers.mat']),'dens_tavg');
  pot_dens_bot = zeros(Nx,Ny);
  for i=1:Nx
    for j=1:Ny
      kmax = find(squeeze(hFacC(i,j,:))>0,1,'last');
      if (isempty(kmax))
        pot_dens_bot(i,j) = NaN;
      else
        pot_dens_bot(i,j) = dens_tavg(i,j,kmax);
      end
    end
  end
  
  %%% Compute depth-integrated density fluxes
  load(fullfile('products',[expname,'_TSfluxes.mat']),'salt_tavg','theta_tavg','uvel_tavg','vvel_tavg');
  
   %%% Pre-compute indices for interpolation  
  secB = zeros(Nsec,1);
  secI = zeros(Nsec,1);
  secH = zeros(Nsec,Nr);
  xidx = zeros(Nsec,1);
  yidx = zeros(Nsec,1);
  for n=1:Nsec
    jm = find(secLats(n)>=YC(1,:),1,'last');
    im = find(secLons(n)>=XC(:,1),1,'last');
    jp = jm+1;
    ip = im+1;
    wp_y = (secLats(n)-YC(1,jm))/(YC(1,jp)-YC(1,jm));
    wm_y = 1-wp_y;
    wp_x = (secLons(n)-XC(im,1))/(XC(ip,1)-XC(im,1));
    wm_x = 1-wp_x;
    if (wp_x > wm_x) %%% Nearest-neighbor interpolation
      xidx(n) = ip;
    else
      xidx(n) = im;
    end
    if (wp_y > wm_y)
      yidx(n) = jp;
    else
      yidx(n) = jm;
    end
    secH(n,:) = squeeze(hFacC(xidx(n),yidx(n),:));
    secB(n) = bathy(xidx(n),yidx(n));
    secI(n) = SHELFICEtopo(xidx(n),yidx(n));

  end
  
  %%% Meshgrid for plotting
  Nz_sec = 600;
  secT = zeros(Nsec,Nz_sec); %%% To store sections of property of interest 
  secS = zeros(Nsec,Nz_sec);
  secD = zeros(Nsec,Nz_sec);
  secV = zeros(Nsec,Nz_sec);
  [ZZ_sec,LO_sec] = meshgrid(1:Nz_sec,secLons); %%% Placeholder
  for n=1:Nsec
    kmin = find(secH(n,:)>0,1,'first');
    kmax = find(secH(n,:)>0,1,'last');
    if (isempty(kmin) || isempty(kmax) || (kmin == kmax))
      continue;
    end
    dz = -secB(n)/Nz_sec;
    ZZ_sec(n,:) = -secI(n)+dz/2:dz:-secB(n);
    T_col = squeeze(theta_tavg(xidx(n),yidx(n),:));
    T_col(secH(n,:)==0) = NaN;
    secT(n,:) = interp1(squeeze(-RC(kmin:kmax)),T_col(kmin:kmax),ZZ_sec(n,:),'linear','extrap');
    S_col = squeeze(salt_tavg(xidx(n),yidx(n),:));
    S_col(secH(n,:)==0) = NaN;
    secS(n,:) = interp1(squeeze(-RC(kmin:kmax)),S_col(kmin:kmax),ZZ_sec(n,:),'linear','extrap');
    D_col = squeeze(dens_tavg(xidx(n),yidx(n),:));
    D_col(secH(n,:)==0) = NaN;
    secD(n,:) = interp1(squeeze(-RC(kmin:kmax)),D_col(kmin:kmax),ZZ_sec(n,:),'linear','extrap');
  
    U_col = squeeze(uvel_tavg(xidx(n),yidx(n),:)); 
    V_col = squeeze(vvel_tavg(xidx(n),yidx(n),:)); 
    if (n == 1)
      dx_sec = (secLons(n+1)-secLons(n))*cosd(secLats(n))
      dy_sec = (secLats(n+1)-secLats(n));
    elseif (n == Nsec)
      dx_sec = (secLons(n)-secLons(n-1))*cosd(secLats(n))
      dy_sec = (secLats(n)-secLats(n-1));
    else
      dx_sec = (secLons(n+1)-secLons(n-1))*cosd(secLats(n))
      dy_sec = (secLats(n+1)-secLats(n-1));
    end
    Vperp_col = (-U_col*dy_sec + V_col*dx_sec) / sqrt(dx_sec^2+dy_sec^2);
    Vperp_col(secH(n,:)==0) = NaN;
    secV(n,:) = interp1(squeeze(-RC(kmin:kmax)),Vperp_col(kmin:kmax),ZZ_sec(n,:),'linear','extrap');

  end
  
  %%% Free up memory
  clear('dens_tavg');
  
  %%% Compute depth-integrated density fluxes
  load(fullfile('products',[expname,'_TSfluxes.mat']),...
    'uvelth_tavg','uvelslt_tavg','salt_tavg','theta_tavg');
  
  salt_u = 0.5*(salt_tavg(:,:,:)+salt_tavg([Nx 1:Nx-1],:,:));
  theta_u = 0.5*(theta_tavg(:,:,:)+theta_tavg([Nx 1:Nx-1],:,:));
  [alpha_u,beta_u] = calcAlphaBeta(salt_u,theta_u,PP);
  dens_u = densjmd95(salt_u,theta_u,PP);
  
  uveldens = ...
    (dens_u).*uvel_tavg ...
    - rhonil.*alpha_u.*(uvelth_tavg-uvel_tavg.*theta_u) ...
    + rhonil.*beta_u.*(uvelslt_tavg-uvel_tavg.*salt_u);
  uveldens_zint = sum(uveldens.*DRF.*hFacW,3);
  uvel_zint = sum(uvel_tavg.*DRF.*hFacW,3);
  
  clear('uvel_tavg','uvelth_tavg','uvelslt_tavg','salt_u',...
    'theta_u','alpha_u','beta_u','dens_u','uveldens');
  
  load(fullfile('products',[expname,'_TSfluxes.mat']),...
    'vvelth_tavg','vvelslt_tavg');
  
  salt_v = 0.5*(salt_tavg(:,:,:)+salt_tavg(:,[Ny 1:Ny-1],:));
  theta_v = 0.5*(theta_tavg(:,:,:)+theta_tavg(:,[Ny 1:Ny-1],:));
  [alpha_v,beta_v] = calcAlphaBeta(salt_v,theta_v,PP);
  dens_v = densjmd95(salt_v,theta_v,PP);
  
  
  vveldens = ...
    (dens_v).*vvel_tavg ...
    - rhonil.*alpha_v.*(vvelth_tavg-vvel_tavg.*theta_v) ...
    + rhonil.*beta_v.*(vvelslt_tavg-vvel_tavg.*salt_v);
  vveldens_zint = sum(vveldens.*DRF.*hFacS,3);
  vvel_zint = sum(vvel_tavg.*DRF.*hFacS,3);
  
  clear('vvel_tavg','vvelth_tavg','vvelslt_tavg','salt_v',...
    'theta_v','alpha_v','beta_v','dens_v','vveldens','theta_tavg','salt_tavg');

end











%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(220)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[837          80        1174        1212]);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DENSITY/DENSITY FLUX PLOT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plotting options
% latMin = min(min(YC));
latMin_b = -83.5;
% latMax = YC(1,end-spongethickness);
latMax_b = -68;
% lonMin = min(min(XC));
lonMin_b = -82;
lonMax_b = -10;
fontsize = 18;
axpos = zeros(4,4);
axpos(1,:) = [0.37 0.52 .53 .45];
cbpos = [0.93 0.52 0.015 .45];
axlabels = {'(a)','(b)','(c)','(d)'};


%%% Set up quivers
uveldens_zint_plot = uveldens_zint - uvel_zint*dens_ref;
vveldens_zint_plot = vveldens_zint - vvel_zint*dens_ref;
ETA = defineMOCgrid (XC,YC,SHELFICEtopo,bathy,false,true);
msk_uveldens = (ETA > 4.5) | ((ETA >3.3) & bathy<-700) | (ETA < 0.3);
uveldens_zint_plot(msk_uveldens) = NaN;
vveldens_zint_plot(msk_uveldens) = NaN;


qfac = 12;

Nxq = floor(Nx/qfac);
Nyq = floor(Ny/qfac);
u_plot = zeros(Nxq,Nyq);
v_plot = zeros(Nxq,Nyq);
x_plot = zeros(Nxq,Nyq);
y_plot = zeros(Nxq,Nyq);
lon_plot = zeros(Nxq,Nyq);
lat_plot = zeros(Nxq,Nyq);
ref_lat = -74;
ref_jidx = find(YC(1,:)>ref_lat,1,'first');
XCq = repmat(cumsum(DXG(:,ref_jidx),1),[1 Ny]);
YCq = cumsum(DYG,2);
for i=1:Nxq
  for j=1:Nyq
    tmp = uveldens_zint_plot(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    u_plot(i,j) = nanmean(tmp(:));
    tmp = vveldens_zint_plot(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    v_plot(i,j) = nanmean(tmp(:));
    tmp = XCq(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    x_plot(i,j) = nanmean(tmp(:));
    tmp = YCq(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    y_plot(i,j) = nanmean(tmp(:));   
    tmp = XC(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    lon_plot(i,j) = nanmean(tmp(:));
    tmp = YC(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    lat_plot(i,j) = nanmean(tmp(:));
  end
end

xMin = XCq(find(XC(:,1)>lonMin_b,1,'first'),1);
yMin = YCq(1,find(YC(1,:)>latMin_b,1,'first'));
xMax = XCq(find(XC(:,1)>lonMax_b,1,'first'),1);
yMax = YCq(1,find(YC(1,:)>latMax_b,1,'first'));

%%% Set up map plot
subplot('Position',axpos(1,:));
set(gca,'Color',[.8 .8 .8]);
bathy_plot = -bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
SHELFICE_plot = SHELFICEtopo;
SHELFICE_plot(SHELFICEtopo-bathy<=0) = NaN;
pcolor(XC,YC,pot_dens_bot);
shading interp;
hold on;
[C,h]=contour(XC,YC,bathy_plot,[100 400 700 1000 3000],'EdgeColor',[.5 .5 .5]);
clabel(C,h,'Color',[.5 .5 .5]);
[C,h]=contour(XC,YC,SHELFICE_plot,[-10 -10],'EdgeColor',[0    0.4470    0.7410],'LineStyle','-','LineWidth',2);
plot(secLons,secLats,'r:','LineWidth',2.5);
axis([lonMin_b lonMax_b latMin_b latMax_b]);
set(gca,'FontSize',fontsize);
ylabel('Latitude','interpreter','latex');
% xlabel('Longitude','interpreter','latex');

%%% Add colorbar and title
h = colorbar;
caxis([27.5 28.2]);
cmap = flip(haxby(50));
% cmap = cmocean('dense',50);
% cmap = cmap(1:14,:);
colormap(gca,cmap);
set(gca,'FontSize',10);
set(h,'Position',cbpos(1,:));
title(h,'kg/m$^3$','Fontsize',fontsize,'interpreter','latex');

%%% Add reference quiver
rq_lat = -78.5;
rq_lon = -25;
rq_amp = 10;
rq_iidx = find(lon_plot(:,1)>rq_lon,1,'first');
rq_jidx = find(lat_plot(1,:)>rq_lat,1,'first');
v_plot(rq_iidx,rq_jidx) = rq_amp;
text(rq_lon+1,rq_lat+0.2,['10 kg/m/s'],'FontSize',fontsize,'interpreter','latex');

%%% Add quivers
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
quiver(x_plot,y_plot,u_plot*cosd(ref_lat)./cosd(lat_plot),v_plot,3,'k-','LineWidth',1);
hold off
set(ax2,'Color','None');
set(ax2,'XLim',[xMin xMax]);
set(ax2,'YLim',[yMin yMax]);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
box off;
set(ax2,'FontSize',fontsize);
set(ax1,'FontSize',fontsize);
% set(ax1,'Color',[.8 .8 .8]);








%%%%%%%%%%%%%%%%%
%%% MAP PANEL %%%
%%%%%%%%%%%%%%%%%

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

%%% Plotting options
clim = [0 5000];
axpos(2,:) = [0 0.5 .3 .3];
icecolor = [186 242 239]/255;


latMin_d = -90;
latMax_d = -60;
lonMin_d = -180;
lonMax_d = 180;

%%% Set up map plot
subplot('Position',axpos(2,:));
axesm('stereo',...
  'fontsize',fontsize,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_d latMax_d], ...
  'MapLonLimit',[lonMin_d lonMax_d], ...   
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
pcolorm(double(LA),double(LO),-double(RTOPO_bathy_plot));
shading interp
colormap(gca,cmocean('deep',50));
caxis(clim);
set(gca,'Position',axpos(2,:));
% hold on;
ax1 = gca;

h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',[0.062000000000001,0.773218142548605,0.011000000000001,0.106198704103672]);
title(h,'m','Fontsize',fontsize,'interpreter','latex');

subplot('Position',axpos(2,:)+[0 0.5 0 0]);
ax2 = axesm('stereo',...
  'fontsize',fontsize,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_d latMax_d], ...
  'MapLonLimit',[lonMin_d lonMax_d], ...   
  'PLineLocation', 10, ...
  'MLineLocation', 60,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');  
axis off;

pcolorm(double(LA),double(LO),double(RTOPO_elev_plot));
hold on;
Nlon = 101;
dLon = (lonMax_b-lonMin_b)/(Nlon-1);
plotm([latMin_b*ones(1,Nlon) latMax_b*ones(1,Nlon) latMin_b],[lonMin_b:dLon:lonMax_b lonMax_b:-dLon:lonMin_b lonMin_b],'k--','LineWidth',2);
Nlon = 101;
dLon = (lonMax_o-lonMin_o)/(Nlon-1);
plotm([latMin_o*ones(1,Nlon) latMax_o*ones(1,Nlon) latMin_o],[lonMin_o:dLon:lonMax_o lonMax_o:-dLon:lonMin_o lonMin_o],'k-','LineWidth',3);
hold off;
shading interp;
% colormap(ax2,[220,243,255]/255);
colormap(ax2,icecolor);
% set(ax2,'Position',axpos-[0.001 0.003 0 -0.004]);
set(ax2,'Position',[-0.0030    0.4950    0.3060    0.3100]);













%%% Plotting options                
fontsize = 18;
labelspacing = 200;
axpos(3,:) = [0.08 0.285 0.82 0.18];
cb_pos = [0.93 0.285 0.015 0.18];
axlabel = '(c)';
clim = [-2 -1];
cmap = cmocean('thermal');

%%%%%%%%%%%%%%%%%%
%%% DENS SLICE %%%
%%%%%%%%%%%%%%%%%%

axes('Position',axpos(3,:));
pcolor(LO_sec,ZZ_sec,secT);
set(gca,'YDir','reverse');
shading interp;
hold on
[C,h] = contour(LO_sec,ZZ_sec,secS,[34.1:.1:35],'EdgeColor',[.8 .8 .8]);
clabel(C,h,'FontSize',fontsize-4,'LabelSpacing',300,'Color',[.8 .8 .8]);
plot(secLons,-secB,'k-','LineWidth',3);
Iidx = find(secI==0,1,'first');
plot(secLons(1:Iidx),-secI(1:Iidx),'k-','LineWidth',3);
hold off;
set(gca,'YLim',[0 630]);
set(gca,'XLim',[startLon endLon]);
caxis(clim);
colormap(gca,cmap);
cbhandle = colorbar;
set(cbhandle,'Position',cb_pos);
title(cbhandle,'$^\circ$C','Fontsize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);
ylabel('Depth (m)','interpreter','latex');
set(gca,'Color',[.8 .8 .8]);
% title(['Instantaneous potential temperature, salinity, ',datestr(datenum('01-Jan-2008')+iter*deltaT/86400)]);


fontsize = 18;
labelspacing = 200;
axpos(4,:) = [0.08 0.055 0.82 0.18];
cb_pos = [0.93 0.055 0.015 0.18];
axlabel = '(d)';
clim = [-0.2 0.2];
cmap = cmocean('balance',30);


%%%%%%%%%%%%%%%%%
%%% VEL SLICE %%%
%%%%%%%%%%%%%%%%%

axes('Position',axpos(4,:));
pcolor(LO_sec,ZZ_sec,secV);
set(gca,'YDir','reverse');
shading interp;
hold on
[C,h] = contour(LO_sec,ZZ_sec,secS,[34.1:.1:35],'EdgeColor','k');
clabel(C,h,'FontSize',fontsize-4,'LabelSpacing',300);
plot(secLons,-secB,'k-','LineWidth',3);
Iidx = find(secI==0,1,'first');
plot(secLons(1:Iidx),-secI(1:Iidx),'k-','LineWidth',3);
hold off;
set(gca,'YLim',[0 630]);
set(gca,'XLim',[startLon endLon]);
caxis(clim);
colormap(gca,cmap);
cbhandle = colorbar;
set(cbhandle,'Position',cb_pos);
title(cbhandle,'m/s','Fontsize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);
xlabel('Longitude','interpreter','latex')
ylabel('Depth (m)','interpreter','latex');
set(gca,'Color',[.8 .8 .8]);
% title(['Instantaneous potential temperature, salinity, ',datestr(datenum('01-Jan-2008')+iter*deltaT/86400)]);






%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.03 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)+0.01 axpos(2,2)-0.01 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.07 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.07 axpos(4,2)-0.05 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');