%%%
%%% paper2_plotDSWcirculation.m
%%%
%%% Plots horizontal circulation and water mass transformation in AABW layer.
%%%

%%%%%%%%%%%%%%%%%
%%% LOAD DATA %%%
%%%%%%%%%%%%%%%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_RTOPO2';
loadexp;
tmin = 1.05*t1year;
tmax = 9.05*t1year;

%%% Load pre-computed thickness fluxes
outfname = [expname,'_AABWcirc_.mat'];
load(fullfile('products',outfname));

%%% Load pre-computed time-mean T and S
outfname = [expname,'_TSfluxes'];
outfname = [outfname,'.mat'];
load(fullfile('products',outfname),'theta_tavg');
load(fullfile('products',outfname),'salt_tavg');
load(fullfile('products',outfname),'tflux_tavg');
load(fullfile('products',outfname),'sflux_tavg');
load(fullfile('products',outfname),'SHIfwFlx_tavg');
load(fullfile('products',outfname),'SHIhtFlx_tavg');

%%% Composite surface heat flux
htFlxDn = tflux_tavg;
htFlxDn(hFacC(:,:,1)==0) = -SHIhtFlx_tavg(hFacC(:,:,1)==0);
htFlxDn(sum(hFacC,3)==0) = NaN;

%%% Mean "sea surface" salinity for salt flux calculation
sst = NaN*ones(Nx,Ny);
sss = NaN*ones(Nx,Ny);
ssp = NaN*ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    kmin = find(hFacC(i,j,:)>0,1,'first');
    if (~isempty(kmin))
      sss(i,j) = salt_tavg(i,j,kmin);
      sst(i,j) = theta_tavg(i,j,kmin);
      ssp = -gravity*rhoConst*RC(kmin)/1e4;
    end
  end
end

%%% Composite virtual surface salt flux
sltFlxDn = sflux_tavg;
sltFlxDn(hFacC(:,:,1)==0) = SHIfwFlx_tavg(hFacC(:,:,1)==0).*sss(hFacC(:,:,1)==0);
sltFlxDn(sum(hFacC,3)==0) = NaN;

%%% Calculate buoyancy flux
[alpha,beta] = calcAlphaBeta(sss,sst,ssp);
Cp = 4e3;
bfluxT = gravity*(alpha.*htFlxDn/Cp/rhoConst);
bfluxS = - gravity*(beta.*sltFlxDn/rhoConst);
bflux = bfluxT + bfluxS;





%%%%%%%%%%%%%%%%%%%%
%%% CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%


%%% Density bounds for water masses
dens_AABW = 27.85;

%%% Density grid indices for water masses
k_AABW = find(dens_levs==dens_AABW);

%%% AABW layer thickness
H_AABW_w = sum(uthic_tavg(:,:,k_AABW:end),3);
H_AABW_s = sum(vthic_tavg(:,:,k_AABW:end),3);

%%% Transports and TWA velocities 
hu_AABW = sum(uflux_tavg(:,:,k_AABW:end),3);
u_AABW = hu_AABW ./ H_AABW_w;
hv_AABW = sum(vflux_tavg(:,:,k_AABW:end),3);
v_AABW = hv_AABW ./ H_AABW_s;

%%% Diapycnal velocity
w_AABW = ( hu_AABW(1:Nx,1:Ny).*DYG(1:Nx,1:Ny) ...
         - hu_AABW([2:Nx 1],1:Ny).*DYG([2:Nx 1],1:Ny) ...
         + hv_AABW(1:Nx,1:Ny).*DXG(1:Nx,1:Ny) ...
         - hv_AABW(1:Nx,[2:Ny 1]).*DYG(1:Nx,[2:Ny 1]) ) ./ RAC;

%%% Remove boundaries
hu_AABW(:,end-spongethickness+1:end) = NaN;
hv_AABW(:,end-spongethickness+1:end) = NaN;
hu_AABW(end-spongethickness+1:end,:) = NaN;
hv_AABW(end-spongethickness+1:end,:) = NaN;
u_AABW(:,end-spongethickness+1:end) = NaN;
v_AABW(:,end-spongethickness+1:end) = NaN;
u_AABW(end-spongethickness+1:end,:) = NaN;
v_AABW(end-spongethickness+1:end,:) = NaN;

qfac = 12*round(Nx/600);

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
    tmp = hu_AABW(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    u_plot(i,j) = nanmean(tmp(:));
    tmp = hv_AABW(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
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
    
  










%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
latMin = min(min(YC));
latMax = YC(1,end-spongethickness);
lonMin = min(min(XC));
lonMax = -10;
xMin = 0;
yMin = 0;
xMax = XCq(find(XC(:,1)>lonMax,1,'first'),1);
yMax = YCq(1,find(YC(1,:)>latMax,1,'first'));
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.05 0.54 .38 .42];
axpos(2,:) = [0.55 0.54 .38 .42];
axpos(3,:) = [0.05 0.05 .38 .42];
axpos(4,:) = [0.55 0.05 .38 .42];
cbpos = zeros(4,4);
cbpos(1,:) = [0.45 0.54 0.015 .42];
cbpos(2,:) = [0.95 0.54 0.015 .42];
cbpos(3,:) = [0.45 0.05 0.015 .42];
cbpos(4,:) = [0.95 0.05 0.015 .42];
axlabels = {'(a)','(b)','(c)','(d)'};








%%% Set up the figure
figure(206)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[382          55        1120         930]);





%%% Set up map plot
subplot('Position',axpos(1,:));
bathy_plot = -bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
pcolor(XC,YC,SHELFICEtopo-bathy);
shading interp;
axis([lonMin lonMax latMin latMax]);
set(gca,'FontSize',fontsize);
ylabel('Latitude','interpreter','latex');
title(['Transport below \sigma_\theta=',num2str(dens_AABW),' kg/m^3']);

%%% Add colorbar and title
h = colorbar;
caxis([0 5000]);
colormap(gca,flip(haxby(50)));
set(gca,'FontSize',10);
set(h,'Position',cbpos(1,:));
title(h,'m','Fontsize',fontsize,'interpreter','latex');

%%% Add reference quiver
rq_lat = -83;
rq_lon = -25;
rq_amp = 50;
rq_iidx = find(lon_plot(:,1)>rq_lon,1,'first');
rq_jidx = find(lat_plot(1,:)>rq_lat,1,'first');
v_plot(rq_iidx,rq_jidx) = rq_amp;
text(rq_lon+2,rq_lat+1,[num2str(rq_amp),' m^2/s'],'FontSize',fontsize);

%%% Add quivers
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
quiver(x_plot,y_plot,u_plot*cosd(ref_lat)./cosd(lat_plot),v_plot,2,'k-');
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





%%% Set up plot
subplot('Position',axpos(2,:));
H_AABW_w(sum(hFacW,3)==0) = NaN;
pcolor(XC,YC,log10(H_AABW_w));
shading interp;

%%% Add colorbar and title
h = colorbar;
colormap(gca,cmocean('amp',10));
caxis([1 log10(3000)]);
set(h,'Position',cbpos(2,:));
title(h,'m','Fontsize',fontsize,'interpreter','latex');
set(h,'YTick',[1 log10(30) 2 log10(300) 3 log10(3000)]);
set(h,'YTickLabel',{'10','30','100','300','1000','3000'});

%%% Add bathymetry contours
hold on;
% [cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
[C,h] = contour(XC,YC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
% hh = clabelm(cs,C);
hh = clabel(C,h,'Color',[0.3 0.3 0.3]);
% set(hh,'fontsize',8,'Color',[0.3 0.3 0.3],'BackgroundColor','none','Edgecolor','none')               
hold off

%%% Labels
set(gca,'FontSize',fontsize);
% xlabel('Longitude','interpreter','latex');
% ylabel('Latitude','interpreter','latex');
axis([lonMin lonMax latMin latMax]);
title(['Layer thickness below \sigma_\theta=',num2str(dens_AABW),' kg/m^3']);
% set(gca,'Color',[.8 .8 .8]);







%%% Set up plot
subplot('Position',axpos(3,:));
w_AABW(sum(hFacC,3)==0) = NaN;
pcolor(XC,YC,w_AABW*1e3);
shading interp;

%%% Add colorbar and title
h = colorbar;
colormap(gca,cmocean('balance',200));
caxis([-.5 .5]);
set(h,'Position',cbpos(3,:));
title(h,'mm/s','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[C,h] = contour(XC,YC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
clabel(C,h);
hold off

%%% Labels
set(gca,'FontSize',fontsize);
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
axis([lonMin lonMax latMin latMax]);
title(['Diapycnal velocity on \sigma_\theta=',num2str(dens_AABW),' kg/m^3']);
% set(gca,'Color',[.8 .8 .8]);








%%% Set up plot
subplot('Position',axpos(4,:));
pcolor(XC,YC,bflux*1e7);
shading interp;
caxis([-1 1]);

%%% Add colorbar and title
h = colorbar;
colormap(gca,cmocean('delta',50));% caxis([-1 1]*1e-3);
set(h,'Position',cbpos(4,:));
title(h,'$\times$10$^{-7}$ m$^2$/s$^3$','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
% [cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
[C,h] = contour(XC,YC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
hh = clabel(C,h);
hold off

%%% Labels
set(gca,'FontSize',fontsize);
xlabel('Longitude','interpreter','latex');
% ylabel('Latitude','interpreter','latex');
axis([lonMin lonMax latMin latMax]);
title(['Ocean surface buoyancy flux']);
% set(gca,'Color',[.8 .8 .8]);








%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.04 axpos(4,2)-0.04 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

