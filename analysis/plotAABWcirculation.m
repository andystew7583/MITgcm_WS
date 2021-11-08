%%%
%%% plotAABWcirculation.m
%%%
%%% Plots horizontal circulation in isopycnal layers.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_notides_RTOPO2';
loadexp;

%%% Load pre-computed data
outfname = [expname,'_AABWcirc_.mat'];
load(fullfile('products',outfname));

%%% Density bounds for water masses
dens_CDW_min = 27.7;
dens_CDW_max = 27.85;
dens_AABW = 27.85;

%%% Density grid indices for water masses
k_CDW_min = find(dens_levs==dens_CDW_min);
k_CDW_max = find(dens_levs==dens_CDW_max) - 1;
k_AASW_max = k_CDW_min - 1;
k_AABW = find(dens_levs==dens_AABW);

%%% AABW layer thickness
H_AABW_w = sum(uthic_tavg(:,:,k_AABW:end),3);
H_AABW_s = sum(vthic_tavg(:,:,k_AABW:end),3);
H_CDW_w = sum(uthic_tavg(:,:,k_CDW_min:k_CDW_max),3);
H_CDW_s = sum(vthic_tavg(:,:,k_CDW_min:k_CDW_max),3);

%%% Transports and TWA velocities in different water mass layers
hu_AASW = sum(uflux_tavg(:,:,1:k_AASW_max),3);
u_AASW = hu_AASW ./ sum(uthic_tavg(:,:,1:k_AASW_max),3);
hv_AASW = sum(vflux_tavg(:,:,1:k_AASW_max),3);
v_AASW = hv_AASW ./ sum(vthic_tavg(:,:,1:k_AASW_max),3);
hu_CDW = sum(uflux_tavg(:,:,k_CDW_min:k_CDW_max),3);
u_CDW = hu_CDW ./ sum(uthic_tavg(:,:,k_CDW_min:k_CDW_max),3);
hv_CDW = sum(vflux_tavg(:,:,k_CDW_min:k_CDW_max),3);
v_CDW = hv_CDW ./ sum(vthic_tavg(:,:,k_CDW_min:k_CDW_max),3);
hu_AABW = sum(uflux_tavg(:,:,k_AABW:end),3);
u_AABW = hu_AABW ./ H_AABW_w;
hv_AABW = sum(vflux_tavg(:,:,k_AABW:end),3);
v_AABW = hv_AABW ./ H_AABW_s;

%%% Mean/eddy decomposition
hu_AABW_mean = sum(uflux_mean(:,:,k_AABW:end),3);
u_AABW_mean = hu_AABW_mean ./ H_AABW_w;
hv_AABW_mean = sum(vflux_mean(:,:,k_AABW:end),3);
v_AABW_mean = hv_AABW_mean ./ H_AABW_s;
hu_AABW_eddy = hu_AABW - hu_AABW_mean;
hv_AABW_eddy = hv_AABW - hv_AABW_mean;
u_AABW_eddy = hu_AABW_eddy ./ H_AABW_w;
v_AABW_eddy = hv_AABW_eddy ./ H_AABW_s;
u_AABW_eddy(H_AABW_w==0) = NaN;
v_AABW_eddy(H_AABW_s==0) = NaN;

hu_CDW_mean = sum(uflux_mean(:,:,k_CDW_min:k_CDW_max),3);
u_CDW_mean = hu_CDW_mean ./ H_CDW_w;
hv_CDW_mean = sum(vflux_mean(:,:,k_CDW_min:k_CDW_max),3);
v_CDW_mean = hv_CDW_mean ./ H_CDW_s;
hu_CDW_eddy = hu_CDW - hu_CDW_mean;
hv_CDW_eddy = hv_CDW - hv_CDW_mean;
u_CDW_eddy = hu_CDW_eddy ./ H_CDW_w;
v_CDW_eddy = hv_CDW_eddy ./ H_CDW_s;
u_CDW_eddy(H_CDW_w==0) = NaN;
v_CDW_eddy(H_CDW_s==0) = NaN;

%%% Cap flow speeds
uabs_AABW_eddy = sqrt(u_AABW_eddy.^2+v_AABW_eddy.^2);
uabs_AABW_eddy_max = 0.3;
u_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) = u_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) .* uabs_AABW_eddy_max./uabs_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max);
v_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) = v_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) .* uabs_AABW_eddy_max./uabs_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max);

uabs_CDW_eddy = sqrt(u_CDW_eddy.^2+v_CDW_eddy.^2);
uabs_CDW_eddy_max = 0.3;
u_CDW_eddy(uabs_CDW_eddy>uabs_CDW_eddy_max) = u_CDW_eddy(uabs_CDW_eddy>uabs_CDW_eddy_max) .* uabs_CDW_eddy_max./uabs_CDW_eddy(uabs_CDW_eddy>uabs_CDW_eddy_max);
v_CDW_eddy(uabs_CDW_eddy>uabs_CDW_eddy_max) = v_CDW_eddy(uabs_CDW_eddy>uabs_CDW_eddy_max) .* uabs_CDW_eddy_max./uabs_CDW_eddy(uabs_CDW_eddy>uabs_CDW_eddy_max);


%%% Attempt at an AABW streamfunction
psi_AABW = zeros(Nx,Ny);
psi_AABW(:,2:Ny) = -cumsum(hu_AABW(:,1:Ny-1).*DYG(:,1:Ny-1),2);
psi_AABW(2:Nx,:) = cumsum(hv_AABW(1:Nx-1,:).*DYG(1:Nx-1,:),1);

%%% Diapycnal velocity
w_AABW = ( hu_AABW(1:Nx,1:Ny).*DYG(1:Nx,1:Ny) ...
         - hu_AABW([2:Nx 1],1:Ny).*DYG([2:Nx 1],1:Ny) ...
         + hv_AABW(1:Nx,1:Ny).*DXG(1:Nx,1:Ny) ...
         - hv_AABW(1:Nx,[2:Ny 1]).*DYG(1:Nx,[2:Ny 1]) ) ./ RAC;

%%% Subsampling options for quiver plots
arrowspacing = 3*round(Nx/304); %%% Scale with grid sizes
xidx = 1:arrowspacing:Nx-arrowspacing;
yidx = 1:arrowspacing:Ny-arrowspacing;

figure(71);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_AABW(xidx,yidx),v_AABW(xidx,yidx));

figure(72);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_AABW(xidx,yidx),hv_AABW(xidx,yidx));

figure(73);
pcolor(XC,YC,w_AABW);
shading interp;
colorbar;
colormap redblue;
caxis([-2 2]*1e-4);

figure(74);
pcolor(XG,YG,psi_AABW/1e6);
shading interp;
colorbar;
colormap redblue;
caxis([-10 10]);


figure(75);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_CDW(xidx,yidx),v_CDW(xidx,yidx));

figure(76);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_CDW(xidx,yidx),hv_CDW(xidx,yidx));


figure(77);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_AASW(xidx,yidx),v_AASW(xidx,yidx));

figure(78);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_AASW(xidx,yidx),hv_AASW(xidx,yidx));


figure(79);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_AABW_eddy(xidx,yidx),v_AABW_eddy(xidx,yidx));
hold on;
[C,h] = contour(XC,YC,bathy,[-4000 -3000 -2000 -1000 -500],'EdgeColor','k');
clabel(C,h);
hold off;

figure(80);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_AABW_eddy(xidx,yidx),hv_AABW_eddy(xidx,yidx));
hold on;
[C,h] = contour(XC,YC,bathy,[-4000 -3000 -2000 -1000 -500],'EdgeColor','k');
clabel(C,h);
hold off;

figure(81);
quiver(XC(xidx,yidx),YC(xidx,yidx),u_AABW_mean(xidx,yidx),v_AABW_mean(xidx,yidx));

figure(82);
quiver(XC(xidx,yidx),YC(xidx,yidx),hu_AABW_mean(xidx,yidx),hv_AABW_mean(xidx,yidx));





















%%% Subsample vectors for quiver plot
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
    tmp = hu_CDW_eddy(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    u_plot(i,j) = nanmean(tmp(:));
    tmp = hv_CDW_eddy(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
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

%%% Plotting options
latMin = -78.5;
latMax = YC(1,end-spongethickness);
lonMin = -65;
lonMax = -10;
xMin = XCq(find(XC(:,1)>lonMin,1,'first'),1);
yMin = YCq(1,find(YC(1,:)>latMin,1,'first'));
xMax = XCq(find(XC(:,1)>lonMax,1,'first'),1);
yMax = YCq(1,find(YC(1,:)>latMax,1,'first'));
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];



figure(83);
clf;
bathy_plot = -bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
pcolor(XC,YC,SHELFICEtopo-bathy);
shading interp;
axis([lonMin lonMax latMin latMax]);
set(gca,'FontSize',fontsize);
ylabel('Latitude','interpreter','latex');
xlabel('Longitude','interpreter','latex');
title(['Eddy-induced transport between \sigma_\theta=',num2str(dens_CDW_min),' kg/m^3 and \sigma_\theta=',num2str(dens_AABW),' kg/m^3']);

%%% Add colorbar and title
h = colorbar;
caxis([0 5000]);
colormap(gca,flip(haxby(50)));
set(gca,'FontSize',10);
% title(h,'m','Fontsize',fontsize,'interpreter','latex');

%%% Add reference quiver
rq_lat = -77;
rq_lon = -21;
rq_amp = 10;
rq_iidx = find(lon_plot(:,1)>rq_lon,1,'first');
rq_jidx = find(lat_plot(1,:)>rq_lat,1,'first');
u_plot(rq_iidx,rq_jidx) = 0;
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
% set(ax1,'Color',[.8 .8 .8]
