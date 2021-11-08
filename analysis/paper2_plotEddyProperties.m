%%%
%%% paper2_plotEddyProperties.m
%%%
%%% Plots properties of mesoscale eddies in our Weddell Sea simulations.
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_notides_RTOPO2';
loadexp;

%%% Load pre-computed EKE and EKE production products
outfname = [expname,'_EKE.mat'];
load(fullfile('products',outfname));

%%% Load pre-computed thickness fluxes
outfname = [expname,'_AABWcirc_.mat'];
load(fullfile('products',outfname));

%%% Load snapshots
% dumpIter = 1774980;
dumpIter = 1753067;
uvel = rdmdsWrapper(fullfile(exppath,'/results/U'),dumpIter);      
vvel = rdmdsWrapper(fullfile(exppath,'/results/V'),dumpIter);
theta = rdmdsWrapper(fullfile(exppath,'/results/T'),dumpIter);


  
  
  








%%%%%%%%%%%%%%%%%%%%
%%% CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%

%%% Compute relative vorticity
uvel(hFacW==0) = NaN;
vvel(hFacS==0) = NaN;
vort = zeros(Nx,Ny,Nr);
vort(:,2:Ny,:) = - (uvel(:,2:Ny,:)-uvel(:,1:Ny-1,:))./repmat(DYC(:,2:Ny),[1 1 Nr]);
vort(2:Nx,:,:) = vort(2:Nx,:,:) + (vvel(2:Nx,:,:)-vvel(1:Nx-1,:,:))./repmat(DXC(2:Nx,:),[1 1 Nr]);
hFacQ = zeros(Nx,Ny,Nr);
for i=2:Nx
  for j=2:Ny
    for k=1:Nr
      hFacQ(i,j,k) = min([hFacW(i,j,k),hFacW(i,j-1,k),hFacS(i,j,k),hFacS(i-1,j,k)]);
    end
  end
end
vort_avg = nansum(vort.*DRF.*hFacQ,3) ./ sum(DRF.*hFacQ,3);
Omega = 2*pi*366/365/86400;
ff = 2*Omega*sind(YG);
  
%%% Construct bottom temp
theta(hFacC==0) = NaN;
% thetabot = theta(:,:,40);
% for i=1:Nx
%   for j=1:Ny
%     if (isnan(thetabot(i,j)))
%       kmax = find(~isnan(theta(i,j,:)),1,'last');
%       if (~isempty(kmax))
%         thetabot(i,j) = theta(i,j,kmax);
%       end
%     end
%   end
% end

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

%%% Remove boundaries
hu_AABW(:,end-spongethickness+1:end) = NaN;
hv_AABW(:,end-spongethickness+1:end) = NaN;
hu_AABW(end-spongethickness+1:end,:) = NaN;
hv_AABW(end-spongethickness+1:end,:) = NaN;
u_AABW(:,end-spongethickness+1:end) = NaN;
v_AABW(:,end-spongethickness+1:end) = NaN;
u_AABW(end-spongethickness+1:end,:) = NaN;
v_AABW(end-spongethickness+1:end,:) = NaN;

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

%%% Cap flow speeds
uabs_AABW_eddy = sqrt(u_AABW_eddy.^2+v_AABW_eddy.^2);
uabs_AABW_eddy_max = 0.2;
u_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) = u_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) .* uabs_AABW_eddy_max./uabs_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max);
v_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) = v_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max) .* uabs_AABW_eddy_max./uabs_AABW_eddy(uabs_AABW_eddy>uabs_AABW_eddy_max);

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
    tmp = u_AABW_eddy(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    u_plot(i,j) = nanmean(tmp(:));
    tmp = v_AABW_eddy(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
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
% latMin = min(min(YC));
latMin = -78.5;
latMax = YC(1,end-spongethickness);
% lonMin = min(min(XC));
lonMin = -65;
lonMax = -10;
xMin = XCq(find(XC(:,1)>lonMin,1,'first'),1);
yMin = YCq(1,find(YC(1,:)>latMin,1,'first'));
xMax = XCq(find(XC(:,1)>lonMax,1,'first'),1);
yMax = YCq(1,find(YC(1,:)>latMax,1,'first'));
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(6,4);
axpos(1,:) = [0.06 0.71 .37 .26];
axpos(2,:) = [0.55 0.71 .37 .26];
axpos(3,:) = [0.06 0.38 .37 .26];
axpos(4,:) = [0.55 0.38 .37 .26];
axpos(5,:) = [0.06 0.05 .37 .26];
axpos(6,:) = [0.55 0.05 .37 .26];
cbpos = zeros(6,4);
cbpos(1,:) = [0.44 axpos(1,2) 0.015 axpos(1,4)];
cbpos(2,:) = [0.93 axpos(2,2) 0.015 axpos(2,4)];
cbpos(3,:) = [0.44 axpos(3,2) 0.015 axpos(3,4)];
cbpos(4,:) = [0.93 axpos(4,2) 0.015 axpos(4,4)];
cbpos(5,:) = [0.44 axpos(5,2) 0.015 axpos(5,4)];
cbpos(6,:) = [0.93 axpos(6,2) 0.015 axpos(6,4)];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};










%%% Set up the figure
figure(212)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[382          55        820         1050]);







%%% Set up pot temp plot
subplot('Position',axpos(1,:));
pcolor(XC,YC,theta(:,:,30));
shading interp;

%%% Add colorbar and title
h = colorbar;
colormap(gca,cmocean('thermal',50));
% colormap(gca,pmkmp(50,'Swtth'));
caxis([-2.5 1]);
set(h,'Position',cbpos(1,:));
% title(h,'$^\circ$C','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[C,h] = contour(XC,YC,SHELFICEtopo-bathy,[3000 4000],'EdgeColor',[.3 .3 .3]); 
clabel(C,h,'Color',[.3 .3 .3]);
[C,h] = contour(XC,YC,SHELFICEtopo-bathy,[500 1000 2000],'EdgeColor','w'); 
clabel(C,h,'Color','w');
hold off

%%% Labels
set(gca,'FontSize',fontsize);
% xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
axis([lonMin lonMax latMin latMax]);
title('Potential temperature at 340m (^oC)');
% set(gca,'Color',[.8 .8 .8]);








%%% Set up vorticity plot
subplot('Position',axpos(2,:));
pcolor(XG,YG,vort_avg./abs(ff));
shading interp;
caxis([-.2 .2]);

%%% Add colorbar and title
h = colorbar;
colormap(gca,cmocean('balance',50));% caxis([-1 1]*1e-3);
set(h,'Position',cbpos(2,:));

%%% Add bathymetry contours
hold on;
% [cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
[C,h] = contour(XC,YC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
hh = clabel(C,h,'Color',[.3 .3 .3]);
hold off

%%% Labels
set(gca,'FontSize',fontsize);
% xlabel('Longitude','interpreter','latex');
% ylabel('Latitude','interpreter','latex');
axis([lonMin lonMax latMin latMax]);
title(['Depth-averaged vortex Rossby number']);
% set(gca,'Color',[.8 .8 .8]);












%%% Set up eddy-induced transport velocity plot
subplot('Position',axpos(3,:));
bathy_plot = -bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
pcolor(XC,YC,SHELFICEtopo-bathy);
shading interp;
axis([lonMin lonMax latMin latMax]);
set(gca,'FontSize',fontsize);
ylabel('Latitude','interpreter','latex');
title(['Eddy-induced velocity below \sigma_\theta=',num2str(dens_AABW),' kg/m^3']);

%%% Add colorbar and title
h = colorbar;
caxis([0 5000]);
colormap(gca,flip(haxby(50)));
set(gca,'FontSize',10);
set(h,'Position',cbpos(3,:));
% title(h,'m','Fontsize',fontsize,'interpreter','latex');

%%% Add reference quiver
rq_lat = -77;
rq_lon = -21;
rq_amp = 0.2;
rq_iidx = find(lon_plot(:,1)>rq_lon,1,'first');
rq_jidx = find(lat_plot(1,:)>rq_lat,1,'first');
u_plot(rq_iidx,rq_jidx) = 0;
v_plot(rq_iidx,rq_jidx) = rq_amp;
text(rq_lon+2,rq_lat+1,[num2str(rq_amp),' m/s'],'FontSize',fontsize);

%%% Add quivers
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
quiver(x_plot,y_plot,u_plot*cosd(ref_lat)./cosd(lat_plot),v_plot,4,'k-');
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
subplot('Position',axpos(4,:));
EKE_zavg(sum(hFacW,3)==0) = NaN;
pcolor(XC,YC,log10(EKE_zavg));
shading interp;

%%% Add colorbar and title
h = colorbar;
colormap(gca,hot(40));
caxis([log10(1e-5) log10(1e-1)]);
set(h,'Position',cbpos(4,:));
% title(h,'m$^2$/s$^2$','Fontsize',fontsize,'interpreter','latex');
set(h,'YTick',[log10(1e-5) log10(1e-4) log10(1e-3) log10(1e-2) log10(1e-1)]);
set(h,'YTickLabel',{'10^-^5','10^-^4','10^-^3','10^-^2','10^-^1'});

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
title(['Depth-averaged eddy kinetic energy (m^2/s^2)']);
% set(gca,'Color',[.8 .8 .8]);





%%% Set up plot
subplot('Position',axpos(5,:));
PEtoEKE_zint = nansum(PEtoEKE.*hFacC.*DRF,3);
PEtoEKE_zint(sum(hFacW,3)==0) = NaN;
pcolor(XC,YC,PEtoEKE_zint);
shading interp;

%%% Add colorbar and title
h = colorbar;
colormap(gca,cmocean('balance',50));
caxis([-1 1]*1e-4);
set(h,'Position',cbpos(5,:));
% title(h,'mm/s','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[C,h] = contour(XC,YC,SHELFICEtopo-bathy,bathycntrs(1:end-1),'EdgeColor',[.3 .3 .3]); 
clabel(C,h,'Color',[.3 .3 .3]);
hold off

%%% Labels
set(gca,'FontSize',fontsize);
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
axis([lonMin lonMax latMin latMax]);
title(['Depth-integrated PE{\rightarrow}EKE (m^3/s^3)']);
% set(gca,'Color',[.8 .8 .8]);








%%% Set up plot
subplot('Position',axpos(6,:));
MKEtoEKE_zint = nansum(MKEtoEKE.*hFacC.*DRF,3);
MKEtoEKE_zint(sum(hFacW,3)==0) = NaN;
pcolor(XC,YC,MKEtoEKE_zint);
shading interp;
caxis([-1 1]*1e-4);

%%% Add colorbar and title
h = colorbar;
colormap(gca,cmocean('balance',50));% caxis([-1 1]*1e-3);
set(h,'Position',cbpos(6,:));
% title(h,'$\times$10$^{-7}$ m$^2$/s$^3$','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
% [cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
[C,h] = contour(XC,YC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
clabel(C,h,'Color',[.3 .3 .3]);
hold off

%%% Labels
set(gca,'FontSize',fontsize);
xlabel('Longitude','interpreter','latex');
% ylabel('Latitude','interpreter','latex');
axis([lonMin lonMax latMin latMax]);
title(['Depth-integrated MKE{\rightarrow}EKE (m^3/s^3)']);
% set(gca,'Color',[.8 .8 .8]);







%%% Inset plot for panel (e)
axes('Position',[0.23 0.2 0.18 0.095]);
pcolor(XC,YC,PEtoEKE_zint);
shading interp;
box on;

%%% Add colorbar and title
h = colorbar;
colormap(gca,cmocean('balance',50));
caxis([-.6 .6]*1e-3);

%%% Add bathymetry contours
hold on;
[C,h] = contour(XC,YC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.3 .3 .3]); 
clabel(C,h,'Color',[.3 .3 .3]);
plot([-40 -32 -32 -40 -40],[-74.8 -74.8 -73.3 -73.3 -74.8],'k-','LineWidth',1);
hold off

%%% Labels
set(gca,'FontSize',fontsize-4);
axis([-40 -32 -74.8 -73.3]);
% set(gca,'Color',[.8 .8 .8]);






%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.04 axpos(4,2)-0.04 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(5,1)-0.04 axpos(5,2)-0.04 0.03 0.03],'String',axlabels{5},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(6,1)-0.04 axpos(6,2)-0.04 0.03 0.03],'String',axlabels{6},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


  
  
  
  
  
  
  
  
