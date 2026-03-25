%%%
%%% paper3_plotDepthIntegratedHeatFluxMaps.m
%%%
%%% Plots maps of depth-integrate heat fluxes, divided into positive and
%%% negative components.
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;

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

%%% Load pre-computed fluxes
load(fullfile('products',outfname),'uvelth_tavg','vvelth_tavg','uvel_tavg','vvel_tavg','theta_tavg','eta','ETA');

%%% Reference surface freezing temperature
theta0 = -1.9;

%%% Compute depth-integrated heat fluxes
uvelth_zint = sum((uvelth_tavg-uvel_tavg.*theta0).*DRF.*hFacW,3);
vvelth_zint = sum((vvelth_tavg-vvel_tavg.*theta0).*DRF.*hFacS,3);

%%% Compute depth-integrated positive/negative zonal heat fluxes
theta_u = 0.5*(theta_tavg(1:Nx,:,:)+theta_tavg([Nx 1:Nx-1],:,:));
uvelth_pos_zint = (uvelth_tavg-uvel_tavg.*theta0);
uvelth_pos_zint(theta_u<theta0) = 0;
uvelth_pos_zint = sum(uvelth_pos_zint.*DRF.*hFacW,3);
uvelth_neg_zint = (uvelth_tavg-uvel_tavg.*theta0);
uvelth_neg_zint(theta_u>=theta0) = 0;
uvelth_neg_zint = sum(uvelth_neg_zint.*DRF.*hFacW,3);
uvelth_eddy_zint = (uvelth_tavg-uvel_tavg.*theta_u);
uvelth_eddy_zint = sum(uvelth_eddy_zint.*DRF.*hFacW,3);
clear('theta_u');

%%% Compute depth-integrated positive/negative meridional heat fluxes
theta_v = 0.5*(theta_tavg(:,1:Ny,:)+theta_tavg(:,[Ny 1:Ny-1],:));
vvelth_pos_zint = (vvelth_tavg-vvel_tavg.*theta0);
vvelth_pos_zint(theta_v<theta0) = 0;
vvelth_pos_zint = sum(vvelth_pos_zint.*DRF.*hFacS,3);
vvelth_neg_zint = (vvelth_tavg-vvel_tavg.*theta0);
vvelth_neg_zint(theta_v>=theta0) = 0;
vvelth_neg_zint = sum(vvelth_neg_zint.*DRF.*hFacS,3);
vvelth_eddy_zint = (vvelth_tavg-vvel_tavg.*theta_v);
vvelth_eddy_zint = sum(vvelth_eddy_zint.*DRF.*hFacS,3);
clear('theta_v');

%%% Free up memory
clear('uvelth_tavg','vvelth_tavg','uvel_tavg','vvel_tavg','theta');






%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
% latMin = min(min(YC));
latMin = -78.5;
% latMax = YC(1,end-spongethickness);
latMax = -72;
% lonMin = min(min(XC));
lonMin = -65;
lonMax = -25;
fontsize = 18;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.26 0.56 .4 .41];
axpos(2,:) = [0.06 0.07 .4 .41];
axpos(3,:) = [0.56 0.07 .4 .41];
cbpos = [0.7 0.56 0.015 .41];
axlabels = {'\textbf{A}','\textbf{B}','\textbf{C}'};


%%% Heat capacity
Cp = 4000;




%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(220)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[986         474        1126         771]);


%%%%%%% PANEL A %%%%%%%

%%% Set up quivers
uvelth_zint_plot = rhoConst*Cp*uvelth_zint/1e6;
vvelth_zint_plot = rhoConst*Cp*vvelth_zint/1e6;
uvelth_zint_plot(ETA>3.3) = NaN;
vvelth_zint_plot(ETA>3.3) = NaN;
uvelth_zint_plot(ETA<0) = NaN;
vvelth_zint_plot(ETA<0) = NaN;

qfac = 25;

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
    tmp = uvelth_zint_plot(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    u_plot(i,j) = nanmean(tmp(:));
    tmp = vvelth_zint_plot(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
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

xMin = XCq(find(XC(:,1)>lonMin,1,'first'),1);
yMin = YCq(1,find(YC(1,:)>latMin,1,'first'));
xMax = XCq(find(XC(:,1)>lonMax,1,'first'),1);
yMax = YCq(1,find(YC(1,:)>latMax,1,'first'));

%%% Set up map plot
subplot('Position',axpos(1,:));
bathy_plot = -bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
pcolor(XC,YC,SHELFICEtopo-bathy);
shading interp;
axis([lonMin lonMax latMin latMax]);
set(gca,'FontSize',fontsize);
ylabel('Latitude','interpreter','latex');
title(['Total heat transport']);

%%% Add colorbar and title
h = colorbar;
caxis([0 5000]);
colormap(gca,flip(haxby(50)));
set(gca,'FontSize',10);
set(h,'Position',cbpos(1,:));
title(h,'m','Fontsize',fontsize,'interpreter','latex');

%%% Add reference quiver
rq_lat = -78.2;
rq_lon = -33;
rq_amp = 20;
rq_iidx = find(lon_plot(:,1)>rq_lon,1,'first');
rq_jidx = find(lat_plot(1,:)>rq_lat,1,'first');
v_plot(rq_iidx,rq_jidx) = rq_amp;
text(rq_lon+1,rq_lat+0.2,[num2str(rq_amp),' MW/m'],'FontSize',fontsize);

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

text(ax1,-34,-75,'FT','FontSize',fontsize);
text(ax1,-48,-77,'BB','FontSize',fontsize);
text(ax1,-60,-74.7,'RD','FontSize',fontsize);






%%%%%%% PANEL B %%%%%%%

%%% Set up quivers
uvelth_zint_plot = rhoConst*Cp*uvelth_pos_zint/1e6;
vvelth_zint_plot = rhoConst*Cp*vvelth_pos_zint/1e6;
uvelth_zint_plot(ETA>3.3) = NaN;
vvelth_zint_plot(ETA>3.3) = NaN;
uvelth_zint_plot(ETA<0) = NaN;
vvelth_zint_plot(ETA<0) = NaN;

qfac = 25;

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
    tmp = uvelth_zint_plot(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    u_plot(i,j) = nanmean(tmp(:));
    tmp = vvelth_zint_plot(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
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

xMin = XCq(find(XC(:,1)>lonMin,1,'first'),1);
yMin = YCq(1,find(YC(1,:)>latMin,1,'first'));
xMax = XCq(find(XC(:,1)>lonMax,1,'first'),1);
yMax = YCq(1,find(YC(1,:)>latMax,1,'first'));

%%% Set up map plot
subplot('Position',axpos(2,:));
bathy_plot = -bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
pcolor(XC,YC,SHELFICEtopo-bathy);
shading interp;
axis([lonMin lonMax latMin latMax]);
set(gca,'FontSize',fontsize);
ylabel('Latitude','interpreter','latex');
xlabel('Longitude','interpreter','latex');
title(['"Warm" heat transport']);

%%% Add reference quiver
rq_lat = -78.2;
rq_lon = -33;
rq_amp = 20;
rq_iidx = find(lon_plot(:,1)>rq_lon,1,'first');
rq_jidx = find(lat_plot(1,:)>rq_lat,1,'first');
v_plot(rq_iidx,rq_jidx) = rq_amp;
text(rq_lon+1,rq_lat+0.2,[num2str(rq_amp),' MW/m'],'FontSize',fontsize);

%%% Add colorbar and title
caxis([0 5000]);
colormap(gca,flip(haxby(50)));
set(gca,'FontSize',10);

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

text(ax1,-34,-75,'FT','FontSize',fontsize);
text(ax1,-48,-77,'BB','FontSize',fontsize);
text(ax1,-60,-74.7,'RD','FontSize',fontsize);






%%%%%%% PANEL C %%%%%%%

%%% Set up quivers
uvelth_zint_plot = rhoConst*Cp*uvelth_neg_zint/1e6;
vvelth_zint_plot = rhoConst*Cp*vvelth_neg_zint/1e6;
uvelth_zint_plot(ETA>3.3) = NaN;
vvelth_zint_plot(ETA>3.3) = NaN;
uvelth_zint_plot(ETA<0) = NaN;
vvelth_zint_plot(ETA<0) = NaN;

qfac = 25;

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
    tmp = uvelth_zint_plot(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
    u_plot(i,j) = nanmean(tmp(:));
    tmp = vvelth_zint_plot(round(qfac*(i-1))+1:round(qfac*i),round(qfac*(j-1))+1:round(qfac*j));
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

v_plot(10,10) = 50;

xMin = XCq(find(XC(:,1)>lonMin,1,'first'),1);
yMin = YCq(1,find(YC(1,:)>latMin,1,'first'));
xMax = XCq(find(XC(:,1)>lonMax,1,'first'),1);
yMax = YCq(1,find(YC(1,:)>latMax,1,'first'));

%%% Set up map plot
subplot('Position',axpos(3,:));
bathy_plot = -bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
pcolor(XC,YC,SHELFICEtopo-bathy);
shading interp;
axis([lonMin lonMax latMin latMax]);
set(gca,'FontSize',fontsize);
xlabel('Longitude','interpreter','latex');
title(['"Cold" heat transport']);

%%% Add colorbar and title
caxis([0 5000]);
colormap(gca,flip(haxby(50)));
set(gca,'FontSize',10);

%%% Add reference quiver
rq_lat = -78.2;
rq_lon = -33;
rq_amp = 20;
rq_iidx = find(lon_plot(:,1)>rq_lon,1,'first');
rq_jidx = find(lat_plot(1,:)>rq_lat,1,'first');
v_plot(rq_iidx,rq_jidx) = rq_amp;
text(rq_lon+1,rq_lat+0.2,[num2str(rq_amp),' MW/m'],'FontSize',fontsize);

%%% Add quivers
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
quiver(x_plot,y_plot,u_plot*cosd(ref_lat)./cosd(lat_plot),v_plot,1.5,'k-','LineWidth',1);
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

text(ax1,-34,-75,'FT','FontSize',fontsize);
text(ax1,-48,-77,'BB','FontSize',fontsize);
text(ax1,-60,-74.7,'RD','FontSize',fontsize);






%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.05 axpos(1,2)-0.03 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.05 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.05 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

