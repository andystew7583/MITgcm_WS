%%%
%%% paper2_plotWaterMasses.m
%%%
%%% Plots simulated water mass distribution and compares with observations
%%%

%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
% expname = 'hires_seq_onethird_RTOPO2';
% expname = 'hires_seq_onesixth_notides_RTOPO2';
expname = 'hires_seq_onetwelfth_RTOPO2';
loadexp;






%%%%%%%%%%%%%%%%%%%
%%%%% OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%

%%% T/S grids
dT = 0.025;
dS = 0.005;
Tmin = -3;
Tmax = 2;
Smin = 33.7;
Smax = 35;
SS = Smin:dS:Smax;
TT = Tmin:dT:Tmax;
NS = length(SS);
NT = length(TT);

%%% Quasi-latitude horizontal grid
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,false);
ETA = repmat(ETA,[1 1 Nr]);












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD REFERENCE STATE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load data filess
outfname = [expname,'_TSfluxes'];
outfname = [outfname,'.mat'];
pt_ref = load(fullfile('products',outfname),'theta_tavg');
pt_ref = pt_ref.theta_tavg;
ss_ref = load(fullfile('products',outfname),'salt_tavg');
ss_ref = ss_ref.salt_tavg;
outfname = [expname,'_SeasonalStrat'];
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

%%% Pressure is just Boussinesq hydrostatic reference pressure
pp_ref = -gravity*rhoConst*repmat(RC,[Nx Ny 1])/1e4;

%%% Remove topography
pt_ref(hFacC==0) = NaN;
ss_ref(hFacC==0) = NaN;
theta_djf(hFacC==0) = NaN;
theta_jja(hFacC==0) = NaN;
salt_djf(hFacC==0) = NaN;
salt_jja(hFacC==0) = NaN;

%%% Potential density
pd_ref = densjmd95(ss_ref,pt_ref,-gravity*rhoConst*repmat(RC(1),[Nx Ny Nr])/1e4);
pd_djf = densjmd95(salt_djf,theta_djf,-gravity*rhoConst*repmat(RC(1),[Nx Ny Nr])/1e4);
pd_jja = densjmd95(salt_jja,theta_jja,-gravity*rhoConst*repmat(RC(1),[Nx Ny Nr])/1e4);















%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% COMPUTATIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Volumes in T/S bins
TSvol_full = binByVolume(ss_ref,pt_ref,[], ...
                  Smin,Smax,dS,Tmin,Tmax,dT, ...
                  RAC,DRF,hFacC);               
                                             
%%% Load CTD data
load_NARE85_IC92;
                







                


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE THE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
fontsize = 12;
framepos = [417    34   864   926];
labelspacing = 200;
axpos = zeros(6,4);
axpos(1,:) = [0.08 0.66 0.39 0.3];
axpos(2,:) = [0.59 0.65 0.37 0.3];
axpos(3,:) = [0.08 0.35 0.38 0.25];
axpos(4,:) = [0.54 0.35 0.38 0.25];
axpos(5,:) = [0.08 0.05 0.38 0.25];
axpos(6,:) = [0.54 0.05 0.38 0.25];
cb1_pos = [0.48 0.66 0.02 0.3];
cb2_pos = [0.93 0.65 0.02 0.15];
cb3_pos = [0.94 0.35 0.02 0.25];
cb4_pos = [0.94 0.05 0.02 0.25];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

figure(203);                
clf;
set(gcf,'Position',framepos);            
                
                



%%%
%%% T/S DIAGRAM
%%%

[TTT,SSS] = meshgrid(TT,SS);
DDD = densjmd95(SSS,TTT,-RC(1)*gravity*rhoConst/1e4*ones(NS,NT)) - 1000;
TSvol_plot = TSvol_full;
TSvol_plot(TSvol_full==0) = NaN;

subplot('Position',axpos(1,:));
pcolor(SSS,TTT,log10(TSvol_plot));
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor',[0.5 0.5 0.5]);
hold off;
clabel(C,h,'Color',[.5 .5 .5]);
% pcolor(SSS,TTT,(TSvol));
shading flat;
xlabel('Salinity (g/kg)');
ylabel('Potential temperature (^oC)');
cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
% colormap(gca,flip(haxby));
colormap(gca,cmocean('amp'));
% colormap(gca,flip(hot));
caxis([10 14]);
set(cbhandle,'YTick',[10 11 12 13 14]);
set(cbhandle,'YTickLabel',{'10^1^0','10^1^1','10^1^2','10^1^3','10^1^4'});
title(cbhandle,'m^3')
axis([33.7 35 -2.7 1.2]);
set(gca,'Position',axpos(1,:));
set(gca,'FontSize',fontsize);
text(34.8,-1.8,'HSSW','FontSize',fontsize-2);
text(34.71,0.6,'WDW','FontSize',fontsize-2);
text(34.75,-2.3,'ISW','FontSize',fontsize-2);
text(34.4,-2,'WW','FontSize',fontsize-2);
text(34.1,-0.5,'AASW','FontSize',fontsize-2);




%%%
%%% MAP 
%%%


%%% Plotting options
latMin = -78;
latMax = -66;
lonMin = -60;
lonMax = -25;
clim = [0 5000];

%%% Set up map plot
subplot('Position',axpos(2,:));
axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)

%%% Add bathymetry contours
bathy_plot = -bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
bathy_plot((SHELFICEtopo<0) & (SHELFICEtopo>bathy)) = NaN;
pcolorm(YC,XC,bathy_plot);

%%% Add sections
hold on;
plotm(NARE85_lat,NARE85_lon,'ko-','MarkerSize',2,'MarkerFaceColor','k');
plotm(IC92_lat,IC92_lon,'ko-','MarkerSize',2,'MarkerFaceColor','k');
hold off
textm(-68.2,-55,'IC92','FontSize',fontsize);
textm(-73.8,-35,'NARE85','FontSize',fontsize);

%%% Add colorbar and title
cbhandle = colorbar;
set(cbhandle,'Position',cb2_pos);
caxis(clim);
colormap(gca,flip(haxby(50)));
set(gca,'FontSize',fontsize);
tightmap;
title(cbhandle,'m','Fontsize',14,'interpreter','latex');
set(gca,'Position',axpos(2,:));










%%%
%%% MODEL, FILCHNER OUTFLOW
%%%

pt_min = -2.1;
pt_max = -1;

%%% Meshgrid for plots
[ZZ,XX]=meshgrid(zz,XC(:,1));
jidx = find(YC(1,:)>-74.75,1,'first')
YC(1,jidx)
dmax = 650;

for i=1:Nx
  kbot = find(hFacC(i,jidx,:)>0,1,'last');
  if (~isempty(kbot))
    ZZ(i,kbot) = -sum(hFacC(i,jidx,1:kbot-1).*DRF(1:kbot-1),3)-hFacC(i,jidx,kbot)*DRF(kbot)/2;
  end
end

%%% Make plots
subplot('Position',axpos(3,:));
colormap jet;
pcolor(XX,-ZZ,squeeze(theta_djf(:,jidx,:)))
shading interp
hold on;
[C,h]=contour(XX,-ZZ,squeeze(pd_djf(:,jidx,:)-1000),[27.4:.1:28.4],'EdgeColor',[.8 .8 .8]);
clabel(C,h,'Color',[.8 .8 .8],'LabelSpacing',labelspacing);
plot(XX(:,1),-bathy(:,jidx),'k-','LineWidth',2);
contour(XX,-ZZ,squeeze(theta_djf(:,jidx,:)),[-1.9 -1.9],'EdgeColor','w','LineWidth',1.5,'LineStyle','--');
hold off;
caxis([pt_min pt_max]);
colormap(gca,cmocean('thermal'));
axis([-37.16 -30.99 0 dmax])
set(gca,'YDir','reverse');
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);
set(gca,'Position',axpos(3,:));
set(gca,'FontSize',fontsize);
text(-37,600,'Model','FontSize',fontsize+2);




%%% 
%%% OBS, FILCHNER OUTFLOW
%%%

subplot('Position',axpos(4,:));
pcolor(NARE85_HH,NARE85_ZZ,NARE85_pt)
shading interp;
hold on;
[C,h]=contour(NARE85_HH,NARE85_ZZ,NARE85_pd,[27.4:.1:28.4],'EdgeColor',[.8 .8 .8]);
clabel(C,h,'Color',[.8 .8 .8],'LabelSpacing',labelspacing);
plot(NARE85_lon,-NARE85_maxdepth,'k-','LineWidth',2);
contour(NARE85_HH,NARE85_ZZ,NARE85_pt,[-1.9 -1.9],'EdgeColor','w','LineWidth',1.5,'LineStyle','--');
for i=1:size(IC92_HH,1)
  plot([NARE85_HH(i,1) NARE85_HH(i,1)],[0 NARE85_ZZ(i,end)],'--','Color',[1 1 1],'LineWidth',0.5);
end
hold off;
set(gca,'YDir','reverse');
colormap(gca,cmocean('thermal'));
caxis([pt_min pt_max]);
set(gca,'Color',[.8 .8 .8]);
axis([-37.16 -30.99 0 dmax]);
set(gca,'Position',axpos(4,:));
set(gca,'FontSize',fontsize);
text(-37,600,'NARE85','FontSize',fontsize+2);

%%% Colorbar for this row
cbhandle = colorbar;
set(cbhandle,'Position',cb3_pos);








%%% 
%%% MODEL, WESTERN WEDDELL
%%%

pt_min = -1.9;
pt_max = 0.5;

%%% Meshgrid for plots
[ZZ,XX]=meshgrid(zz,XC(:,1));
jidx = find(YC(1,:)>-67.62,1,'first')
YC(i,jidx)
dmax = 3000;

for i=1:Nx
  kbot = find(hFacC(i,jidx,:)>0,1,'last');
  if (~isempty(kbot))
    ZZ(i,kbot) = -sum(hFacC(i,jidx,1:kbot-1).*DRF(1:kbot-1),3)-hFacC(i,jidx,kbot)*DRF(kbot)/2;
    ZZ(i,kbot) = bathy(i,jidx);
  end
end

%%% Make plots
subplot('Position',axpos(5,:));
colormap jet;
pcolor(XX,-ZZ,squeeze(theta_jja(:,jidx,:)))
shading interp
hold on;
[C,h]=contour(XX,-ZZ,squeeze(pd_jja(:,jidx,:)-1000),[27.4:.1:27.7 27.75:.05:27.8 27.85],'EdgeColor',[.8 .8 .8]);
clabel(C,h,'Color',[.8 .8 .8],'LabelSpacing',labelspacing);
plot(XX(:,1),-bathy(:,jidx),'k-','LineWidth',2);
hold off;
caxis([pt_min pt_max]);
colormap(gca,cmocean('thermal'));
axis([-57.52 -52.07 0 dmax])
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);
set(gca,'Position',axpos(5,:));
set(gca,'FontSize',fontsize);
text(-57.4,2800,'Model','FontSize',fontsize+2);





%%%
%%% OBS, WESTERN WEDDELL
%%%

subplot('Position',axpos(6,:));
pcolor(IC92_HH,IC92_ZZ,IC92_pt)
shading interp;
hold on;
[C,h]=contour(IC92_HH,IC92_ZZ,IC92_pd,[27.4:.1:27.7 27.75:.05:28.4],'EdgeColor',[.8 .8 .8]);
clabel(C,h,'Color',[.8 .8 .8],'LabelSpacing',labelspacing);
plot(IC92_lon,-IC92_maxdepth,'k-','LineWidth',2);
for i=1:size(IC92_HH,1)
  plot([IC92_HH(i,1) IC92_HH(i,1)],[0 IC92_ZZ(i,end)],'--','Color',[1 1 1],'LineWidth',0.5);
end
hold off;
set(gca,'YDir','reverse');
xlabel('Longitude')
colormap(gca,cmocean('thermal'));
caxis([pt_min pt_max]);
set(gca,'Color',[.8 .8 .8]);
axis([-57.52 -52.07 0 dmax]);
set(gca,'Position',axpos(6,:));
set(gca,'FontSize',fontsize);
text(-57.4,2800,'IC92','FontSize',fontsize+2);

cbhandle = colorbar;
set(cbhandle,'Position',cb4_pos);










%%%
%%% PLOT LABELS
%%%

annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1) axpos(2,2) 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.05 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.05 axpos(4,2)-0.04 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(5,1)-0.05 axpos(5,2)-0.05 0.03 0.03],'String',axlabels{5},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(6,1)-0.05 axpos(6,2)-0.05 0.03 0.03],'String',axlabels{6},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
