
%%%
%%% plotBaroSF
%%%
%%% Plots time-mean barotropic streamfunction.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_RTOPO2';
loadexp;





%%% Set up the figure
figure(201)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34   681   926]);
% fontsize = 18;














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DOMAIN BATHYMETRY PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load pre-computed streamfunction
load(fullfile('products',[expname,'_PsiBT.mat']));

%%% Plotting options
latMin = min(min(YC));
latMax = YC(1,end);
lonMin = min(min(XC));
lonMax = XC(end,1);
clim = [0 5000];

%%% Set up map plot
subplot(3,1,1)
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
pcolorm(YC,XC,bathy_plot);

msk_iceshelf = 0*bathy;
msk_iceshelf((SHELFICEtopo>=0) | (SHELFICEtopo<=bathy)) = NaN;
edge_cavity = 0*bathy;
edge_cavity((SHELFICEtopo>=0) | (SHELFICEtopo<=bathy)) = 1;
msk_obcs = 5000*ones(size(bathy));
msk_obcs((YC<YC(1,end-spongethickness+1)) & (XC<XC(end-spongethickness+1,1))) = NaN;
msk_obcs(SHELFICEtopo<=bathy) = NaN;
edge_obcs = ones(size(bathy));
edge_obcs((YC<YC(1,end-spongethickness+1)) & (XC<XC(end-spongethickness+1,1))) = 0;
edge_obcs(SHELFICEtopo<=bathy) = 0;
hold on;
pcolorm(YC,XC,msk_iceshelf,'FaceAlpha',0.3);
contourm(YC,XC,edge_cavity,[0.5 0.5],'EdgeColor',[1 1 1],'LineWidth',1.5);
pcolorm(YC,XC,msk_obcs,'FaceAlpha',0.5);
contourm(YC,XC,edge_obcs,[0.5 0.5],'EdgeColor',[0 0 0],'LineWidth',1.5);
hold off;

%%% Add colorbar and title
h = colorbar;
caxis(clim);
colormap(flip(haxby(50)));
set(gca,'FontSize',10);
set(h,'Position',[0.92 0.68 0.02 .29])
tightmap;
title(h,'m','Fontsize',14,'interpreter','latex');
set(gca,'Position',[0.07 0.68 .85 .26]);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DOMAIN DEPTH-MEAN SPEED PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load pre-computed streamfunction
load(fullfile('products',[expname,'_PsiBT.mat']));


%%% Plotting options
latMin = min(min(YC));
latMax = YC(1,end-spongethickness);
lonMin = min(min(XC));
lonMax = XC(end-spongethickness,1);
clim = [0 0.25];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];

%%% Set up map plot
subplot(3,1,2)
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

%%% Plot streamfunction in color contours
uu_zavg = sum(uu.*hFacW.*DRF,3)./sum(hFacW.*DRF,3);
vv_zavg = sum(vv.*hFacS.*DRF,3)./sum(hFacS.*DRF,3);
pcolorm(YC,XC,sqrt(uu_zavg.^2+vv_zavg.^2));
shading interp
colormap(gca,hot(50));
caxis(clim);

%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.92 0.33 0.02 .26])
tightmap;
title(h,'m/s','Fontsize',14,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.75 .75 .75]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.75 .75 .75],'BackgroundColor','none','Edgecolor','none')       
        
%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',[0.07 0.33 .85 .29]);
hold off







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CAVITY STREAMFUNCTION PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load pre-computed streamfunction
load(fullfile('products',[expname,'_PsiBT.mat']));

%%% Adjust grids/streamfunction to focus on cavity
rem = YC(1,:)>-72;
YG(:,rem)=[];
XG(:,rem)=[];
Psi(:,rem)=[];
bathy(:,rem)=[];
SHELFICEtopo(:,rem)=[];   
YC(:,rem)=[];
XC(:,rem)=[]; 
Psi_plot = Psi/1e6;
Psi_plot(bathy>=SHELFICEtopo)=NaN;
Psi_plot(Psi_plot==0)=NaN;

%%% Plotting options
latMin = min(min(YC));
latMax = -72;
lonMin = -80;
lonMax = -30;
clim = [-3 3];
colorcntrs = [-3:0.25:3];
linecntrs = [-3:0.25:3];
bathycntrs = [0 250 500 1000 2000 3000 4000];

%%% Set up map plot
subplot(3,1,3);
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



%%% Plot streamfunction in color contours
contourfm(YG,XG,Psi_plot,colorcntrs,'EdgeColor','none');
set(gca,'Position',[0.07 0.01 .85 .26]);
colormap(gca,redblue(50));

%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.92 0.01 0.02 .23])
tightmap;
title(h,'Sv','Fontsize',14,'interpreter','latex');

%%% Add contour lines on selected streamlines
hold on
contourm(YG,XG,Psi_plot,linecntrs,'EdgeColor',[.5,.5,.5],'LineWidth',.6)


%%% Add bathymetry contours
% bathy_plot = bathy;
% bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
% [cs,C] = contourm(YC,XC,bathy_plot,bathycntrs,'EdgeColor','black'); 
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor','black'); 
hh = clabelm(cs,C);
set (hh,'fontsize',8,'BackgroundColor','none','Edgecolor','none')
         
hold on
        
% contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
        
caxis(clim)

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
h = title('Barotropic Stream Function (Sv) (FRIS Cavity)','interpreter','latex','Fontsize',18);
set(h,'Position',[0 .9 0])

hold off









annotation('textbox',[0.02 0.65 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.35 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.05 0.05 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

