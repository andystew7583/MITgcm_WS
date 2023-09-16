
%%%
%%% plotModelSetup
%%%
%%% Plots the model domain showing grid size and an instantaneous salinity snapshot.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;





%%% Set up the figure
figure(201)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34   681   926]);
% fontsize = 18;
fontsize = 12;













%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GRID SPACING PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load pre-computed streamfunction
load(fullfile('products',[expname,'_PsiBT.mat']));

%%% Plotting options
latMin = min(min(YC));
latMax = YC(1,end);
lonMin = min(min(XC));
lonMax = XC(end,1);
% clim = [min(min(DXG)) max(max(DXG))];
% clim = [.5 1.6];
clim = [0 5000];
gridcntrs = [0.5:.1:1.6];

%%% Set up map plot
subplot(2,1,1)
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

gridsize_plot = max(DXG,DYG)/1000;
gridsize_plot(SHELFICEtopo-bathy<=0) = NaN;
% pcolorm(YC,XC,gridsize_plot);

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

%%% Add grid size contours
hold on;
[cs,C] = contourm(YC,XC,gridsize_plot,gridcntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')       

%%% Add colorbar and title
h = colorbar;
caxis(clim);
colormap(flip(haxby(50)));
% colormap(cmocean('amp',11));
set(gca,'FontSize',10);
set(h,'Position',[0.92 0.55 0.02 .4])
tightmap;
title(h,'m','Fontsize',14,'interpreter','latex');
set(gca,'Position',[0.07 0.55 .85 .4]);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DOMAIN DEPTH-MEAN SPEED PANEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load pre-computed streamfunction
% load(fullfile('products',[expname,'_PsiBT.mat']));
iter = 1985760;
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

%%% Plotting options
latMin = min(min(YC));
latMax = YC(1,end-spongethickness);
lonMin = min(min(XC));
lonMax = XC(end-spongethickness,1);
clim = [34.2 35];
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];

%%% Set up map plot
subplot(2,1,2)
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
% uu_zavg = sum(uu.*hFacW.*DRF,3)./sum(hFacW.*DRF,3);
% vv_zavg = sum(vv.*hFacS.*DRF,3)./sum(hFacS.*DRF,3);
% pcolorm(YC,XC,sqrt(uu_zavg.^2+vv_zavg.^2));
pcolorm(YC,XC,FF);
shading interp
colormap(gca,cmocean('haline',16));
caxis(clim);

%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.92 0.05 0.02 .4]);
tightmap;
title(h,'g/kg','Fontsize',14,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')       
        
%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',[0.07 0.05 .85 .4]);
hold off








annotation('textbox',[0.02 0.55 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.05 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


