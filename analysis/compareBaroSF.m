%%%
%%% plotBaroSF
%%%
%%% Plots time-mean barotropic streamfunction.
%%%

%%% Load experiment
expdir = '../experiments';
expname_tides = 'hires_seq_onethird_RTOPO2';
expname_notides = 'hires_seq_onethird_notides_RTOPO2';
expname = expname_tides;
loadexp; %%% To get grids

%%% Options
show_cavity = true;


%%% Load pre-computed streamfunction
load(fullfile('products',[expname_tides,'_PsiBT.mat']));
Psi_tides = Psi;
load(fullfile('products',[expname_notides,'_PsiBT.mat']));
Psi_notides = Psi;


%%% Modify streamfunction if need be, depending on domain to be plotted
if (show_cavity)
  
  %%% if wanting to look at ice shelf cavity in particular -->
  rem = YC(1,:)>-70;

  YG(:,rem)=[];
  XG(:,rem)=[];
  Psi_tides(:,rem)=[];
  Psi_notides(:,rem)=[];
  bathy(:,rem)=[];
  SHELFICEtopo(:,rem)=[];    
  YC(:,rem)=[];
  XC(:,rem)=[]; 
  
  Psi_tides_plot = Psi_tides/1e6;
  Psi_tides_plot(bathy>=SHELFICEtopo)=NaN;
  Psi_tides_plot(Psi_tides_plot==0)=NaN;
  Psi_notides_plot = Psi_notides/1e6;
  Psi_notides_plot(bathy>=SHELFICEtopo)=NaN;
  Psi_notides_plot(Psi_notides_plot==0)=NaN;
  
else
  
  %%% Full domain
  Psi_tides_plot = Psi_tides/1e6;
  Psi_notides_plot = Psi_notides/1e6;
  
end


%%% Set up the figure
figure(2)
clf
scrsz = get(0,'ScreenSize');
% fontsize = 18;





%%% Set plotting range based on whether we need to show cavity only or
%%% whole domain
if (show_cavity)
  latMin = min(min(YC));
  latMax = -70;
  lonMin = -80;
  lonMax = -20;
  clim = [-1 1];
  colorcntrs = [-1:0.1:1];
  linecntrs = [-1:0.1:1];
  bathycntrs = [0 250 500 1000 2000 3000 4000];
else 
  latMin = min(min(YC));
  latMax = max(max(YC));
  lonMin = min(min(XC));
  lonMax = max(max(XC));
  clim = [-5 5];
  colorcntrs = [-4.75:.5:4.75];
  linecntrs = [-4.5:1:4.5];
%   bathycntrs = [-5000:1000:-1000 -500 -200 -100];
  bathycntrs = [0 250 500 1000 2000 3000 4000];
end



%%% Set up map plot
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
contourfm(YG,XG,Psi_tides_plot-Psi_notides_plot,colorcntrs,'EdgeColor','none');
set(gca,'Position',[0.07 0.03 .85 .92]);
colormap redblue(50);

%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.92 0.1 0.02 .8])
tightmap;
title(h,'Sv','Fontsize',20,'interpreter','latex');

%%% Add contour lines on selected streamlines
hold on
contourm(YG,XG,Psi_tides_plot-Psi_notides_plot,linecntrs,'EdgeColor',[.5,.5,.5],'LineWidth',.6)
hold on


%%% Add bathymetry contours
% bathy_plot = bathy;
% bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
% [cs,C] = contourm(YC,XC,bathy_plot,bathycntrs,'EdgeColor','black'); 
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor','black'); 
hh = clabelm(cs,C);
set (hh,'fontsize',10,'BackgroundColor','none','Edgecolor','none')
         
hold on
        
% contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
        
caxis(clim)

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
h = title('Barotropic Stream Function (Sv) (FRIS Cavity)','interpreter','latex','Fontsize',18);
set(h,'Position',[0 .9 0])

hold on



%%%%%%%%%%%%%%%%

% contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
