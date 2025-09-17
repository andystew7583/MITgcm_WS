%%%
%%% paper3_plotEKEprod.m
%%%
%%% Plots EKE and its production rates
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;

%%% Load shelf heat budget diagnostics
outfname = [expname,'_ShelfHeatBudget.mat'];
load(fullfile('products',outfname),'usq_eddy_int','vsq_eddy_int','PEtoEKE_int','MKEtoEKE_int');



%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.2 0.43 .54 .54];
axpos(2,:) = [0.01 0.01 .5 .45];
axpos(3,:) = [0.47 0.01 .55 .45];
cbpos = zeros(2,4);
cbpos(1,:) = [0.7 0.56 0.01 .1];
cbpos(2,:) = [0.46 0.1 0.01 .1];
axlabels = {'\textbf{A}','\textbf{B}','\textbf{C}'};
rho0 = 1027;
Cp = 4000;
colororder = get(gca,'ColorOrder');
linewidth = 1.5;

%%% Plotting range for salinity figure
latMin_b = min(min(YC));
latMax_b = YC(1,end-spongethickness);
lonMin_b = min(min(XC));
lonMax_b = XC(end-spongethickness,1);







%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(208)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34  1000  926]);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EDDY KINETIC ENERGY MAP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EKE_plot = mean(0.5*(usq_eddy_int+vsq_eddy_int),3) ./ sum(hFacC.*DRF,3);
EKE_plot(sum(hFacC,3)==0) = NaN;

%%% Plotting options
clim = [0 1e-2];
cmap = cmocean('amp',50);
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];


%%% Set up map plot
axes('Position',axpos(1,:)+[.1 0 0 0]);
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

%%% Plot surface heat flux
pcolorm(YC,XC,(EKE_plot));
shading interp
caxis(clim);
colormap(gca,cmap);


%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',cbpos(1,:))
tightmap;
title(h,'m$^2$/s$^2$','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.05 .05 .05],'BackgroundColor','none','Edgecolor','none')       
% Nlon = 101;
% dLon = (lonMax_c-lonMin_c)/(Nlon-1);
% plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'r--','LineWidth',2);
hold off;

%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos(1,:));
hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BAROCLINIC PRODUCTION MAP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PEtoEKE_plot = mean(PEtoEKE_int,3) ./ sum(hFacC.*DRF,3);
PEtoEKE_plot(sum(hFacC,3)==0) = NaN;

%%% Plotting options
clim = [-1 1];
cmap = cmocean('balance',50);
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];


%%% Set up map plot
axes('Position',axpos(2,:)+[.1 0 0 0]);
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

%%% Plot surface heat flux
pcolorm(YC,XC,PEtoEKE_plot*1e7);
shading interp
caxis(clim);
colormap(gca,cmap);


%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',cbpos(2,:))
tightmap;
title(h,'10$^{-7}$ m$^2$/s$^3$','Fontsize',fontsize,'interpreter','latex');

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.05 .05 .05],'BackgroundColor','none','Edgecolor','none')       
% Nlon = 101;
% dLon = (lonMax_c-lonMin_c)/(Nlon-1);
% plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'r--','LineWidth',2);
hold off;

%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos(2,:));
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BAROTROPIC PRODUCTION MAP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MKEtoEKE_plot = mean(MKEtoEKE_int,3) ./ sum(hFacC.*DRF,3);
MKEtoEKE_plot(sum(hFacC,3)==0) = NaN;

%%% Plotting options
clim = [-1 1];
cmap = cmocean('balance',50);
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];


%%% Set up map plot
axes('Position',axpos(3,:)+[.1 0 0 0]);
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

%%% Plot surface heat flux
pcolorm(YC,XC,MKEtoEKE_plot*1e7);
shading interp
caxis(clim);
colormap(gca,cmap);

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.05 .05 .05],'BackgroundColor','none','Edgecolor','none')       
% Nlon = 101;
% dLon = (lonMax_c-lonMin_c)/(Nlon-1);
% plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'r--','LineWidth',2);
hold off;

%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos(3,:));
hold off





%%% Add panel labels
annotation('textbox',[axpos(1,1)+0.06 axpos(1,2)+0.08 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)+0.06 axpos(2,2)+0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)+0.06 axpos(3,2)+0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

figure1 = gcf;

% Create textbox
annotation(figure1,'textbox',...
  [0.65 0.47038228941685 0.392 0.0242980561555076],...
  'String',{'Depth-averaged barotropic production'},...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
  [0.15 0.47038228941685 0.33 0.0242980561555073],...
  'String',{'Depth-averaged baroclinic production'},...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
  [0.432 0.959583153347735 0.2685 0.0242980561555075],...
  'String',{'Depth-averaged eddy kinetic energy'},...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none');


