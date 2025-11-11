%%%
%%% animEddiesVsMelt.m
%%%
%%% Makes a movie of the vorticity and the ice shelf melt to visualize their correspondence.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Read experiment data
loadexp;

%%% Vertical grid spacing matrix
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

%%% Diagnostic indix corresponding to instantaneous velocity
diagnum = length(diag_frequency);
diagnum = 1;

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(diagnum);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round(nIter0 + (1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
nDumps = length(dumpIters);

%%% Physical parameters
rho_i = 920;
Sref = 35;

%%% Initialize movie
figure(1);
set(gcf,'Color','w');
set(gcf,'Position',[156         412        2117         840]);
M = moviein(nDumps);


%%% TODO
% dumpStart = 2105280;
% dumpStep = 86400/60;
% nDumps = 366;
% dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;

%%% TODO
% dumpStart = 1578240;
% dumpStep = 43200/60;
% nDumps = 731;
% dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;

% DXC3D = repmat(DXC(2:Nx,:),[1 1 Nr]);
% DYC3D = repmat(DYC(:,2:Ny),[1 1 Nr]);
% vort = zeros(Nx,Ny,Nr);
Omega = 2*pi*366/365/86400;
ff = 2*Omega*sind(YG);
  
%%% Loop through iterations
% for n=88:length(dumpIters)
% for n = 1:88
% for n = 730:length(dumpIters)
% for n = 1:1526
for n = 2120:length(dumpIters)
% for n=250:300
 
  tt(n) =  dumpIters(n)*deltaT/86400;
  tt(n)
  
  %%% Attempt to load either instantaneous velocities or their squares
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),dumpIters(n));
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/SIuice_12hourly'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/SIvice_12hourly'),dumpIters(n));
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_daily'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_daily'),dumpIters(n)); 
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n));      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n)); 
  salt = rdmdsWrapper(fullfile(exppath,'/results/SALT_inst'),dumpIters(n));   
  SHIfwFlx = rdmdsWrapper(fullfile(exppath,'/results/SHIfwFlx_inst'),dumpIters(n)); 
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/U'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/V'),dumpIters(n)); 
  if (isempty(uvel) || isempty(vvel))   
    break;
  end

  
  %%% Plot the vorticity  

  %%% Vorticity on a z-level
  vort = zeros(Nx,Ny);
  zlev = 1; %%% Surface
  % zlev = 39; %%% 500m
%   zlev = 44;
  uvel = uvel(:,:,zlev);
  vvel = vvel(:,:,zlev);
  uvel(hFacW(:,:,zlev)==0) = NaN;
  vvel(hFacS(:,:,zlev)==0) = NaN;
  vort(:,2:Ny) = - (uvel(:,2:Ny)-uvel(:,1:Ny-1))./DYC(:,2:Ny);
  vort(2:Nx,:) = vort(2:Nx,:) + (vvel(2:Nx,:)-vvel(1:Nx-1,:))./DXC(2:Nx,:);
  
  %%% Divergence on a z-level
  % vort = zeros(Nx,Ny);
  % zlev = 39; %%% 500m
  % uvel(hFacW==0) = NaN;
  % vvel(hFacS==0) = NaN;
  % vort(:,1:Ny-1) = (vvel(:,2:Ny,zlev)-vvel(:,1:Ny-1,zlev))./DYG(:,1:Ny-1);
  % vort(1:Nx-1,:) = vort(2:Nx,:) + (uvel(2:Nx,:,zlev)-uvel(1:Nx-1,:,zlev))./DXG(1:Nx-1,:);

  %%% Barotropic vorticity
%   vort = zeros(Nx,Ny);
%   ubt = sum(uvel.*DZ.*hFacW,3) ./ sum(DZ.*hFacW,3);
%   vbt = sum(vvel.*DZ.*hFacS,3) ./ sum(DZ.*hFacS,3);
%   vort(:,2:Ny) = - (ubt(:,2:Ny)-ubt(:,1:Ny-1))./DYC(:,2:Ny);
%   vort = vort + (vbt([2:Nx 1],:)-vbt(:,:))./DXC; 

%   vort = 0*vort;
%   uvel(hFacW==0) = NaN;
%   vvel(hFacS==0) = NaN;  
%   vort(:,2:Ny,:) = - (uvel(:,2:Ny,:)-uvel(:,1:Ny-1,:))./DYC3D;
%   vort(2:Nx,:,:) = vort(2:Nx,:,:) + (vvel(2:Nx,:,:)-vvel(1:Nx-1,:,:))./DXC3D;
%   vort = nanmean(vort,3); %%% Approximate depth-average

  %%% Convert to melt rate
  meltrate = -SHIfwFlx/rho_i*t1year;
  meltrate((hFacC(:,:,1)>0) | (sum(hFacC,3)==0)) = NaN;

  % latMin = YC(1,spongethickness+1);
  latMin = YC(1,1);
%   latMin = -78.5;
  latMax = YC(1,end-spongethickness);
  lonMin = XC(spongethickness+1,1);
%   lonMin = -62;
  lonMax = XC(end-spongethickness,1);

  % lonMin = -24;
  % lonMax = -20;
  % latMin = -74;
  % latMax = -72;
  
  clf;    
  set(gcf,'color','w');
  meltlim = [-10 10];
  % meltlim = [-40 40];
  
  
  %%%%% Plot vorticity %%%%%

  %%% Plotting options
  % clim = [-.4 .4];
  clim = [-1 1];
  bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
  axpos = [0.02 0 0.45 0.9];

  %%% Set up map plot
  subplot('Position',axpos);
  axesm('eqaconicstd',...
    'fontsize',13,...
    'Grid','on', ...    
    'Frame','off', ...
    'MapLatLimit',[latMin latMax], ...
    'MapLonLimit',[lonMin lonMax], ... 
    'MapParallels',[-85 -65], ...
    'PLineLocation', 1, ...
    'MLineLocation', 2,...
    'MeridianLabel', 'on', ...
    'ParallelLabel', 'on');    %, ...
  %           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
  axis off;
  setm(gca,'MLabelParallel',-20)
  
  %%% Plot vorticity in color contours
  pcolorm(YG,XG,vort./abs(ff));
  shading flat
  shading interp
  colormap(gca,cmocean('delta',40));
  caxis(clim);
  
  
  tightmap;

  %%% Add colorbar and title
  h = colorbar;
  set(h,'FontSize',12);
  set(h,'Position',[0.02 0.43 0.01 .12]);
  title(h,'$\zeta/|f|$','Fontsize',14,'interpreter','latex');

  %%% Add bathymetry contours
  hold on;
  thick_plot = SHELFICEtopo-bathy;
  [cs,C] = contourm(YC,XC,thick_plot,bathycntrs,'EdgeColor',[.25 .25 .25]); 
  hh = clabelm(cs,C);
  set(hh,'fontsize',8,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')     
  hold off;
  
  hold on;
  Nlon = 101;
  dLon = (lonMax-lonMin)/(Nlon-1);
  plotm([latMin*ones(1,Nlon) latMax*ones(1,Nlon) latMin],[lonMin:dLon:lonMax lonMax:-dLon:lonMin lonMin],'k-','LineWidth',1);
  hold off
          
  %%% Add axis labels
  xlabel('Longitude','interpreter','latex');
  ylabel('Latitude','interpreter','latex');
  set(gca,'Position',axpos);
  ax = gca;
  








  %%% Add plot of ice shelf melt rate
  subplot('Position',[0 0 0.01 0.01])
  ax2 = axesm('eqaconicstd',...
    'fontsize',13,...
    'Grid','on', ...    
    'Frame','off', ...
    'MapLatLimit',[latMin latMax], ...
    'MapLonLimit',[lonMin lonMax], ... 
    'MapParallels',[-85 -65], ...
    'PLineLocation', 0, ...
    'MLineLocation', 0,...
    'MeridianLabel', 'on', ...
    'ParallelLabel', 'on');   
  axis off;
  setm(ax2,'MLabelParallel',-20)

  meltrate(meltrate==0) = NaN;
  pcolorm(YC,XC,meltrate);
  shading flat
  colormap(ax2,cmocean('balance'));

  tightmap;

   %%% Add colorbar and title
  h = colorbar;
  set(h,'FontSize',12);
  set(h,'Position',[0.48 0.43 0.01 .12]);
  title(h,'Melt rate (m/yr)','Fontsize',14,'interpreter','latex');

  hold on;
  Nlon = 101;
  dLon = (lonMax-lonMin)/(Nlon-1);
  plotm([latMin*ones(1,Nlon) latMax*ones(1,Nlon) latMin],[lonMin:dLon:lonMax lonMax:-dLon:lonMin lonMin],'k-','LineWidth',1);
  hold off

  drawnow;
  set(ax2,'Position',get(ax,'Position'));
  caxis(meltlim);
  drawnow;











 %%
  
  %%%%% Plot subsurface salinity %%%%%
  
  salt_plot = sum(salt(:,:,1:21).*hFacC(:,:,1:21).*DRF(:,:,1:21),3)./sum(hFacC(:,:,1:21).*DRF(:,:,1:21),3);
  salt_plot(hFacC(:,:,1)==0) = NaN;


  %%% Plotting options
  clim = [33.4 34.65];
  bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
  axpos = [0.52 0 0.45 0.9];

  %%% Set up map plot
  subplot('Position',axpos);
  axesm('eqaconicstd',...
    'fontsize',13,...
    'Grid','on', ...    
    'Frame','off', ...
    'MapLatLimit',[latMin latMax], ...
    'MapLonLimit',[lonMin lonMax], ... 
    'MapParallels',[-85 -65], ...
    'PLineLocation', 1, ...
    'MLineLocation', 2,...
    'MeridianLabel', 'on', ...
    'ParallelLabel', 'on');    %, ...
  %           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
  axis off;
  setm(gca,'MLabelParallel',-20)
  
  %%% Plot vorticity in color contours
  pcolorm(YG,XG,salt_plot);
  shading interp
  colormap(gca,cmocean('haline',40));
  caxis(clim);
  
  
  tightmap;

  %%% Add colorbar and title
  h = colorbar;
  set(h,'FontSize',12);
  set(h,'Position',[0.95 0.43 0.01 .12]);
  title(h,'$S$ (g/kg)','Fontsize',14,'interpreter','latex');

  %%% Add bathymetry contours
  hold on;
  thick_plot = SHELFICEtopo-bathy;
  [cs,C] = contourm(YC,XC,thick_plot,bathycntrs,'EdgeColor',[.25 .25 .25]); 
  hh = clabelm(cs,C);
  set(hh,'fontsize',8,'Color',[.25 .25 .25],'BackgroundColor','none','Edgecolor','none')     
  hold off;
  
  hold on;
  Nlon = 101;
  dLon = (lonMax-lonMin)/(Nlon-1);
  plotm([latMin*ones(1,Nlon) latMax*ones(1,Nlon) latMin],[lonMin:dLon:lonMax lonMax:-dLon:lonMin lonMin],'k-','LineWidth',1);
  hold off
          
  %%% Add axis labels
  xlabel('Longitude','interpreter','latex');
  ylabel('Latitude','interpreter','latex');
  set(gca,'Position',axpos);
  ax = gca;
  








  %%% Add plot of ice shelf melt rate
  subplot('Position',[0 0 0.01 0.01])
  ax2 = axesm('eqaconicstd',...
    'fontsize',13,...
    'Grid','on', ...    
    'Frame','off', ...
    'MapLatLimit',[latMin latMax], ...
    'MapLonLimit',[lonMin lonMax], ... 
    'MapParallels',[-85 -65], ...
    'PLineLocation', 0, ...
    'MLineLocation', 0,...
    'MeridianLabel', 'on', ...
    'ParallelLabel', 'on');   
  axis off;
  setm(ax2,'MLabelParallel',-20)

  meltrate(meltrate==0) = NaN;
  pcolorm(YC,XC,meltrate);
  shading flat
  colormap(ax2,cmocean('balance'));

  tightmap;

  hold on;
  Nlon = 101;
  dLon = (lonMax-lonMin)/(Nlon-1);
  plotm([latMin*ones(1,Nlon) latMax*ones(1,Nlon) latMin],[lonMin:dLon:lonMax lonMax:-dLon:lonMin lonMin],'k-','LineWidth',1);
  hold off

  drawnow;
  set(ax2,'Position',get(ax,'Position'));
  caxis(meltlim);
  drawnow;

  thedate = datenum('2008-01-00')+floor(tt(n));
  annotation('TextBox',[0.46 0.93 0.2 0.05],'EdgeColor','None','String',datestr(thedate),'interpreter','latex','FontSize',16);


  M(n) = getframe(gcf);  
   
end