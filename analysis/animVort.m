%%%
%%% animVort.m
%%%
%%% Makes a movie of the vorticity.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Read experiment data
% loadexp;

%%% Vertical grid spacing matrix
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

%%% Diagnostic indix corresponding to instantaneous velocity
diagnum = length(diag_frequency);
% diagnum = 1;

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(diagnum);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round(nIter0 + (1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
nDumps = length(dumpIters);

%%% Initialize movie
figure(1);
set(gcf,'Color','w');
M = moviein(nDumps);

%%% TODO
% dumpStart = 2105280;
% dumpStep = 86400/60;
% nDumps = 366;
% dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;

%%% TODO
dumpStart = 1578240;
dumpStep = 43200/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;

% DXC3D = repmat(DXC(2:Nx,:),[1 1 Nr]);
% DYC3D = repmat(DYC(:,2:Ny),[1 1 Nr]);
% vort = zeros(Nx,Ny,Nr);
Omega = 2*pi*366/365/86400;
ff = 2*Omega*sind(YG);
  
%%% Loop through iterations
for n=1:length(dumpIters)
% for n=250:300
 
  tt(n) =  dumpIters(n)*deltaT/86400;
  tt(n)
  
  %%% Attempt to load either instantaneous velocities or their squares
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),dumpIters(n));
  uvel = rdmdsWrapper(fullfile(exppath,'/results/SIuice_12hourly'),dumpIters(n)) ;      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/SIvice_12hourly'),dumpIters(n));
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_daily'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_daily'),dumpIters(n)); 
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),dumpIters(n));      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));    
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/U'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/V'),dumpIters(n)); 
  if (isempty(uvel) || isempty(vvel))   
    break;
  end
  
  %%% Plot the vorticity  

  %%% Vorticity on a z-level
  vort = zeros(Nx,Ny);
  zlev = 1;
  uvel(hFacW(:,:,zlev)==0) = NaN;
  vvel(hFacS(:,:,zlev)==0) = NaN;
  vort(:,2:Ny) = - (uvel(:,2:Ny,zlev)-uvel(:,1:Ny-1,zlev))./DYC(:,2:Ny);
  vort(2:Nx,:) = vort(2:Nx,:) + (vvel(2:Nx,:,zlev)-vvel(1:Nx-1,:,zlev))./DXC(2:Nx,:);
  
  %%% Divergence on a z-level
%   vort = zeros(Nx,Ny);
%   zlev = 1;
%   uvel(hFacW==0) = NaN;
%   vvel(hFacS==0) = NaN;
%   vort(:,1:Ny-1) = (vvel(:,2:Ny,zlev)-vvel(:,1:Ny-1,zlev))./DYG(:,1:Ny-1);
%   vort(1:Nx-1,:) = vort(2:Nx,:) + (uvel(2:Nx,:,zlev)-uvel(1:Nx-1,:,zlev))./DXG(1:Nx-1,:);

%   %%% Barotropic vorticity
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

%   latMin = min(min(YC));
  latMin = -78.5;
  latMax = YC(1,end-spongethickness);
%   lonMin = min(min(XC));
  lonMin = -62;
  lonMax = XC(end-spongethickness,1);
  
  clf;    
  set(gcf,'color','w');
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
  pcolorm(YG,XG,vort./abs(ff));
  shading interp
  colormap(cmocean('balance'));

  hold on;
  [cs,C]=contourm(YC,XC,SHELFICEtopo-bathy,[100 500 1000 2000 3000 4000],'EdgeColor','k');
  chandle = clabelm(cs,C);
  set(chandle,'fontsize',12,'Color','k','BackgroundColor','none','Edgecolor','none') 
  hold off;

  %%% Add colorbar and title
  h = colorbar;  
  set(h,'Position',[0.92 0.33 0.01 .26])  
  caxis([-.3 .3]);
  
  set(gca,'FontSize',16);  
  thedate = datenum('2008-01-01')+floor(tt(n));
  h = title(['\zeta/|f| on ',datestr(thedate)]);
%   h = title(['\zeta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
%   h = title(['\delta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
  %  h = title(['\zeta/f at t= ',num2str(tt(n)/365,'%.1f'),' years']);  
  set(gca,'Position',[0 0 0.95 0.95])
  tightmap;
  drawnow;
  set(h,'Position',get(h,'Position')+[0 +0.0075 0])


  M(n) = getframe(gcf);  
   
end