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
diagnum = 69;

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
figure(2);
set(gcf,'Color','w');
M = moviein(nDumps);
icecolor = [186 242 239]/255;

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

DXC3D = repmat(DXC(2:Nx,:),[1 1 Nr]);
DYC3D = repmat(DYC(:,2:Ny),[1 1 Nr]);
% vort = zeros(Nx,Ny,Nr);
Omega = 2*pi*366/365/86400;
ff = 2*Omega*sind(YG);
g = 9.81;
  
%%% Loop through iterations
% for n=88:length(dumpIters)
% for n = 1:88
% for n = 60:127
for n=1:length(dumpIters)
 
  tt(n) =  dumpIters(n)*deltaT/86400;
  tt(n)
  
  %%% Attempt to load either instantaneous velocities or their squares
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),dumpIters(n)) ;      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),dumpIters(n));
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/SIuice_12hourly'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/SIvice_12hourly'),dumpIters(n));
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_daily'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_daily'),dumpIters(n)); 
  % uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n));      
  % vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n));    
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/U'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/V'),dumpIters(n)); 
  % eta = rdmdsWrapper(fullfile(exppath,'/results/ETAN_inst'),dumpIters(n));
  % eta = rdmdsWrapper(fullfile(exppath,'/results/SIarea_inst'),dumpIters(n));
  % eta = rdmdsWrapper(fullfile(exppath,'/results/ETAN'),dumpIters(n));
  % if (isempty(uvel) || isempty(vvel))   
  %   break;
  % end
  
  %%% Plot the vorticity  

  %%% Vorticity on a z-level
%   vort = zeros(Nx,Ny);
%   zlev = 1;
% %   zlev = 25;
%   % zlev = 44;
%   uvel = uvel(:,:,zlev);
%   vvel = vvel(:,:,zlev);
%   uvel(hFacW(:,:,zlev)==0) = NaN;
%   vvel(hFacS(:,:,zlev)==0) = NaN;
%   vort(:,2:Ny) = - (uvel(:,2:Ny)-uvel(:,1:Ny-1))./DYC(:,2:Ny);
%   vort(2:Nx,:) = vort(2:Nx,:) + (vvel(2:Nx,:)-vvel(1:Nx-1,:))./DXC(2:Nx,:);


  %%% Geostrophic relative vorticity
  % eta(hFacC(:,:,1)==0) = NaN;
  % vvel = 0*eta;
  % uvel = 0*eta;
  % vvel(2:end,1:end) = g * (diff(eta,1,1) ./ DXC(2:end,1:end)) ./ ff(2:end,1:end);
  % uvel(1:end,2:end) = -g * diff(eta,1,2) ./ DYC(1:end,2:end) ./ ff(1:end,2:end);   
  % vort = 0*eta;
  % vort(2:end-1,2:end-1) = (vvel(3:end,2:end-1) - vvel(2:end-1,2:end-1)) ./ DXG(2:end-1,2:end-1) ...
  %                       - (uvel(2:end-1,3:end) - uvel(2:end-1,2:end-1)) ./ DYG(2:end-1,2:end-1);
    
  
  %%% Okubo-Weiss on a z-level
 
  % zlev = 1;
  % uvel = uvel(:,:,zlev);
  % vvel = vvel(:,:,zlev);
  % uvel(hFacW(:,:,zlev)==0) = NaN;
  % vvel(hFacS(:,:,zlev)==0) = NaN;
  % dudx = (uvel([2:Nx 1],:) - uvel(1:Nx,:)) ./ DXG; %%% du/dx on cell centers
  % dvdy = (vvel(:,[2:Ny 1]) - vvel(:,1:Ny)) ./ DYG; %%% dv/dy on cell centers  
  % dvdx = (vvel(1:Nx,:) - vvel([Nx 1:Nx-1],:)) ./ DXC; %%% dv/dx on cell corners
  % dudy = (uvel(:,1:Ny) - uvel(:,[Ny 1:Ny-1])) ./ DYC; %%% du/du on cell corners
  % Sn = dudx - dvdy;
  % Ss = dudy + dvdx;
  % omega = dvdx - dudy;
  % OW = Ss.^2 - omega.^2 + 0.25*(Sn(1:Nx,1:Ny).^2 + Sn([Nx 1:Nx-1],1:Ny).^2 + Sn(1:Nx,[Ny 1:Ny-1]).^2 + Sn([Nx 1:Nx-1],[Ny 1:Ny-1]).^2);
  
  
  %%% Divergence on a z-level
%   vort = zeros(Nx,Ny);
%   zlev = 44;
%   uvel(hFacW==0) = NaN;
%   vvel(hFacS==0) = NaN;
%   vort(:,1:Ny-1) = (vvel(:,2:Ny,zlev)-vvel(:,1:Ny-1,zlev))./DYG(:,1:Ny-1);
%   vort(1:Nx-1,:) = vort(2:Nx,:) + (uvel(2:Nx,:,zlev)-uvel(1:Nx-1,:,zlev))./DXG(1:Nx-1,:);

  %%% Barotropic vorticity
  % vort = zeros(Nx,Ny);
  % ubt = sum(uvel.*DZ.*hFacW,3) ./ sum(DZ.*hFacW,3);
  % vbt = sum(vvel.*DZ.*hFacS,3) ./ sum(DZ.*hFacS,3);
  % vort(:,2:Ny) = - (ubt(:,2:Ny)-ubt(:,1:Ny-1))./DYC(:,2:Ny);
  % vort = vort + (vbt([2:Nx 1],:)-vbt(:,:))./DXC; 

  %%% Depth-averaged vorticity
  vort = 0*uvel;
  uvel(hFacW==0) = NaN;
  vvel(hFacS==0) = NaN;  
  vort(:,2:Ny,:) = - (uvel(:,2:Ny,:)-uvel(:,1:Ny-1,:))./DYC3D;
  vort(2:Nx,:,:) = vort(2:Nx,:,:) + (vvel(2:Nx,:,:)-vvel(1:Nx-1,:,:))./DXC3D;
  vort = nansum(vort.*DRF,3) ./ nansum(~isnan(vort).*DRF,3); %%% Approximate depth-average

  %%% Water column thickness to plot
  thick_plot =SHELFICEtopo-bathy;
  % thick_plot(sum(hFacC,3)==0) = NaN;
  % vort(sum(hFacC,3)==0) = NaN;

  % latMin = min(min(YC));
  % latMin = YC(1,spongethickness+1);
  % latMin = -75.5;
  latMin = -78.5;
  latMax = YC(1,end-spongethickness);
  % lonMin = min(min(XC));
  % latMax = -72;
  % lonMin = XC(spongethickness+1,1);
  % lonMin = -44;
  lonMin = -62;
  lonMax = XC(end-spongethickness,1);
  % lonMax = -28;

  % % Desired geographic window (what you *mean*)
  % latlim = [latMin latMax];
  % lonlim = [lonMin lonMax];
  % 
  % % Choose map "center" and desired clockwise rotation
  % lat0 = mean(latlim);
  % lon0 = mean(lonlim);
  % az   = 30;                 % clockwise, degrees
  
  clf;    
  set(gcf,'color','w');
  ax = axesm('eqaconicstd',...
  'fontsize',18,...
  'Grid','on', ...    
  'Frame','off', ...  
  'MapParallels',[-85 -65], ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ... 
  'PLineLocation', 1, ...
  'MLineLocation', 4,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
            % 'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
  % 'Origin', [0.5*(latMin+latMax) 0.5*(lonMin+lonMax) 0], ...


  axis off;
  setm(gca,'MLabelParallel',-20)


  pcolorm(YG,XG,vort./abs(ff));
  % pcolorm(YG,XG,OW./ff.^2);
  shading interp
  colormap(cmocean('balance'));

  

  hold on;
  
  % [cs,C]=contourm(YC,XC,thick_plot,[100 500 1000 2000 3000 4000],'EdgeColor','k');
  [cs,C]=contourm(YC,XC,thick_plot,[500 1000 2000 3000 4000],'EdgeColor','k');
  chandle = clabelm(cs,C);
  set(chandle,'fontsize',14,'Color','k','BackgroundColor','none','Edgecolor','none') 
  hold off;

  %%% Add colorbar and title
  h = colorbar;  
  set(h,'Position',[0.92 0.33 0.01 .26])  
  % caxis([-.8 .8]);
  caxis([-.4 .4]);
  
  set(gca,'FontSize',18);    
  tightmap;
  drawnow;
  set(gca,'Position',[0.04 0 0.89 0.9])  




  %%% Add Bathymetry
  axpos = get(gca,'Position');
  subplot('Position',[0 0 0.01 0.01])  
  ax3 = axesm('eqaconicstd',...
  'fontsize',18,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 0, ...
  'MLineLocation', 0,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
  axis off;
  setm(ax3,'MLabelParallel',-20)
  set(ax3,'Position',axpos);
  
  land_plot = ones(Nx,Ny);
  land_plot(~isnan(vort)) = NaN;
  pcolorm(YG,XG,land_plot);
  shading interp
  colormap(ax3,[.8 .8 .8]);
  tightmap;
  drawnow;
  set(ax3,'Position',axpos+[0 0 0 0]);




  %%% Add ice shelves
  axpos = get(gca,'Position');
  subplot('Position',[0 0 0.01 0.01])  
  ax2 = axesm('eqaconicstd',...
  'fontsize',18,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 0, ...
  'MLineLocation', 0,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
  axis off;
  setm(ax2,'MLabelParallel',-20)
  set(ax2,'Position',axpos);
  
  ice_plot = SHELFICEtopo;
  ice_plot((ice_plot == 0) | (SHELFICEtopo==bathy)) = NaN;
  pcolorm(YG,XG,ice_plot);
  shading interp
  colormap(ax2,icecolor);
  tightmap;
  drawnow;
  set(ax2,'Position',axpos+[0 0 0 0]);








  thedate = datenum('2008-01-00')+floor(tt(n));
  %%% Add title
    h = title(['\zeta/|f| on ',datestr(thedate)],'FontSize',18);
%   h = title(['\zeta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
%   h = title(['\delta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
  %  h = title(['\zeta/f at t= ',num2str(tt(n)/365,'%.1f'),' years']);  
  set(h,'Position',get(h,'Position')+[0 0.01 0])



  M(n) = getframe(gcf);  
   
end