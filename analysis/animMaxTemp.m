%%%
%%% animMaxTemp.m
%%%
%%% Makes a movie of the maximum temperature in the water column.
%%%
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

g = 9.81;
  
%%% Loop through iterations
% for n=88:length(dumpIters)
% for n = 1:88
% for n = 60:127
for n=1:length(dumpIters)
 
  tt(n) =  dumpIters(n)*deltaT/86400;
  tt(n)
  
  %%% Attempt to load either instantaneous velocities or their squares
  theta = rdmdsWrapper(fullfile(exppath,'/results/THETA_12hourly'),dumpIters(n)) ; 
  theta(hFacC==0) = NaN;
  thetamax = max(theta(:,:,16:end),[],3);


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


  pcolorm(YC,XC,thetamax);  
  shading interp
  cmap1 = cmocean('thermal',50);
  cmap1 = cmap1(1:40,:);
  cmap2 = cmocean('thermal',400);
  cmap2 = cmap2(321:400,:);
  cmap = [cmap1 ; cmap2];
  colormap(ax,cmap);

  

  hold on;
  
  [cs,C]=contourm(YC,XC,thick_plot,[500 1000 2000 3000 4000],'EdgeColor','k');
  chandle = clabelm(cs,C);
  set(chandle,'fontsize',14,'Color','k','BackgroundColor','none','Edgecolor','none') 
  hold off;

  %%% Add colorbar and title
  h = colorbar;  
  set(h,'Position',[0.92 0.33 0.01 .26])  
  caxis([-2 1]);
  
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
  land_plot(sum(hFacC,3)>0) = NaN;
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
    h = title(['Potential temperature (^oC) on ',datestr(thedate)],'FontSize',18);
%   h = title(['\zeta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
%   h = title(['\delta/|f| at t= ',num2str(round(tt(n)),'%3d'),' days']);
  %  h = title(['\zeta/f at t= ',num2str(tt(n)/365,'%.1f'),' years']);  
  set(h,'Position',get(h,'Position')+[0 0.01 0])



  M(n) = getframe(gcf);  
   
end