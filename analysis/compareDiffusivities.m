%%%
%%% compareDiffusivities.m
%%%
%%% Compares results of submesoscale-resolving simulations with different
%%% biharmonic diffusivities.
%%%

%%% To set paths
setExpname;

%%% Select output snapshot to use for diagnostics
n = 250;

%%% List of experiments to make plots for
expnames = {...
  ... % 'hires_nest_oneseventysecond_notides_RTOPO2_old', ...
  'hires_nest_oneseventysecond_notides_RTOPO2_visctest', ...
  'hires_nest_oneseventysecond_notides_RTOPO2_visctest1', ...
  'hires_nest_oneseventysecond_notides_RTOPO2_visctest2', ...
  'hires_nest_oneseventysecond_notides_RTOPO2_visctest3', ...
  'hires_nest_oneseventysecond_notides_RTOPO2_visctest4'};
K4 = [ ... %0 ...
  7e4 2.3e5 7e5 2.3e6 7e6];
lastiter = [ ... %522720 
  542160 544320 557280 548640 540000]*20/86400*2;

%%% Main loop over experiments
for m=1:length(expnames)
 
  %%% Read experiment data
  expname = expnames{m};
  loadexp;

  %%% Physical parameters
  Omega = 2*pi*366/365/86400;
  ff = 2*Omega*sind(YG);
  g = 9.81;

  %%% Diagnostic indix corresponding to instantaneous velocity
  diagnum = length(diag_frequency);
  
  %%% This needs to be set to ensure we are using the correct output
  %%% frequency
  diagfreq = diag_frequency(diagnum);


  
  %%% Frequency of diagnostic output
  dumpFreq = abs(diagfreq);
  nDumps = round(nTimeSteps*deltaT/dumpFreq);
  dumpIters = round(nIter0 + (1:nDumps)*dumpFreq/deltaT);
  dumpIters = dumpIters(dumpIters >= nIter0);
  nDumps = length(dumpIters);


  tt =  dumpIters(n)*deltaT/86400;
  
  %%% Read velocity snapshot
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n));      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n));      
  wvel = rdmdsWrapper(fullfile(exppath,'/results/WVEL_inst'),dumpIters(n));      
  salt = rdmdsWrapper(fullfile(exppath,'/results/SALT_inst'),dumpIters(n));    
  
  uvel(hFacW==0) = NaN;
  vvel(hFacS==0) = NaN;
  
  %%% Vorticity on a z-level
  vort = zeros(Nx,Ny);
  zlev = 1;
  vort(:,2:Ny) = - (uvel(:,2:Ny,zlev)-uvel(:,1:Ny-1,zlev))./DYC(:,2:Ny);
  vort(2:Nx,:) = vort(2:Nx,:) + (vvel(2:Nx,:,zlev)-vvel(1:Nx-1,:,zlev))./DXC(2:Nx,:);
  
  %%% Divergence on a z-level
  divg = zeros(Nx,Ny);
  zlev = 1; 
  divg(:,1:Ny-1) = (vvel(:,2:Ny,zlev)-vvel(:,1:Ny-1,zlev))./DYG(:,1:Ny-1);
  divg(1:Nx-1,:) = divg(2:Nx,:) + (uvel(2:Nx,:,zlev)-uvel(1:Nx-1,:,zlev))./DXG(1:Nx-1,:);

  
  %%% Make vorticity plot
  thedate = datenum('2008-01-00')+floor(tt);
  titlestr = ['\zeta/|f| on ',datestr(thedate),', K4 = ',num2str(K4(m),'%.1e'),' m^4/s'];
  outfname = ['vort_K4=',num2str(K4(m),'%.1e'),'_iter=',num2str(n),'.png'];
  map_plot(vort./abs(ff),titlestr,outfname,bathy,SHELFICEtopo,XC,YC,XG,YG);

  %%% Make divergence plot
  titlestr = ['\delta/|f| on ',datestr(thedate),', K4 = ',num2str(K4(m),'%.1e'),' m^4/s'];
  outfname = ['divg_K4=',num2str(K4(m),'%.1e'),'_iter=',num2str(n),'.png'];
  map_plot(divg./abs(ff),titlestr,outfname,bathy,SHELFICEtopo,XC,YC,XG,YG);
  

  
  %%% Extract ROI for spectra
  xidx = 793:1081;
  yidx = 671:970;
  uu = uvel(xidx,yidx,1);
  vv = vvel(xidx,yidx,1);
  ww = wvel(xidx,yidx,10);  
  ss = salt(xidx,yidx,1);
  qq = vort(xidx,yidx);

  xx = 1:length(xidx);
  
  %%% Detrend
  for j=1:size(uu,2)
    uu(:,j) = uu(:,j) - (uu(end,j)-uu(1,j)).*(xx'-xx(1))/(xx(end)-xx(1)+1);
  end
  for j=1:size(uu,2)
    vv(:,j) = vv(:,j) - (vv(end,j)-vv(1,j)).*(xx'-xx(1))/(xx(end)-xx(1)+1);
  end
  for j=1:size(ww,2)
    ww(:,j) = ww(:,j) - (ww(end,j)-ww(1,j)).*(xx'-xx(1))/(xx(end)-xx(1)+1);
  end
  for j=1:size(ss,2)
    ss(:,j) = ss(:,j) - (ss(end,j)-ss(1,j)).*(xx'-xx(1))/(xx(end)-xx(1)+1);
  end
  for j=1:size(qq,2)
    qq(:,j) = qq(:,j) - (qq(end,j)-qq(1,j)).*(xx'-xx(1))/(xx(end)-xx(1)+1);
  end

  %%% Transform to Fourier space
  uu_fft = fft(uu);
  vv_fft = fft(vv);
  ww_fft = fft(ww);
  ss_fft = fft(ss);
  qq_fft = fft(qq);
  ww_fft(1,:,:) = 0;
  ss_fft(1,:,:) = 0;

  %%% Plot horizotal energy spectrum
  figure(100);
  if (m == 1)
    clf;
  end
  % loglog(abs(mean(uu_fft(1:end/2,:),2)).^2+abs(mean(vv_fft(1:end/2,:),2)).^2);
  loglog(mean(abs(uu_fft(1:end/2,:)).^2+abs(vv_fft(1:end/2,:)).^2,2));
  if (m == 1)
    hold on;
    legstr1 = {};
  end
  legstr1 = {legstr1{:},['K4 = ',num2str(K4(m),'%.1e')]};
  if (m == length(expnames))
    hold off;
    set(gca,'FontSize',18);
    xlabel('Wavenumber');
    ylabel('E(k)')
    legend(legstr1,'Location','SouthWest');
    hold on;
    loglog([1e1 1e2],(1e2/1e1*[1e1 1e2]).^-3,'k--');
    text(1e1,1e2,'k^-^3');
    hold off;
    title(['Horizontal velocity, zonal power spectrum, ',thedate]);
    print('-dpng','-r300',fullfile('Figures','diffusivityComparison',['HorizontalEnergySpectrum_n=',num2str(n),'.png']));
  end

  %%% Plot vertical energy spectrum 
  figure(101);
  if (m == 1)
    clf;
  end
  % loglog(abs(mean(ww_fft(1:end/2,:),2)).^2);
  loglog(mean(abs(ww_fft(1:end/2,:)).^2,2));
  if (m == 1)
    hold on;
    legstr2 = {};
  end
  legstr2 = {legstr2{:},['K4 = ',num2str(K4(m),'%.1e')]};
  if (m == length(expnames))
    hold off;
    set(gca,'FontSize',18);
    xlabel('Wavenumber');
    ylabel('E(k)')
    legend(legstr2,'Location','SouthWest');
    hold on;
    loglog([1e1 1e2],2e-3*(1/1e1*[1e1 1e2]).^-1,'k--');
    text(1e1,2e-3,'k^-^1');
    hold off;
    title(['Vertical velocity, zonal power spectrum, ',thedate]);
    print('-dpng','-r300',fullfile('Figures','diffusivityComparison',['VerticalEnergySpectrum_n=',num2str(n),'.png']));
  end

  %%% Plot salinity variance spectrum 
  figure(102);
  if (m == 1)
    clf;
  end
  loglog(mean(abs(ss_fft(1:end/2,:)).^2,2));
  if (m == 1)
    hold on;
    legstr3 = {};
  end
  legstr3 = {legstr3{:},['K4 = ',num2str(K4(m),'%.1e')]};
  if (m == length(expnames))
    hold off;
    set(gca,'FontSize',18);
    xlabel('Wavenumber');
    ylabel('S(k)')
    legend(legstr3,'Location','SouthWest');
    hold on;
    loglog([1e1 1e2],1e0*(1/1e1*[1e1 1e2]).^-2.5,'k--');
    text(1e1,1e0,'k^-^2^.^5');
    hold off;
    title(['Salinity, zonal power spectrum, ',thedate]);
    print('-dpng','-r300',fullfile('Figures','diffusivityComparison',['SalinityVarianceSpectrum_n=',num2str(n),'.png']));
  end

  %%% Plot enstrophy power spectrum 
  figure(103);
  if (m == 1)
    clf;
  end
  loglog(mean(abs(qq_fft(1:end/2,:)).^2,2));
  if (m == 1)
    hold on;
    legstr4 = {};
  end
  legstr4 = {legstr4{:},['K4 = ',num2str(K4(m),'%.1e')]};
  if (m == length(expnames))
    hold off;
    set(gca,'FontSize',18);
    xlabel('Wavenumber');
    ylabel('Z(k)')
    legend(legstr4,'Location','SouthWest');
    hold on;
    loglog([1e1 1e2],1e-5*(1/1e1*[1e1 1e2]).^-1,'k--');
    text(1e1,1e-5,'k^-^1');
    hold off;
    title(['Enstrophy, zonal power spectrum, ',thedate]);
    print('-dpng','-r300',fullfile('Figures','diffusivityComparison',['EnstrophyPowerSpectrum_n=',num2str(n),'.png']));
  end

end








%%%
%%% Convenience function to generate a plot of vorticity or divergence
%%% fields.
%%%
function map_plot (vort,titlename,outfname,bathy,SHELFICEtopo,XC,YC,XG,YG)

  lonMin = -36;
  lonMax = -32;
  latMin = -77;
  latMax = -76;

  Nx = size(vort,1);
  Ny = size(vort,2);

  figure(1);
  set(gcf,'Color','w');
  icecolor = [186 242 239]/255;

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
  'ParallelLabel', 'on');  


  axis off;
  setm(gca,'MLabelParallel',-20)


  pcolorm(YG,XG,vort);
  shading flat
  colormap(cmocean('balance'));

  

  hold on;
  
  thick_plot = -bathy;
  
  [cs,C]=contourm(YC,XC,thick_plot,[500 1000 2000 3000 4000],'EdgeColor','k');
  chandle = clabelm(cs,C);
  set(chandle,'fontsize',14,'Color','k','BackgroundColor','none','Edgecolor','none') 
  hold off;

  %%% Add colorbar and title
  h = colorbar;  
  set(h,'Position',[0.92 0.33 0.01 .26])  
  caxis([-1 1]);
  
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
  'ParallelLabel', 'off');
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
  'ParallelLabel', 'off');  
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

  %%% Add title
  h = title(titlename,'FontSize',18);
  % set(h,'Position',get(h,'Position')+[0 0.01 0])

  print('-dpng','-r300',fullfile('Figures','diffusivityComparison',outfname));

end