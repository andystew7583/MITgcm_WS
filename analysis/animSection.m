%%%
%%% animSection.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of properties along a user-defined section. Requires
%%% that the experiment be generated using the 'newexp' package. 
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = false;

%%% Read experiment data
loadexp;

% %%% Define grid to extract data sections - front of RL and SWIT
% startLat = -74.2;
% endLat = -72.2;
% startLon = -26;
% endLon = -17;
% Nsec = 401;
% dLat = (endLat-startLat)/(Nsec-1);
% dLon = (endLon-startLon)/(Nsec-1);
% secLats = startLat:dLat:endLat;
% secLons = startLon:dLon:endLon;

%%% Define grid to extract data sections - front of RL and SWIT
stat_lat = [-83.35, -81.83, -81.38, -81.06, -80.38, -78.77, -77.24, -75.45, -74.28 -73.00];
stat_lon = [-60.25, -57.17, -51.36, -43.83, -41.43, -39.07, -37.33, -33.22, -32.20 -31.80];
Nsec = 401;
dLon = (stat_lon(end)-stat_lon(1))/(Nsec-1);
secLons = stat_lon(1):dLon:stat_lon(end);
secLats = interp1(stat_lon,stat_lat,secLons,'spline');
startLat = secLats(1);
endLat = secLats(end);

%%% Select diagnostic variable to animate

diagnum = 5;
outfname =diag_fileNames{1,diagnum};
% outfname = 'THETA_12hourly';
% outfname = 'Eta';

%%% Data index in the output data files
outfidx = 1;

%%% Specify color range
set_crange = 1;

%%% Year to start animnation
startyr = 4;

switch (diagnum)
  
  case 4   

    crange = [-2.5 1]; %/%% Filchner temp
    cmap = cmocean('thermal',20);

  case 5

    crange = [32 34.8]; %%% salinity    
    cmap = cmocean('haline',30);

end
% crange = [-3 3]; %%%temp
% crange = [33.4 34.65]; %%% salinity
% crange = [0 10]; %%%% for KPP hbl
% crange = [0 1]; %%% For sea ice area
% crange = [-.6 .6]; %%% For velocities or stresses
% crange = [-1 1]*1e-4; %%% For freshwater fluxes
% crange =[-100 100]; %%% Qnet
% crange = [-300 300]; %%% swnet
% crange =[-.1 .1]; %%% SFLUX
% crange =[-1 1]*3e-3; %%% WVEL
% % crange = [0 3];
% crange = [-0.5 0.5];
% crange = [0.3 1]; %% SSH
% crange = [-0.3 0.3]; %%% Melt rate 

% cmap = pmkmp(100,'Swtth');
% cmap = cmocean('ice',100);
% cmap = haxby(56);
% cmap = jet(200);
% cmap = redblue(20);
% cmap = cmocean('amp',100);
% %
% titlestr = 'Salinity (g/kg)';
% titlestr = 'Surface salinity (g/kg)';
% titlestr = 'Temperature ($^\circ$C)';
% titlestr = 'Zonal velocity (m/s)';
% titlestr = 'Surface temperature ($^\circ$C)';
% titlestr = 'Sea ice concentration';
titlestr = '';

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
years = 2007:1:2015;


%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics
% nIter0 = 1278720;
% deltaT = 200
% nIter0 = 587520;
% nIter0 = 394509; 

%%% Default
dumpFreq = abs(diag_frequency(diagnum));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% For daily/12-hourly outputs
% dumpStart = 1578240;
% % dumpStart = 1945440;
% dumpStep = 86400/2/60;
% nDumps = 731;
% dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;







 %%% Pre-compute indices for interpolation
secA = zeros(Nsec,Nr); %%% To store sections of property of interest 
secB = zeros(Nsec,1);
secI = zeros(Nsec,1);
secH = zeros(Nsec,Nr);
xidx = zeros(Nsec,1);
yidx = zeros(Nsec,1);
for n=1:Nsec
  jm = find(secLats(n)>=YC(1,:),1,'last');
  im = find(secLons(n)>=XC(:,1),1,'last');
  jp = jm+1;
  ip = im+1;
  wp_y = (secLats(n)-YC(1,jm))/(YC(1,jp)-YC(1,jm));
  wm_y = 1-wp_y;
  wp_x = (secLons(n)-XC(im,1))/(XC(ip,1)-XC(im,1));
  wm_x = 1-wp_x;
  if (wp_x > wm_x) %%% Nearest-neighbor interpolation
    xidx(n) = ip;
  else
    xidx(n) = im;
  end
  if (wp_y > wm_y)
    yidx(n) = jp;
  else
    yidx(n) = jm;
  end
  secH(n,:) = squeeze(hFacC(xidx(n),yidx(n),:));
  secB(n) = bathy(xidx(n),yidx(n));
  secI(n) = SHELFICEtopo(xidx(n),yidx(n));
end

%%% Meshgrid for plotting
[ZZ,LA] = meshgrid(RC,secLats);
for n=1:Nsec
  kmin = find(secH(n,:)>0,1,'first');
  kmax = find(secH(n,:)>0,1,'last');
  if (isempty(kmin) || isempty(kmax))
    continue;
  end
  if (kmin == kmax)
    secA(n,kmin) = NaN;
    continue;
  end
  % ZZ(n,kmin) = RF(kmin);% - (1-sechFac(n,kmin))*DRF(kmin);
  ZZ(n,kmin) = RF(kmin+1)+0.5*secH(n,kmin)*DRF(kmin);
  % ZZ(n,kmax) = RF(kmax+1);% - (1-sechFac(n,kmax))*DRF(kmax);
  ZZ(n,kmax) = RF(kmax)-0.5*secH(n,kmax)*DRF(kmax);
  
end









%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 14;
plotloc = [0.08  0.07  0.84  0.9];
framepos = [100   462   810   520];
cb_pos = [0.95 0.07 0.01 0.9];    
icecolor = [186 242 239]/255;

%%% Set up the figure
handle = figure(22);
set(handle,'Position',framepos);
M = moviein(length(dumpIters));

for n=startyr*12:length(dumpIters) 

  dumpIters(n);
    
  t = dumpIters(n)*deltaT/t1year;
  n
  
  tyears(n) = t;
  tdays(n) = dumpIters(n)*deltaT/t1day;
  
  A = rdmdsWrapper(fullfile(exppath,'results',outfname),dumpIters(n));
  if (isempty(A))
%     continue;
    error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
  end     
  
  if (~isempty(find(isnan(A))))
    error(['Found NaNs at iter=',num2str(dumpIters(n))]);
    break
  end
  
    



  %%% Extract data along defined sections 
  A(hFacC==0) = NaN;
  for m=1:Nsec
    secA(m,:) = squeeze(A(xidx(m),yidx(m),:));        
  end
  
  
  
  clf;
  pcolor(LA,-ZZ,secA);
  set(gca,'Position',plotloc);
  set(gca,'YDir','reverse');
  shading interp;
  hold on
  plot(secLats,-secB,'k-','LineWidth',3);
  Iidx = find(secI==0,1,'first');
  plot(secLats(1:Iidx),-secI(1:Iidx),'k-','LineWidth',3);
  hold off;
  set(gca,'YLim',[0 2000]);
  set(gca,'XLim',[startLat endLat]);  
  colormap(gca,cmap);
  cbhandle = colorbar;
  set(cbhandle,'Position',cb_pos);
  set(gca,'FontSize',fontsize);
  xlabel('Latitude')
  ylabel('Depth (m)');
  set(gca,'Color',[.8 .8 .8]);
  % title(['Instantaneous potential temperature, salinity, ',datestr(datenum('01-Jan-2008')+iter*deltaT/86400)]);
  if (set_crange)  
    caxis(crange);
  end  

  ax1 = gca;
  ax2 = axes('Position',get(ax1,'Position'));
  msk = ones(Nsec,Nr);
  msk(ZZ<repmat(secI,[1 Nr])) = NaN;
  pcolor(LA,-ZZ,msk);
  shading interp;
  colormap(ax2,icecolor);
  set(ax2,'XTick',[]);
  set(ax2,'YTick',[]);
  set(ax2,'YLim',[0 2000]);
  set(ax2,'YDir','reverse');
  set(ax2,'XLim',[startLat endLat]);
  set(ax2,'Color','None')

 
  %%% Finish the plo
  if (~isempty(titlestr))
%     title([titlestr,', $t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');
    thedate = datenum('2008-01-01')+floor(tdays(n));
%     annotation('textbox',[0.3 0.95 0.5 0.05],'String',[titlestr,', ',months{mod(n-1,12)+1},' ',num2str(years(floor((n-1)/12)+1))],'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');   
    annotation('textbox',[0.3 0.95 0.5 0.05],'String',[titlestr,', ',datestr(thedate)],'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');   
  end

  
  M(n) = getframe(gcf);
  
end





