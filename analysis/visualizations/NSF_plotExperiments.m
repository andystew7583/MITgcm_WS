
%%%
%%% NSF_plotExperiments.m
%%%
%%% Plots hydrography in our proposed perturbation experiments, for our NSF
%%% proposal in January 2025.
%%%

%%% Options
expdir = '../experiments';
expnames = {'WC_onethird_ref','WC_onethird_dpyc-150_strat1e-4','WC_onethird_dpyc-150_strat3e-4'};
% iters = [1774980,591660,2366640]; %%% Dec 2015
iters = [1736632,553312,2328292]; %%% May 2015

expnames = {'WC_onethird_ref','WC_onethird_dpyc-150_strat1e-4','WC_onethird_dpyc-150_strat1e-4_aU0.5_aV0.5_ws'};
iters = [1736632,553312,1144972]; %%% May 2015


expnames = {'WC_onethird_ref','WC_onethird_dpyc-150_strat1e-4','WC_onethird_dpyc-100_strat2e-4'};
iters = [1736632,553312,1736632]; %%% May 2015

%%% Define grid to extract data sections
startLat = -82.3;
endLat = -72;
startLon = -62;
endLon = -40;
Nsec = 201;
dLat = (endLat-startLat)/(Nsec-1);
dLon = (endLon-startLon)/(Nsec-1);
secLats = startLat:dLat:endLat;
secLons = startLon:dLon:endLon;

%%% High-level plot options
axlabels = {'\textbf{A}','\textbf{B}','\textbf{C}','\textbf{D}','\textbf{E}','\textbf{F}'};
axpos = zeros(6,4);
axpos(1,:) = [0.05 0.67 0.45 0.28];
axpos(2,:) = [0.52 0.7 0.4 0.25];
axpos(3,:) = [0.05 0.345 0.45 0.28];
axpos(4,:) = [0.52 0.375 0.4 0.25];
axpos(5,:) = [0.05 0.03 0.45 0.28];
axpos(6,:) = [0.52 0.05 0.4 0.25];
cb_pos = [0.95 0.05 0.015 0.9];

%%% Set up the figure
figure(201)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34  1000  926]);

%%% Loop through and make the plots
for m = 1:3   
  ax = subplot('Position',axpos(2*m-1,:)+[.1 0 0 0]);
  plotTempMap(expdir,expnames{m},iters(m),axpos(2*m-1,:),secLats,secLons);
  ax = subplot('Position',axpos(2*m,:));
  plotTempSlice(expdir,expnames{m},iters(m),axpos(2*m,:),secLats,secLons);
end

cbhandle = colorbar(ax);
set(cbhandle,'Position',cb_pos);
title(cbhandle,'$^\circ$C','Fontsize',14,'interpreter','latex');

%%% Add panel labels
annotation('textbox',[axpos(1,1)+0.05 axpos(1,2)-0.0 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.05 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)+0.05 axpos(3,2)-0.0 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.05 axpos(4,2)-0.05 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(5,1)+0.05 axpos(5,2)-0.0 0.03 0.03],'String',axlabels{5},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(6,1)-0.05 axpos(6,2)-0.05 0.03 0.03],'String',axlabels{6},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


figure1 = gcf;

% Create textbox
annotation(figure1,'textbox',...
  [0.013 0.455097192224623 0.03 0.0300000000000002],'String','\textbf{(ii)}',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',20,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.00800000000000002 0.802829373650109 0.03 0.0299999999999999],...
  'String','\textbf{(i)}',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',20,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',[0.012 0.159200863930887 0.03 0.03],...
  'String','\textbf{(iii)}',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',20,...
  'FitBoxToText','off');





%%%%%%%%%%%%%%%%%%%%%%%
%%% TEMPERATURE MAP %%%
%%%%%%%%%%%%%%%%%%%%%%%

function plotTempMap (expdir,expname,iter,axpos,secLats,secLons)

  %%% Load experiment parameters
  loadexp;
  
  %%% Plotting range 
  latMin_b = min(min(YC));
  latMax_b = -70;
  lonMin_b = -78;
  lonMax_b = -25; 
  
  %%% Plotting options
  fontsize = 13;
  clim = [-2.8 -0.8];
  cmap = cmocean('balance',32);
  cmap = cmap(8:23,:);  
  bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
    
  %%% Load temperature snapshot
  T = rdmdsWrapper(fullfile(exppath,'results','THETA'),iter);
  
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
      % FF(i,j) = wp(i,j)*T(i,j,kp(i,j)) + wn(i,j)*T(i,j,kn(i,j));
      FF(i,j) = T(i,j,kmax(i,j));
    end
  end
  FF(sum(hFacC,3)==0) = NaN;
  
  %%% Free up memory
  clear('T');
  
  
  
  
  
  
  
  
  
  
  

  
  %%% Set up map plot
  axesm('eqaconicstd',...
    'fontsize',fontsize,...
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
  
  %%% Plot temperature in color contours
  pcolorm(YC,XC,FF);
  shading interp
  caxis(clim);
  colormap(gca,cmap);
 
  %%% Plot transect location for VHF figure
  hold on;
  plotm(secLats,secLons,'k-.','LineWidth',1);
  hold off;
  
  %%% Add bathymetry contours
  hold on;
  [cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
  hh = clabelm(cs,C);
  set(hh,'fontsize',8,'Color',[.05 .05 .05],'BackgroundColor','none','Edgecolor','none')       
  % Nlon = 101;
  % dLon = (lonMax_c-lonMin_c)/(Nlon-1);
  % plotm([latMin_c*ones(1,Nlon) latMax_c*ones(1,Nlon) latMin_c],[lonMin_c:dLon:lonMax_c lonMax_c:-dLon:lonMin_c lonMin_c],'r--','LineWidth',2);
  hold off;  
  tightmap;

  %%% Add axis labels
  xlabel('Longitude','interpreter','latex');
  ylabel('Latitude','interpreter','latex');
  set(gca,'Position',axpos);
  hold off


end











%%%%%%%%%%%%%%%%%%
%%% TEMP SLICE %%%
%%%%%%%%%%%%%%%%%%
function plotTempSlice (expdir,expname,iter,axpos,secLats,secLons)

  %%% Load experiment parameters
  loadexp;

  %%% Plotting options                
  fontsize = 14;
  labelspacing = 200;  
  icecolor = [186 242 239]/255;
  clim = [-2.8 -0.8];
  cmap = cmocean('balance',32);
  cmap = cmap(8:23,:);  

  %%% Extract data along defined sections
  S = rdmdsWrapper(fullfile(exppath,'results','SALT'),iter);
  T = rdmdsWrapper(fullfile(exppath,'results','THETA'),iter);
  S(hFacC==0) = NaN;
  T(hFacC==0) = NaN;
  Nsec = length(secLats);
  secS = zeros(Nsec,Nr);
  secT = zeros(Nsec,Nr);
  secB = zeros(Nsec,1);
  secI = zeros(Nsec,1);
  secH = zeros(Nsec,Nr);
  for n=1:Nsec
    jm = find(secLats(n)>=YC(1,:),1,'last');
    im = find(secLons(n)>=XC(:,1),1,'last');
    jp = jm+1;
    ip = im+1;
    wp_y = (secLats(n)-YC(1,jm))/(YC(1,jp)-YC(1,jm));
    wm_y = 1-wp_y;
    wp_x = (secLons(n)-XC(im,1))/(XC(ip,1)-XC(im,1));
    wm_x = 1-wp_x;
    if (wp_x > wm_x)
      xidx = ip;
    else
      xidx = im;
    end
    if (wp_y > wm_y)
      yidx = jp;
    else
      yidx = jm;
    end
    secS(n,:) = squeeze(S(xidx,yidx,:));
    secT(n,:) = squeeze(T(xidx,yidx,:));  
    secH(n,:) = squeeze(hFacC(xidx,yidx,:));
    secB(n) = bathy(xidx,yidx);
    secI(n) = SHELFICEtopo(xidx,yidx);
  end
  
  
  
  
  
  
  
  
  
  
  
                  

  
  
  
  
  [ZZ,LA] = meshgrid(RC,secLats);
  for n=1:Nsec
    kmin = find(secH(n,:)>0,1,'first');
    kmax = find(secH(n,:)>0,1,'last');
    if (isempty(kmin) || isempty(kmax))
      continue;
    end
    if (kmin == kmax)
      secT(n,kmin) = NaN;
      secS(n,kmin) = NaN;
      continue;
    end
    ZZ(n,kmin) = RF(kmin);% - (1-sechFac(n,kmin))*DRF(kmin);
    ZZ(n,kmax) = RF(kmax+1);% - (1-sechFac(n,kmax))*DRF(kmax);
    
  end
  
  
  % axes('Position',axpos);
  pcolor(LA,-ZZ,secT);
  set(gca,'YDir','reverse');
  shading interp;
  hold on
  [C,h] = contour(LA,-ZZ,secS,[33:.2:34 34.1:.1:35],'EdgeColor','k');
  clabel(C,h,'FontSize',12);
  plot(secLats,-secB,'k-','LineWidth',3);
  Iidx = find(secI==0,1,'first');
  plot(secLats(1:Iidx),-secI(1:Iidx),'k-','LineWidth',3);
  hold off;
  set(gca,'YLim',[0 1400]);
  set(gca,'XLim',[-82.3 (-73)]);
  caxis(clim);
  colormap(gca,cmap);  
  set(gca,'FontSize',fontsize);
  xlabel('Latitude')
  ylabel('Depth (m)');
  set(gca,'Color',[.8 .8 .8]);
  set(gca,'Position',axpos);
  % title(['Instantaneous potential temperature, salinity, ',datestr(datenum('01-Jan-2008')+iter*deltaT/86400)]);
  
  ax1 = gca;
  ax2 = axes('Position',get(ax1,'Position'));
  msk = ones(Nsec,Nr);
  msk(ZZ<repmat(secI,[1 Nr])) = NaN;
  pcolor(LA,-ZZ,msk);
  shading interp;
  colormap(ax2,icecolor);
  set(ax2,'XTick',[]);
  set(ax2,'YTick',[]);
  set(ax2,'YLim',[0 1400]);
  set(ax2,'YDir','reverse');
  set(ax2,'XLim',[-82.3 (-73)]);
  set(ax2,'Color','None')


end
