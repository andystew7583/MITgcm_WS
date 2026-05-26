%%%
%%% plotCurtainComparison.m
%%%
%%% Compares simulations with and without a "curtain" across the Filchner
%%% Trough.
%%%

%%% Pointer to experiment directory
expdir = '../experiments';

%%% Pointer to storage directory for output .mat files
proddir = './products_WCbatch';

%%% Load pre-computed diagnostics
load(fullfile(proddir,"curtain_diags.mat"));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT MELT TIME SERIES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
colororder = get(gca,'ColorOrder');
for m = 1:length(expnames)
  expname = expnames{m}; 
  plot (diags_curtain{m}.tt,smooth(diags_curtain{m}.FRISmelt,24),'LineWidth',2,'Color',colororder(m,:));
  if (m == 1)
    hold on;
  end
end
for m = 1:length(expnames)
  expname = expnames{m};  
  plot (diags_curtain{m}.tt,diags_curtain{m}.FRISmelt,'LineWidth',0.5,'Color',colororder(m,:));  
end
hold off;
xlabel('Time (years)')
ylabel('$F_{\mathrm{melt}}$ (Gt/yr)','interpreter','latex','FontSize',16);
set(gca,'FontSize',16);
axis([0 36 0 600]);
legend(titles,'Location','Best');

%%% Print to file
print('-dpng','-r150',fullfile('figures','curtain',['melt_timeseries.png']));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT FILCHNER TROUGH TRANSECTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define grid to extract data sections 
stat_lat = [-83.35, -81.83, -81.38, -81.06, -80.38, -78.77, -77.24, -75.45, -74.28 -73.00];
stat_lon = [-60.25, -57.17, -51.36, -43.83, -41.43, -39.07, -37.33, -33.22, -32.20 -31.80];
Nsec = 401;
dLon = (stat_lon(end)-stat_lon(1))/(Nsec-1);
secLons = stat_lon(1):dLon:stat_lon(end);
secLats = interp1(stat_lon,stat_lat,secLons,'spline');
startLat = secLats(1);
endLat = secLats(end);






for m = 1:length(expnames)
  
  expname = expnames{m};
  Nx = size(diags_curtain{m}.theta_tavg,1);
  Ny = size(diags_curtain{m}.theta_tavg,2);
  Nr = size(diags_curtain{m}.theta_tavg,3);

  %%% Pre-compute indices for interpolation
  xidx = zeros(Nsec,1);
  yidx = zeros(Nsec,1);
  XC = diags_curtain{m}.XC;
  YC = diags_curtain{m}.YC;
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
  
  end
  
  %%% Remove dry grid cells
  T =  diags_curtain{m}.theta_tavg;
  S =  diags_curtain{m}.salt_tavg;
  T(diags_curtain{m}.hFacC==0) = NaN;
  S(diags_curtain{m}.hFacC==0) = NaN;

  %%% Extract data along defined sections   
  secT = zeros(Nsec,Nr); 
  secS = zeros(Nsec,Nr); 
  secB = zeros(Nsec,1);
  secI = zeros(Nsec,1);
  secH = zeros(Nsec,Nr);
  for n=1:Nsec
    secT(n,:) = squeeze(T(xidx(n),yidx(n),:));
    secS(n,:) = squeeze(S(xidx(n),yidx(n),:));
    secH(n,:) = squeeze(diags_curtain{m}.hFacC(xidx(n),yidx(n),:));
    secB(n) = diags_curtain{m}.bathy(xidx(n),yidx(n));
    secI(n) = diags_curtain{m}.SHELFICEtopo(xidx(n),yidx(n));
  end

  %%% Meshgrid for plotting
  [ZZ,LA] = meshgrid(diags_curtain{m}.RC,secLats);
  for n=1:Nsec
    kmin = find(secH(n,:)>0,1,'first');
    kmax = find(secH(n,:)>0,1,'last');
    if (isempty(kmin) || isempty(kmax))
      continue;
    end
    if (kmin == kmax)    
      continue;
    end
    ZZ(n,kmin) = diags_curtain{m}.RF(kmin+1)+0.5*secH(n,kmin)*diags_curtain{m}.DRF(kmin);
    ZZ(n,kmax) = diags_curtain{m}.RF(kmax)-0.5*secH(n,kmax)*diags_curtain{m}.DRF(kmax);
    
  end





  %%% Plotting options
  scrsz = get(0,'ScreenSize');
  fontsize = 14;
  plotloc = [0.08  0.075  0.84  0.9];
  framepos = [100   462   810   520];
  cb_pos = [0.95 0.07 0.01 0.9];    
  icecolor = [186 242 239]/255;
  cmap = cmocean('thermal',100);
  
  %%% Set up the figure
  handle = figure(m);
  set(handle,'Position',framepos);
  
  clf;
  pcolor(LA,-ZZ,secT);
  set(gca,'Position',plotloc);
  set(gca,'YDir','reverse');
  shading interp;
  hold on
  % plot(secLats,-secB,'k-','LineWidth',3);
  % Iidx = find(secI==0,1,'first');
  % plot(secLats(1:Iidx),-secI(1:Iidx),'k-','LineWidth',3);
  [C,h] = contour(LA,-ZZ,secS,'EdgeColor',[1 1 1]);
  clabel(C,h,'Color','w','FontSize',14);
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
  caxis([-2.5 0]);

  ax1 = gca;
  ax2 = axes('Position',get(ax1,'Position'));
  msk = ones(Nsec,Nr);
  msk(ZZ<repmat(secI,[1 Nr])) = NaN;
  pcolor(LA,-ZZ,msk);
  shading flat;
  colormap(ax2,icecolor);
  set(ax2,'XTick',[]);
  set(ax2,'YTick',[]);
  set(ax2,'YLim',[0 2000]);
  set(ax2,'YDir','reverse');
  set(ax2,'XLim',[startLat endLat]);
  set(ax2,'Color','None')
 
  %%% Add title 
  titlestr = titles{m};
  annotation('textbox',[0.1 0.9 0.5 0.05],'String',[titlestr],'interpreter','latex','FontSize',fontsize+4,'LineStyle','None');   
  
  %%% Print to file
  print('-dpng','-r150',fullfile('figures','curtain',[expname,'_FTsection.png']));


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% MAP WATER COLUMN THICKNESS %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  thic_plot =diags_curtain{m}.SHELFICEtopo-diags_curtain{m}.bathy;
  thic_plot(thic_plot == 0) = NaN;
  handle = figure(10+m);
  set(handle,'Position',framepos);
  pcolor(diags_curtain{m}.XC,diags_curtain{m}.YC,thic_plot);
  set(gca,'Position',plotloc);
  shading interp;
  colormap(flip(haxby(30)));
  hold on;
  plot(secLons,secLats,'k--');
  hold off;
  axis([-80 -20 -83 -72])
  caxis([0 3000]);
  set(gca,'FontSize',fontsize);  
  cbhandle = colorbar;
  set(cbhandle,'Position',cb_pos);
  xlabel('Longitude');
  ylabel('Latitude');
  titlestr = titles{m};
  annotation('textbox',[0.1 0.9 0.5 0.05],'String',[titlestr],'interpreter','latex','FontSize',fontsize+4,'LineStyle','None');   

  %%% Print to file
  print('-dpng','-r150',fullfile('figures','curtain',[expname,'_columnthickness.png']));

end
 


