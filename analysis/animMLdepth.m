%%%
%%% animMLdepth.m
%%%
%%% Animates maps of mixed layer depth in MITgcm output.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = false;

%%% Read experiment data
% loadexp;

%%% Vertical level to which potential density should be referenced
k_pot_dens = 1;

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 1;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 1;

%%% Set true to plot the field in the topmost wet cell at each horizontal
%%% location
topplot = 0;

%%% Set true to plot the field in the middle of the water column at each
%%% horizotal location
midplot = 0;

%%% Set true for a zonal average
yzavg = 0;

%%% Layer to plot in the y/z plane
%%%for 1/3 DEGREE yzlayer = 126; 1/6 = 404;
yzlayer = 322;

% load ../newexp/ELEV.mat
% 
% 
% Mice_elev=Mice_elev';


%%% Specify color range
set_crange = 1;
crange = [0 600];
cmap = cmocean('amp',30);

% titlestr = 'Salinity (g/kg)';
% titlestr = 'Surface salinity (g/kg)';
% titlestr = 'Temperature ($^\circ$C)';
% titlestr = 'Surface temperature ($^\circ$C)';
titlestr = 'Mixed layer depth (m)';
% titlestr = '';

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
years = 2007:1:2015;


%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics
% nIter0 = 1278720;
% deltaT = 200
% nIter0 = 587520;
% nIter0 = 394509; 

%%% For 1/24 run with tides
% dumpFreq = abs(diag_frequency(diagnum));
% nDumps = round(endTime/dumpFreq);
% dumpIters = round((1:nDumps)*dumpFreq/deltaT);
% dumpIters = dumpIters(dumpIters > nIter0);

%%% For daily/12-hourly outputs
dumpStart = 1578240;
% dumpStart = 1945440;
dumpStep = 86400/2/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;




%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 14;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
%   plotloc = [0.0855    0.0888    0.7916    0.8624];  
  plotloc = [0.04  0.02   0.84  0.9];
%   framepos = [100   306   936   676];
%   framepos = [100   462   810   520];
  framepos = [441         253        1529         994];
end

%%% Set up the figure
handle = figure(22);
set(handle,'Position',framepos);
M = moviein(length(dumpIters));

Amean = [];
Amax = [];

% for n = 15*12:length(dumpIters)
% for n = 1:length(dumpIters)
% for n = 1:366
% for n=5*12
% for n=7*12:8*12
% for n = 34a
% for n=48:length(dumpIters)
% for n=2:length(dumpIters)
for n =1:length(dumpIters)
% for n = 400
  dumpIters(n);
    
  t = dumpIters(n)*deltaT/t1year;
  
  
  tyears(n) = t;
  tdays(n) = dumpIters(n)*deltaT/t1day;
  
  suff = '_12hourly';
  T = rdmdsWrapper(fullfile(exppath,'results',['THETA',suff]),dumpIters(n));          
  S = rdmdsWrapper(fullfile(exppath,'results',['SALT',suff]),dumpIters(n));          
  if (isempty(S) || isempty(T))
    error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
  end    
  A = densjmd95(S,T,-RC(k_pot_dens)*gravity*rhonil/1e4*ones(Nx,Ny,Nr))-1000;

  
  if (~isempty(find(isnan(A))))
    error(['Found NaNs at iter=',num2str(dumpIters(n))]);
    break
  end
  
  MLD = zeros(Nx,Ny);
  for i=1:Nx
    for j=1:Ny
      if (sum(hFacC(i,j,:),3)==0)
        continue;
      end

      %%% Extract hFacs for this column
      hFacC_column = squeeze(hFacC(i,j,:));

      %%% Find index of first and last wet grid cells
      idx_top = find(squeeze(hFacC(i,j,:))>0,1,'First');
      if (isempty(idx_top))
        MLD(i,j) = NaN;
        continue;
      end
      idx_bot = find(hFacC_column>0,1,'Last');

      %%% Density of this column
      dens_column = squeeze(A(i,j,:));
      
      %%% Find elevation of top of ocean and bottom of ocean
      z_top = RF(idx_top+1) + hFacC(i,j,idx_top)*DRF(idx_top);
      z_bot = RF(idx_bot) - hFacC(i,j,idx_bot)*DRF(idx_bot);

      
      %%% Reference point is 10m below ocean surface
      z_ref = z_top - 10; 

      %%% Find index of reference point for top of ML
      idx_ref = 0;
      z_tmp = z_top;
      z_prev = z_top;
      z_next = z_top;
      for k = idx_top:idx_bot
        z_prev = z_next;
        if (k == idx_top)
          z_tmp = 0.5*(z_top+RF(k+1));
        elseif (k == idx_bot)
          z_tmp = RF(k) - 0.5*DRF(k)*hFacC_column(k);          
        else
          z_tmp = RC(k);
        end
        if (z_tmp < z_ref)
          idx_ref = k;
          z_next = z_tmp;
          break;
        end
      end
      if (idx_ref == 0)
        error('Could not find k_ref - this should never happen');
      end

      %%% Linearly interpolate to find density 10m below top of ocean
      if (idx_ref == idx_top)
        dens_ref = dens_column(idx_top);
      else
        dens_ref = ( dens_column(idx_ref)*(z_prev-z_ref) + dens_column(idx_ref-1)*(z_ref-z_next) ) / (z_prev-z_next);
      end
      
      %%% ML depth defined by a 0.03 kg/m^3 density difference from 10m
      %%% below top of ocean
      rho_ML = dens_ref + 0.03;
      
      %%% Find first point exceeding the ML depth criterion
      idx_mld = idx_ref - 1 + find(dens_column(idx_ref:idx_bot)>rho_ML,1,'First');

      %%% If there are no such points then the stratification is too weak,
      %%% so the MLD is the water column thickness
      if (isempty(idx_mld))

        z_ML = z_bot;        

      %%% Otherwise, linearly interpolate to find MLD
      else      

        if (idx_mld == idx_bot)
          z_next = 0.5*(RF(idx_mld)+z_bot);
        else
          z_next = RC(idx_mld);
        end

        if ((idx_mld - 1) == idx_top)
          z_prev = 0.5*(RF(idx_mld)+z_top);
        else
          z_prev = RC(idx_mld-1);
        end

        z_ML = z_prev + (rho_ML - dens_column(idx_mld-1)) * (z_next-z_prev) / (dens_column(idx_mld) - dens_column(idx_mld-1));
        
      end

      %%% Compute mixed layer depth
      MLD(i,j) = z_top - z_ML;

    end
  end
    
%     clf;
%     axes('FontSize',fontsize);
  set(gcf,'color','w');
%     FF(FF==0) = NaN;
% %     contourf(XC,YC,FF,100,'EdgeColor','None');  
%     pcolor(XC,YC,FF);
%     shading interp;        
%     hold on;
% %     p = surface(XC,YC,Mice_elev),shading interp;
% %     p.FaceColor = [.5 .5 .5];
% %     p.FaceAlpha = .5;
% %     hold on
% %     [cs,C] = contour(XC,YC,bathy,[-5000:1000:-1000 -500 -200],'EdgeColor','k');
%     
% %     hh = clabel(cs,C);
%     hold off;
%     xlabel('x (km)');
%     ylabel('y (km)');
%     xlabel('Longitude','interpreter','latex');
%     ylabel('Latitude','interpreter','latex');    
%   handle=colorbar;
%   set(handle,'FontSize',fontsize);

  latMin = min(min(YC));
  latMax = YC(1,end-spongethickness);
  lonMin = min(min(XC));
%     lonMax = XC(end-spongethickness,1);   
  lonMax = XC(end,1);
  
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
  pcolorm(YC,XC,MLD);
%     pcolorm(YC,XC,FF/920/35*86400*365);
  shading interp
  
  hold on;
  [cs,C]=contourm(YC,XC,SHELFICEtopo-bathy,[100 500 1000 2000 3000 4000],'EdgeColor','k');
  chandle = clabelm(cs,C);
  set(chandle,'fontsize',12,'Color','k','BackgroundColor','none','Edgecolor','none') 
  hold off;

  %%% Add colorbar and title
  h = colorbar;
  set(gca,'FontSize',14);
  set(h,'Position',[0.92 0.33 0.02 .26])
  tightmap;


 


 
  %%% Finish the plot
  colormap(cmap);
  if (~isempty(titlestr))
%     title([titlestr,', $t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');
    thedate = datenum('2008-01-01')+floor(tdays(n));
%     annotation('textbox',[0.3 0.95 0.5 0.05],'String',[titlestr,', ',months{mod(n-1,12)+1},' ',num2str(years(floor((n-1)/12)+1))],'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');   
    annotation('textbox',[0.4 0.95 0.5 0.05],'String',[titlestr,', ',datestr(thedate)],'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');   
  end
  if (set_crange)  
    caxis(crange);
  end
  
  set(gca,'Position',plotloc);
  set(gca,'FontSize',fontsize);
%   annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  M(n) = getframe(gcf);
  
end





