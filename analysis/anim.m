%%%
%%% anim.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = false;

%%% Read experiment data
% loadexp;

%%% Select diagnostic variable to animate
% diagnum = 75;
diagnum = 2;
% diagnum = 7;
outfname =diag_fileNames{1,diagnum};

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 1;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 0;

%%% Set true to plot the field in the topmost wet cell at each horizontal
%%% location
topplot = 0;

%%% Set true to plot the field in the middle of the water column at each
%%% horizontal location
midplot = 0;

%%% Set true for a zonal average
yzavg = 0;

%%% Layer to plot in the y/z plane
%%%for 1/3 DEGREE yzlayer = 126; 1/6 = 404;
yzlayer = 600;

% load ../newexp/ELEV.mat
% 
% 
% Mice_elev=Mice_elev';


%%% Specify color range
set_crange = 1;


% crange = [-2.2 -1.6]; %/%% Filchner temp
% crange = [-3 1]; %/%%temp
% crange = [34.2 35.0]; %%% salinity
% crange = [0 10]; %%%% for KPP hbl
% crange = [0 1]; %%% For sea ice area
% crange = [-.6 .6]; %%% For velocities or stresses
% crange = [-1 1]*1e-4; %%% For freshwater fluxes
% crange =[-100 100]; %%% Qnet
% crange = [-300 300]; %%% swnet
crange = [0 3];
% crange = [-0.2 0.2];
% crange = [0.2 1.2]; %% SSH

% cmap = pmkmp(100,'Swtth');
% cmap = cmocean('haline',100);
% cmap = cmocean('thermal',100);
% cmap = cmocean('ice',100);
% cmap = haxby;
% cmap = jet(200);
% cmap = redblue(100);
cmap = cmocean('amp',100);

% titlestr = 'Bottom salinity (g/kg)';
% titlestr = 'Sea ice concentration';
titlestr = '';

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
years = 2007:1:2015;


%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics
% nIter0 = 1278720;
% deltaT = 200
% nIter0 = 587520;
% nIter0 = 394509; %%% For 1/24 run with tides
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

%%% Mesh grids for plotting
if (botplot || topplot || midplot || ~xyplot)
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
end
if (~xyplot)  
  %%% Create mesh grid with vertical positions adjusted to sit on the bottom
  %%% topography and at the surface
  [ZZ,YY] = meshgrid(zz,yy);  
  for j=1:Ny
    if (yzavg)
      hFacC_col = squeeze(hFacC(:,j,:));    
      hFacC_col = max(hFacC_col,[],1);    
    else
      hFacC_col = squeeze(hFacC(yzlayer,j,:))';
    end
    zz_topface = zz(kmin(i,j))-(0.5-hFacC_col(kmin(i,j)))*delR(kmin(i,j));
    zz_botface = zz(kmax(i,j))+(0.5-hFacC_col(kmax(i,j)))*delR(kmax(i,j));
    ZZ(j,kmin(i,j)) = zz_topface;
    if (kmax(i,j)>1)
      ZZ(j,kmax(i,j)) = zz_botface;
    end
  end
  
end

% for i = 1:Nx
%     for j = 1:Ny
%         if SHELFICEtopo(i,j)-bathy(i,j)<=0
%             
%             Mice_elev(i,j)=NaN;
%             SHELFICEtopo(i,j)=NaN;
%             bathy(i,j)=NaN;
%         end
%     end
% end
% 
% bathy(y<0)=NaN;
% SHELFICEtopo(SHELFICEtopo==0) = NaN;
% Mice_elev(Mice_elev==0) = NaN;


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
  framepos = [100   462   810   520];
end

%%% Set up the figure
handle = figure(22);
set(handle,'Position',framepos);
M = moviein(length(dumpIters));

Amean = [];
Amax = [];

% for n = 15*12:length(dumpIters)
% for n = 1:length(dumpIters)
% for n = 20
% for n=5*12
% for n=7*12:8*12
% for n = 34
% for n=48:length(dumpIters)
% for n=2:length(dumpIters)
for n =700:length(dumpIters)
  dumpIters(n);
    
  t = dumpIters(n)*deltaT/t1year;
  
  
  tyears(n) = t;
  A = rdmdsWrapper(fullfile(exppath,'results',outfname),dumpIters(n));
  if (isempty(A))
    continue;
%     error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
  end     
  
  if (~isempty(find(isnan(A))))
    error(['Found NaNs at iter=',num2str(dumpIters(n))]);
    break
  end
  
  
  
  

   [idx1 idx2] = find(isnan(A));

   if (size(A,3)>1)
     A(hFacC==0) = NaN;
   end
   
   Amax(n)= nanmax(nanmax(A(:,:,1)));
   ['Max value: ',num2str(nanmax(A(:)))]
   ['Min value: ',num2str(nanmin(A(:)))]
   
%%% x/y plot
 if (xyplot)
    
    if (botplot)      
      FF = zeros(Nx,Ny);
      for i=1:Nx
        for j=1:Ny
          FF(i,j) = A(i,j,kmax(i,j));
        end
      end
      
    elseif (topplot)
        FF = zeros(Nx,Ny);
        for i=1:Nx
          for j=1:Ny
            FF(i,j) = A(i,j,kmin(i,j));
          end
        
        end
    elseif (midplot)
          FF = zeros(Nx,Ny);
          for i=1:Nx
            for j=1:Ny
              FF(i,j) = wp(i,j)*A(i,j,kp(i,j)) + wn(i,j)*A(i,j,kn(i,j));
            end
          end
        
    else
          FF = squeeze(A(:,:,xylayer,outfidx));        
          FF(hFacC(:,:,xylayer)==0) = NaN;
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
    pcolorm(YC,XC,FF);
    shading interp

    %%% Add colorbar and title
    h = colorbar;
    set(gca,'FontSize',10);
    set(h,'Position',[0.92 0.33 0.02 .26])
    tightmap;



    
  %%% y/z zonally-averaged plot
 
 else
    
    
%     if (yzavg)
      Ayz = ((squeeze(A(yzlayer,:,:))));    
%     else
%       Ayz = squeeze(A(yzlayer,:,:,outfidx));
%     end
%     Ayz =  squeeze(Ayz(jrange,:,:));
    Ayz(Ayz==0)=NaN;
    jrange = 1:Ny;
    [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz,200,'EdgeColor','None');              
%     pcolor(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:));
    shading interp;
    hold on;
    [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz,20,'EdgeColor','k');
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[-2:0.5:12],'EdgeColor','k');
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[-0.1:0.005:0.1],'EdgeColor','k');
%     if (yzavg)
%       h = plot(yy,min(bathy,[],1)/1000,'k','LineWidth',3);  
%     else
      h = plot(yy,SHELFICEtopo(yzlayer,:)/1000,'k','LineWidth',3);  
      h = plot(yy,bathy(yzlayer,:)/1000,'k','LineWidth',3);  
%     end
    hold off;
    xlabel('Offshore $y$ (latitude)','interpreter','latex','fontsize',15);
    ylabel('Height $z$ (km)','interpreter','latex','fontsize',15);
    set(gca,'YLim',[-1.3 0]);
    set(gca,'XLim',[-83 -64]);    
    handle=colorbar;
    set(handle,'FontSize',fontsize);

 end

 


 
  %%% Finish the plot
  colormap(cmap);
  if (~isempty(titlestr))
%     title([titlestr,', $t=',num2str(tyears(n),'%.1f'),'$ years'],'interpreter','latex');
    annotation('textbox',[0.3 0.95 0.5 0.05],'String',[titlestr,', ',months{mod(n-1,12)+1},' ',num2str(years(floor((n-1)/12)+1))],'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');   
  end
  if (set_crange)  
    caxis(crange);
  end
  
  set(gca,'Position',plotloc);
  set(gca,'FontSize',fontsize);
%   annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  M(n) = getframe(gcf);
  
end





