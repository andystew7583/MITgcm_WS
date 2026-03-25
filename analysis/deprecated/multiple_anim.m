%%%%% anim multiple plots in one screen
%%% Reads diagnostic output(s) from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = false;

%%% Read experiment data
loadexp;

%%% Select diagnostic variable to animate
diagnum_1 = 17;
diagnum_2 = 16;
outfname_1 =diag_fileNames{1,diagnum_1};
outfname_2 =diag_fileNames{1,diagnum_2};

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer_1 = 40;  
xylayer_2 = 40;
%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 1;

%%% Set true to plot the field in the topmost wet cell at each horizontal
%%% location
topplot = 0;

%%% Set true to plot the field in the middle of the water column at each
%%% horizontal location
midplot = 0;

%%% Set true for a zonal average
yzavg = 0;

%%% Layer to plot in the y/z plane
yzlayer_2 = 80;
yzlayer_1 = 97;

%%% Specify color range
set_crange = 1;
crange_2 = [-3 1];
% crange = [0 199];
% crange = [0 5]; %%% For sea ice thicknes
% crange_2 = [-0.2 0.2]; %%% For velocities or stresses
% crange = [-1 1]*1e-4; %%% For freshwater fluxes
crange_1 = [34.4 34.9]; %%% salinity
%crange_2 = [-100 100]; %%% Qnet




%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics
dumpFreq = abs(diag_frequency(diagnum_1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((0:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);




%%% Mesh grids for plotting
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
if (~xyplot)  
  %%% Create mesh grid with vertical positions adjusted to sit on the bottom
  %%% topography and at the surface
  [ZZ,YY] = meshgrid(zz,yy);  
  [Z2,Y2] = meshgrid(zz,yy);
  for j=1:Ny
    if (yzavg)
      hFacC_col_1 = squeeze(hFacC(:,j,:));    
      hFacC_col_1 = max(hFacC_col_1,[],1); 
      hFacC_col_2 = squeeze(hFacC(:,j,:));    
      hFacC_col_2 = max(hFacC_col_2,[],1); 
    else
      hFacC_col_1 = squeeze(hFacC(yzlayer_1,j,:))';
      hFacC_col_2 = squeeze(hFacC(yzlayer_2,j,:))';
    end
    zz_topface_1 = zz(kmin(i,j))-(0.5-hFacC_col_1(kmin(i,j)))*delR(kmin(i,j));
    zz_topface_2 = zz(kmin(i,j))-(0.5-hFacC_col_1(kmin(i,j)))*delR(kmin(i,j));
    zz_topface_1 = zz(kmax(i,j))+(0.5-hFacC_col_2(kmax(i,j)))*delR(kmax(i,j));
    zz_botface_2 = zz(kmax(i,j))+(0.5-hFacC_col_2(kmax(i,j)))*delR(kmax(i,j));

    ZZ(j,kmin(i,j)) = zz_topface_1;
    Z2(j,kmin(i,j)) = zz_topface_2;
    
    
    if (kmax(i,j)>1)
      Z2(j,kmax(i,j)) = zz_botface_2;
      ZZ(j,kmax(i,j)) = zz_botface_1;
    end
  end
  
end

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 18;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.15 0.15 0.7 0.76];
  framepos = [100    500   800  800];
end

%%% Set up the figure
fig = figure(1);
set(fig,'Position',plotloc);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');

f = moviein(nDumps);

Amean = [];
tyears = [];
% for n=	1:length(dumpIters)
for n=	1:10
% for n=1:3
dumpIters(n);
    
  t = dumpIters(n)*deltaT/86400/365;
  
  
  tyears(n) = t;
  A = rdmdsWrapper(fullfile(exppath,'results',outfname_1),dumpIters(n));
  if (isempty(A))
    error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
  end     
  
  if (~isempty(find(isnan(A))))
    break
  end
  
  V = rdmdsWrapper(fullfile(exppath,'results',outfname_2),dumpIters(n));
  if (isempty(V))
    error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
  end     
  
  if (~isempty(find(isnan(V))))
    break
  end
  
  
  ['Max value: ',num2str(max(max(max(A(A~=0)))))]
  ['Min value: ',num2str(min(min(min(A(A~=0)))))]
  
  ['Max value: ',num2str(max(max(max(V(V~=0)))))]
  ['Min value: ',num2str(min(min(min(V(V~=0)))))]
       
%   A(hFac==0) = NaN;
  
  
%   DX = repmat(delX',[1 Ny Nr]);
%   DY = repmat(delY,[Nx 1 Nr]);
%   DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);  
%   Amean(n) = sum(sum(sum(A.*DX.*DY.*DZ.*hFacC)))/sum(sum(sum(DX.*DY.*DZ.*hFacC)))
%   
%   Axy = sum(A.*DY.*DZ.*hFacW,3);

   [idx1 idx2] = find(isnan(A));

  
%%% x/y plot
 if (xyplot)
    
    if (botplot)      
      FF = zeros(Nx,Ny);
      VV = zeros(Nx,Ny);

      for i=1:Nx
        for j=1:Ny
          FF(i,j) = A(i,j,kmax(i,j));
          VV(i,j) = V(i,j,kmax(i,j));

        end
      end
    elseif (topplot)
        FF = zeros(Nx,Ny);
        for i=1:Nx
          for j=1:Ny
            FF(i,j) = A(i,j,kmin(i,j));
            VV(i,j) = V(i,j,kmin(i,j));
          end
        end
    elseif (midplot)
          FF = zeros(Nx,Ny);
          for i=1:Nx
            for j=1:Ny
              FF(i,j) = wp(i,j)*A(i,j,kp(i,j)) + wn(i,j)*A(i,j,kn(i,j));
              VV(i,j) = wp(i,j)*V(i,j,kp(i,j)) + wn(i,j)*V(i,j,kn(i,j));

            end
          end
       
     else
          FF = squeeze(A(:,:,xylayer_1,outfidx));  
          VV = squeeze(V(:,:,xylayer_2,outfidx));        
        
     end
        
    
    FF(FF==0) = NaN;
    VV(VV==0) = NaN;
    
    subplot(2,1,1);
    contourf(XC,YC,FF,100,'EdgeColor','None');
%     set(gca,'OuterPosition',[0.55 0.065 0.475 0.28]);
    pcolor(XC,YC,FF);
    shading interp;  
    handle=colorbar;
    colormap jet(100);
    hold on;
    contour(XC,YC,bathy,[-5000:1000:-1000 -500],'EdgeColor','k');
    hold off;
    xlabel('x (km)');
    ylabel('y (km)');
    xlabel('Longitude','interpreter','latex');
    ylabel('Latitude','interpreter','latex');
    if (set_crange)  
     caxis(crange_1);
    end
    title(['$t=',num2str(tyears(n)*365,'%.1f'),'$ days'],'interpreter','latex');

    hold on
 
    subplot(2,1,2);
    contourf(XC,YC,VV,100,'EdgeColor','None');
    pcolor(XC,YC,VV);
    shading interp;  
    colormap jet(100);
    handle=colorbar;
    hold on;
    contour(XC,YC,bathy,[-5000:1000:-1000 -500],'EdgeColor','k');
    hold off;
    xlabel('x (km)');
    ylabel('y (km)');
    xlabel('Longitude','interpreter','latex');
    ylabel('Latitude','interpreter','latex');
    if (set_crange)  
     caxis(crange_2);
    end

    
  else

    
    A(A==0) = NaN;
    V(V==0) = NaN;
    if (yzavg)
      Ayz = squeeze(nanmean(squeeze(A(:,:,:,outfidx))));  
      Ayz_2 = squeeze(nanmean(squeeze(V(:,:,:,outfidx))));    

    else
      Ayz = squeeze(A(yzlayer_1,:,:,outfidx));
      Ayz_2 = squeeze(A(yzlayer_2,:,:,outfidx));

    end
    
    jrange = 1:Ny;

    
    subplot(2,1,1);
    pcolor(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:));
    shading interp;
    hold on;
    [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),10,'EdgeColor','k');
    if (yzavg)
      h = plot(yy,min(bathy,[],1)/1000,'k','LineWidth',3);  
    else
      h = plot(yy,SHELFICEtopo(yzlayer_1,:)/1000,'k','LineWidth',3);  
      h = plot(yy,bathy(yzlayer_1,:)/1000,'k','LineWidth',3);  
    end
    hold off;
    xlabel('Offshore $y$ (km)','interpreter','latex');
    ylabel('Height $z$ (km)','interpreter','latex');
    set(gca,'YLim',[-2 0]);
    
     %%% Finish the plot
    handle=colorbar;
    colormap jet(200);
    set(handle,'FontSize',fontsize);
    title(['$t=',num2str(tyears(n)*365,'%.1f'),'$ days'],'interpreter','latex');
    if (set_crange)  
     caxis(crange_1);
    end
    set(gca,'Position',plotloc);

    set(gca,'FontSize',fontsize);

  
    subplot(2,1,2);
    if (yzavg)
      Ayz_2 = squeeze(nanmean(squeeze(V(:,:,:,outfidx))));    
    else
      Ayz_2 = squeeze(A(yzlayer_2,:,:,outfidx));
    end
      pcolor(Y2(jrange,:),Z2(jrange,:)/1000,Ayz_2(jrange,:));
      shading interp;
      hold on;
      [C,h]=contour(Y2(jrange,:),Z2(jrange,:)/1000,Ayz_2(jrange,:),10,'EdgeColor','k');

    if (yzavg)
       h = plot(yy,min(bathy,[],1)/1000,'k','LineWidth',3);  
    else
       h = plot(yy,SHELFICEtopo(yzlayer_2,:)/1000,'k','LineWidth',3);  
       h = plot(yy,bathy(yzlayer_2,:)/1000,'k','LineWidth',3);  
    end
    hold off;
    xlabel('Offshore $y$ (km)','interpreter','latex');
    ylabel('Height $z$ (km)','interpreter','latex');
    set(gca,'YLim',[-2 0]);
  
    
   %%% Finish the plot
   handle=colorbar;
   colormap jet(200);
   set(handle,'FontSize',fontsize);
   if (set_crange)  
    caxis(crange_2);
   end
   set(gca,'FontSize',fontsize);
 end
 
  f(n) = getframe(gcf);

end

%   close all
%   [h, w, p] = size(f(1).cdata);  % use 1st frame to get dimensions
%   hf = figure; 
% % resize figure based on frame's w x h, and place at (150, 150)
%   set(hf, 'position', [150 150 w h]);
%   axis off
%   movie(hf,f);
%   mplay(f)
  
  
  
  



