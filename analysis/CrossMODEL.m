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
loadexp;

%%% Select diagnostic variable to animate
diagnum = 12;
outfname =diag_fileNames{1,diagnum};
dumpFreq = abs(diag_frequency(15));
deltaT = 400;
% nIter0 = 0;
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);


%%% Data index in the output data files
outfidx = 1;
%%%%%%%%%%%
%%%%%NEXT
Nyear_s =19;
Nyear_f = 27;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;
Nyears = Nyear_f-Nyear_s;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Saltt = NaN(Nx,Ny,Nr,nDump_Finish-nDump_start);
Tempt = NaN(Nx,Ny,Nr,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  SALT = rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n+nDump_start));
  Saltt(:,:,:,n) = (SALT);
  Temp = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n+nDump_start));
  Tempt(:,:,:,n) = (Temp);
end

SALT_season = zeros(size(Saltt,1),size(Saltt,2),Nr,Nyears,4);
Temp_season = zeros(size(Tempt,1),size(Tempt,2),Nr,Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,:,n,1) = mean(Saltt(:,:,:,(n-1)*12+[12 1 2]),4);
  SALT_season(:,:,:,n,2) = mean(Saltt(:,:,:,(n-1)*12+(3:5)),4);
  SALT_season(:,:,:,n,3) = mean(Saltt(:,:,:,(n-1)*12+(6:8)),4);
  SALT_season(:,:,:,n,4) = mean(Saltt(:,:,:,(n-1)*12+(9:11)),4);
  Temp_season(:,:,:,n,1) = mean(Tempt(:,:,:,(n-1)*12+[12 1 2]),4);
  Temp_season(:,:,:,n,2) = mean(Tempt(:,:,:,(n-1)*12+(3:5)),4);
  Temp_season(:,:,:,n,3) = mean(Tempt(:,:,:,(n-1)*12+(6:8)),4);
  Temp_season(:,:,:,n,4) = mean(Tempt(:,:,:,(n-1)*12+(9:11)),4);
end
TEMP=squeeze(mean(Temp_season,4));
TEMP=TEMP(:,:,:,3);
SALT=squeeze(mean(SALT_season,4));
SALT=SALT(:,:,:,3);
%%% Set true for a zonal average
yzavg = 1;

%%% Layer to plot in the y/z plane

yzlayer = 82;

set_crange = 1;
crange = [-2 .5]; %%%temp
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
  
  
Ayz = (nanmean(squeeze(TEMP(:,:,:,outfidx))));    

Ayz2 = squeeze(nanmean(squeeze(SALT(:,:,:,outfidx))));   




% SALT=squeeze(SALT(yzlayer,:,:));
% SALT(SALT==0)=NaN;
% TEMP = squeeze(TEMP(yzlayer,:,:));
% TEMP(TEMP==0)=NaN;


%%% Calculate potential density

pdc = densmdjwf(SALT,TEMP,500*ones(Nx,Ny,Nr)) - 1000;
pdc(pdc<20)=0;
pdc = squeeze(pdc(yzlayer,:,:));
pdc(pdc==0)=NaN;
hold on

% [C,f] = contour(SS_gridc,PT_gridc,pdc,32.0:.1:35.0,'EdgeColor','k');





    jrange = 1:Ny;
%     [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,TEMP,200,'EdgeColor','None');              
%     pcolor(YY,ZZ/1000,Ayz);
%     shading interp;
    hold on;
    [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,pdc,[25:.04:43],'EdgeColor','k');
    clabel(C,h)
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[-2:0.5:12],'EdgeColor','k');
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[-0.1:0.005:0.1],'EdgeColor','k');
%     if (yzavg)
%       h = plot(yy,min(bathy,[],1)/1000,'k','LineWidth',3);  
%     else
      h = plot(yy,SHELFICEtopo(yzlayer,:)/1000,'k','LineWidth',3);  
      h = plot(yy,bathy(yzlayer,:)/1000,'k','LineWidth',3);  
%     end
    hold off;
    xlabel('Latitude ($^\circ$S)','interpreter','latex','fontsize',15);
    ylabel('Height $z$ (km)','interpreter','latex','fontsize',15);
    set(gca,'YLim',[-1 0]);
  
 hold on
 


 
  %%% Finish the plot
%   handle=colorbar;
%   colormap jet(200);
%   crange = [-3 .5];
%   set(handle,'FontSize',15);
% %   title(['$t=',num2str(tyears(n)*365,'%.1f'),'$ days'],'interpreter','latex');
%   if (set_crange)  
%     caxis(crange);
%   end
%   
%   set(gca,'Position',plotloc);
  set(gca,'FontSize',16); 
  
 title('Potential Density referenced to 500m, 55W','fontsize',14,'interpreter','latex');

  
  
  
  