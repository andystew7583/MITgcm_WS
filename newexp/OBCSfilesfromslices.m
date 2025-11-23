%%%setting initial.bin files for mitgcm parameters for momemtum/tracers





%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%% SET-UP %%%%%%%%%%%%%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%  

defineGrid

        
%%% Matlab utilities 

addpath ../newexp_utils
addpath ../utils/matlab
addpath ../analysis/

%%% Physical parameters
gravity = 9.81;
rho0 = 1027;
beta = (densjmd95(35,-1.5,1)-densjmd95(34,-1.5,1)) / rho0; %%% Reference haline contraction coefficient


%%%% Setting min salinity at eastern boundary
min_salt_EB = true; %%% Originally set true
min_salt_NB = true; %%% Originally set true
% min_salt = 34;
min_salt = 34.15*ones(1,12); %%% Originally 34.15
% min_salt = [33.9200, 33.5328, 33.5238, 33.7174, 34.0604, 34.0973, 34.0736 34.0862, 34.0818, 34.0774, 34.0728, 34.0680]; %%% Min. salinities from Kapp Norvegia dataset (Hattermann 2018 JPO)
% min_salt = min_salt + 0.2;

%%% Set true to reset temperature to freezing for fresh waters in a way that
%%% depends piecewise-linearly on salinity, approximately consistent
%%% with Kapp Norvegia dataset
reset_temp_linearly_EB = true;
salt_freeze = 34.2;
salt_max = 34.7;
temp_max = 1.5;
temp_freeze = -1.8;
  

%%% Eastern boundary continental shelf properties
set_shelf_properties = true; %%% Originally set true
% shelf_salt = 34;
% shelf_salt = 34.15; %%% Originally set to 34.15, no longer used
% bathy_max = -400; %%% Limits of bathymetry at eastern boundary over which
bathy_min = -800; %%% to feather modification of shelf temperature, originally -400 to -600. Bathy_max no longer used.
% bathy_max = -600; %%% Limits of bathymetry at eastern boundary over which
% bathy_min = -1200; %%% to feather modification of shelf temperature

%%% Set true to modify boundary sea ice concentrations
mod_bdry_iceflux = false;
mod_bdry_icethic = false;
mod_bdry_iceconc = false;

%%% Shifting pycnocline at eastern boundary
shift_pyc = true; %%% Originally set true
shift_lat_min = -71; %%% Original value -70
shift_lat_ref = -70; %%% 11/2025 ALS: Added to extend pyc shift into cavity
shift_lat_max = -68.5; %%% Original value -68.5
shift_depth_max = 200; %%% Original value 200
shift_N2 = 1e-7; %%% Stratification above the shifted pycnocline







%%%%% Input file path for OBCS SOSE generated files
OBCS_storage_dir = './OBCS';

%%% Load monthy-binned SOSE boundary data from .mat file
load(fullfile(OBCS_storage_dir,'OBCS.mat'));






%%%
%%% read in bathymetry and shelf ice
%%%

fid = fopen(fullfile(inputconfigdir,bathyFile),'r','b');
h = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);

fid = fopen(fullfile(inputconfigdir,SHELFICEtopoFile),'r','b');
icedraft = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);

hFacC = zeros(Nx,Ny,Nr);
fid = fopen(fullfile(inputconfigdir,'hFacC.bin'),'r','b');
for k=1:Nr
     hFacC(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


icedraft_east = icedraft(end,:);

bathy_east = h(end,:);

clear h icedraft







%%%%%%%%% SOSE POINTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('sosegrid.mat','XC','YC','RC');
[SOSElonC,SOSElatC] = switchLons (XC,YC,xmin,xmax,ymin,ymax);


[XX,ZZ] = meshgrid(SOSElonC,RC);
[YY,ZY] = meshgrid(SOSElatC,RC);



%%%%%%%%%%% OUR GRID %%%%%%%%%%%

[XM,ZM] = meshgrid(xmc,mrc);
[YM,GR] = meshgrid(ymc,mrc);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Northern Boundary Interp %%%

OBNt = interpBdyData2(squeeze(theta_obcs_north),XX,ZZ,XM,ZM);
OBNs = interpBdyData2(squeeze(salt_obcs_north),XX,ZZ,XM,ZM);
OBNu = interpBdyData2(squeeze(uvel_obcs_north),XX,ZZ,XM,ZM);
OBNv = interpBdyData2(squeeze(vvel_obcs_north),XX,ZZ,XM,ZM);

landidx = find(squeeze(salt_obcs_north(:,:,1,1))==0);
OBNeta = interpBdyData1(squeeze(PHIHYD_obcs_north(:,:,1,:)/g),SOSElonC,xmc,landidx);
OBNh = interpBdyData1(squeeze(SIThick_obcs_north),SOSElonC,xmc,landidx);
OBNa = interpBdyData1(squeeze(SIArea_obcs_north),SOSElonC,xmc,landidx);
OBNuice = interpBdyData1(squeeze(SIUvel_obcs_north),SOSElonC,xmc,landidx);
OBNvice = interpBdyData1(squeeze(SIVvel_obcs_north),SOSElonC,xmc,landidx);
OBNsl = 0*OBNh;








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Eastern Boundary Interp %%%

OBEt = interpBdyData2(squeeze(theta_obcs_east),YY,ZY,YM,GR);
OBEs = interpBdyData2(squeeze(salt_obcs_east),YY,ZY,YM,GR);
OBEu = interpBdyData2(squeeze(uvel_obcs_east),YY,ZY,YM,GR);
OBEv = interpBdyData2(squeeze(vvel_obcs_east),YY,ZY,YM,GR);


landidx = find(squeeze(salt_obcs_east(:,:,1,1))==0);
OBEeta = interpBdyData1(squeeze(PHIHYD_obcs_east(:,:,1,:)/g),SOSElatC,ymc,landidx);
OBEh = interpBdyData1(squeeze(SIThick_obcs_east),SOSElatC,ymc,landidx);
OBEa = interpBdyData1(squeeze(SIArea_obcs_east),SOSElatC,ymc,landidx);
OBEuice = interpBdyData1(squeeze(SIUvel_obcs_east),SOSElatC,ymc,landidx);
OBEvice = interpBdyData1(squeeze(SIVvel_obcs_east),SOSElatC,ymc,landidx);
OBEsl = 0*OBEh;


OBCS_data_dir_snow = '../data/SOSEdata/13-17';
XC = ncread('../data/SOSEdata/13-17/SIhsnow.nc','XC'); %%% 2013-2017 solution is on a different grid
YC = ncread('../data/SOSEdata/13-17/SIhsnow.nc','YC');
XC = repmat(XC,[1 length(YC)]);
YC = repmat(YC',[size(XC,1) 1]);
[SOSElonC,SOSElatC] = switchLons (XC,YC,xmin,xmax,ymin,ymax);
landidx = 1:find(max(squeeze(SIHsnow_obcs_east),2)>0,1,'first')-1;
OBEsn = interpBdyData1(squeeze(SIHsnow_obcs_east),SOSElatC,ymc,landidx);
landidx = 1:find(max(squeeze(SIHsnow_obcs_north),2)>0,1,'first')-1;
OBNsn = interpBdyData1(squeeze(SIHsnow_obcs_north),SOSElonC,xmc,landidx);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Western Boundary Interp %%%

% OBWt = interpBdyData2(squeeze(theta_obcs_west),YY,ZY,YM,GR);
% OBWs = interpBdyData2(squeeze(salt_obcs_west),YY,ZY,YM,GR);
% OBWu = interpBdyData2(squeeze(uvel_obcs_west),YY,ZY,YM,GR);
% OBWv = interpBdyData2(squeeze(vvel_obcs_west),YY,ZY,YM,GR);

% landidx = find(squeeze(salt_obcs_west(:,:,1,1))==0);
% OBWeta = interpBdyData1(squeeze(PHIHYD_obcs_west(:,:,1,:)/g),SOSElatC,ymc,landidx);
% OBWh = interpBdyData1(squeeze(SIThick_obcs_west),SOSElatC,ymc,landidx);
% OBWa = interpBdyData1(squeeze(SIArea_obcs_west),SOSElatC,ymc,landidx);
% OBWuice = interpBdyData1(squeeze(SIUvel_obcs_west),SOSElatC,ymc,landidx);
% OBWvice = interpBdyData1(squeeze(SIVvel_obcs_west),SOSElatC,ymc,landidx);
% OBWsl = 0*OBWh;
% OBWsn = 0.3*ones(size(OBWh));





%%% Set minimum salinity and accompanying freezing temperature at eastern boundary
if (min_salt_EB)

  for i=1:size(OBEs,3)
    OBEs_tmp = OBEs(:,:,i);
    OBEs_tmp(OBEs_tmp<min_salt(i)) = min_salt(i);
    OBEs(:,:,i) = OBEs_tmp;
  end

  %%% Confines the temperature below a line in T/S space, inspired by the
  %%% Kapp Norvegia water masses from Hattermann 2018 JPO
  if (reset_temp_linearly_EB)
    OBEt(OBEs <= salt_freeze) = temp_freeze;
    saltcond = (OBEs > salt_freeze) & (OBEs < salt_max); 
    temp_limit = temp_freeze + (temp_max-temp_freeze)*(OBEs - salt_freeze)/(salt_max-salt_freeze);
    OBEt(saltcond & (OBEt > temp_limit)) = temp_limit(saltcond & (OBEt > temp_limit));
  else
    % for i=1:size(OBEt,3)
    %   OBEt_tmp = OBEt(:,:,i);
    %   OBEs_tmp = OBEs(:,:,i);
    %   OBEt_tmp(OBEs_tmp==min_salt(i)) = freezingTemp(min_salt(i),0);
    %   OBEt(:,:,i) = OBEt_tmp;
    % end
  end

end



%%% Shift the stratification at the eastern boundary to deepen the
%%% pycnocline
if (shift_pyc)
  
  %%% Look over monthly climatologies
  for n = 1:size(OBEs,3)
   
    
    %%% Loop through latitudinal indices
    for j = 1:Ny    
      
      %%% Change in depth varies with latitude
      shift_depth = (shift_lat_max - ymc(j)) / (shift_lat_max - shift_lat_ref) * shift_depth_max;
      shift_depth = max(min(shift_depth,shift_depth_max),0);
      
      %%% Find index of first gridpoint below the shift depth;
      kmin = find(zz<-(shift_depth+dz(1)/2),1,'first');
      OBEs_tmp = OBEs(j,:,n);
      OBEt_tmp = OBEt(j,:,n);

      %%% Everything below the shift depth is linearly interpolated from
      %%% the water column properties above
      OBEs_tmp(kmin:end) = interp1(zz,OBEs(j,:,n),zz(kmin:end)+shift_depth,'linear');
      OBEt_tmp(kmin:end) = interp1(zz,OBEt(j,:,n),zz(kmin:end)+shift_depth,'linear');

      %%% If there is unstable stratification, adjust kmin downward (at
      %%% most 10 times, arbitrarily)
      for ktmp = 1:10
        N2loc = -(g/rho0)*(densjmd95(OBEs_tmp(kmin),OBEt_tmp(kmin),-zzf(kmin+1)) - densjmd95(OBEs_tmp(kmin+1),OBEt_tmp(kmin+1),-zzf(kmin+1)))/(zz(kmin)-zz(kmin+1));
        if (N2loc < shift_N2)
          kmin = kmin + 1;
        else
          break;
        end
      end

      %%% Everything above this gridpoint is simply replaced a linear
      %%% salinity stratification and uniform temperature
      OBEs_tmp(1:kmin-1) = OBEs_tmp(kmin) - shift_N2/(gravity*beta)*(zz(1:kmin-1)-zz(kmin));
      OBEt_tmp(1:kmin-1) = OBEt_tmp(kmin);     
      
      %%% Replace OB properties
      OBEs(j,:,n) = OBEs_tmp;
      OBEt(j,:,n) = OBEt_tmp;
      
    end
    
    %%% Density (approximate) for thermal wind calculation
    OBEd = densmdjwf(OBEs(:,:,n),OBEt(:,:,n),repmat(-zz,[Ny 1]));
      
    %%% Recalculate zonal velocity so that it's in thermal wind balance  
    for j = 1:Ny    
      
      %%% Index of deepest wet grid cell
      kbot = find(squeeze(hFacC(end,j,:))>0,1,'last');
      
      %%% Calculate thermal wind where we have shifted the pycnocline
      if (~isempty(kbot) && ymc(j)<shift_lat_max)
        OBEu_tmp = OBEu(j,:,n);
        f0 = 2*Omega*sind(ymc(j));
        deltay = (ymc(j+1)-ymc(j-1))* 2*pi/360 * Rp;
        for k = kbot-1:-1:1 %%% Preserve velocity in lowest grid cell
          OBEu_tmp(k) = OBEu_tmp(k+1) + (zz(k)-zz(k+1)) *  g/(rho0*f0) * (OBEd(j+1,k)-OBEd(j-1,k)) / deltay;
        end
        OBEu(j,:,n) = OBEu_tmp;
      end
        
    end
  end
  
end


%%% Modify eastern boundary stratification to set properties on continental
%%% shelf
%%% 11/2025 ALS: Moved to after shift_pyc code so cavity properties derive
%%% from shifted pycnocline
if (set_shelf_properties)
  % for i = 1:size(OBEt,3)

    % for j=1:Ny
    %   for k = 1:Nr
      %   Pressure = -rho0*g*zz(1);
      %   Tf = freezingTemp(OBEs(j,k,i),Pressure);
      % 
      %   if (bathy_east(j) >= bathy_max)      
      %     OBEt(j,k,i) = Tf;
      %     OBEs(j,k,i) = shelf_salt;
      %   else
      %     if (bathy_east(j) >= bathy_min)
      %       OBEt(j,k,i) = ( Tf * (bathy_east(j)-bathy_min) + OBEt(j,k,i) * (bathy_max-bathy_east(j)) ) / (bathy_max - bathy_min);
      %       OBEs(j,k,i) = ( shelf_salt * (bathy_east(j)-bathy_min) + OBEs(j,k,i) * (bathy_max-bathy_east(j)) ) / (bathy_max - bathy_min);       
      %     end
      %   end
      % end
      
    %   end
    % end    
  % end

  %%% Uniform extrapolation into the cavity for depths greater than
  %%% bathy_min
  jshelf = find((icedraft_east>bathy_east) & (bathy_east>=bathy_min),1,'last');
  OBEt(1:jshelf,:,:) = repmat(OBEt(jshelf+1,:,:),[length(1:jshelf) 1 1]);
  OBEs(1:jshelf,:,:) = repmat(OBEs(jshelf+1,:,:),[length(1:jshelf) 1 1]);

  %%% Replace any salinities lower than the lowest salinities outside the
  %%% cavity with the lowest salinity outside the cavity. This avoids very
  %%% low salinities appearing in the cavity
  jicefront = find(icedraft_east<0,1,'last')+1;
  for i=1:size(OBEt,3)
    OBEs_tmp = OBEs(:,:,i);
    OBEs_tmp(OBEs_tmp<OBEs_tmp(jicefront,1)) = OBEs_tmp(jicefront,1);
    OBEs(:,:,i) = OBEs_tmp;
  end
end






%%% Set minimum salinity and accompanying freezing temperature at northern boundary
if (min_salt_NB)

   for i=1:size(OBNs,3)
    OBNt_tmp = OBNt(:,:,i);
    OBNs_tmp = OBNs(:,:,i);
    OBNs_tmp(OBNs_tmp<min_salt(i)) = min_salt(i);   
    OBNt_tmp(OBNs_tmp==min_salt(i)) = freezingTemp(min_salt(i),0);
    OBNs(:,:,i) = OBNs_tmp;
    OBNt(:,:,i) = OBNt_tmp;
   end

end




%%% Modify eastern boundary sea ice velocity to control sea ice influx 
if (mod_bdry_iceflux)
  OBEuice = OBEuice*2; %%% Doubled ice speed comes closer to interior ice speeds
end
if (mod_bdry_icethic)
  OBEa(:,:) = 0.9; %%% Uniform ice concentration along eastern boundary
  OBNa(:,:) = 0.9; %%% Uniform ice concentration along northern boundary
end  
if (mod_bdry_iceconc)
  OBEh(:,:) = 1; %%% Uniform ice thickness along eastern boundary
  OBNh(:,:) = 1; %%% Uniform ice thickness along northern boundary
end
    
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    WRITING FILES %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




data = OBNt;
writeDataset(data,fullfile(inputconfigdir,OBNtFile),ieee,prec);
clear data



data = OBEt;
writeDataset(data,fullfile(inputconfigdir,OBEtFile),ieee,prec);
clear data



data = OBNs;
writeDataset(data,fullfile(inputconfigdir,OBNsFile),ieee,prec);
clear data


data = OBEs;
writeDataset(data,fullfile(inputconfigdir,OBEsFile),ieee,prec);
clear data



data = OBNu;
writeDataset(data,fullfile(inputconfigdir,OBNuFile),ieee,prec);
clear data



data = OBEu;
writeDataset(data,fullfile(inputconfigdir,OBEuFile),ieee,prec);
clear data


data = OBNv;
writeDataset(data,fullfile(inputconfigdir,OBNvFile),ieee,prec);
clear data



data = OBEv;
writeDataset(data,fullfile(inputconfigdir,OBEvFile),ieee,prec);
clear data




data = OBNa;
writeDataset(data,fullfile(inputconfigdir,OBNaFile),ieee,prec);
clear data


data = OBEa;
writeDataset(data,fullfile(inputconfigdir,OBEaFile),ieee,prec);
clear data




data = OBNh;
writeDataset(data,fullfile(inputconfigdir,OBNhFile),ieee,prec);
clear data



data = OBEh;
writeDataset(data,fullfile(inputconfigdir,OBEhFile),ieee,prec);
clear data




data = OBNsn;
writeDataset(data,fullfile(inputconfigdir,OBNsnFile),ieee,prec);
clear data


data = OBEsn;
writeDataset(data,fullfile(inputconfigdir,OBEsnFile),ieee,prec);
clear data




data = OBNuice;
writeDataset(data,fullfile(inputconfigdir,OBNuiceFile),ieee,prec);
clear data


data = OBEuice;
writeDataset(data,fullfile(inputconfigdir,OBEuiceFile),ieee,prec);
clear data




data = OBNvice;
writeDataset(data,fullfile(inputconfigdir,OBNviceFile),ieee,prec);
clear data



data = OBEvice;
writeDataset(data,fullfile(inputconfigdir,OBEviceFile),ieee,prec);
clear data




data = OBNsl;
writeDataset(data,fullfile(inputconfigdir,OBNslFile),ieee,prec);
clear data



data = OBEsl;
writeDataset(data,fullfile(inputconfigdir,OBEslFile),ieee,prec);
clear data




data = OBNeta;
writeDataset(data,fullfile(inputconfigdir,OBNetaFile),ieee,prec);
clear data



data = OBEeta;
writeDataset(data,fullfile(inputconfigdir,OBEetaFile),ieee,prec);
clear data













%%% Plotting options
framepos = [1000         142        1018        1072];
fontsize = 13;
[ZZ,YY] = meshgrid(zz,ymc);

%%% Generate model plots
for m = 1:size(OBEt,3)

  OBEs_plot = OBEs(:,:,m);
  OBEs_plot(squeeze(hFacC(end,:,:))==0) = NaN;

  %%% Plot salinity section
  figure(11);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  pcolor(YY(:,:),-ZZ(:,:)/1000,OBEs_plot);
  shading interp;

  set(gca,'YDir','reverse');
  caxis([33.6 34.7]);
  if (m >= 10)
    xlabel('Latitude','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  set(gca,'XLim',[-71 -67]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('haline',20));
    set(handle,'Position',[.93 .05 .01 .9]);
  end
  


  OBEt_plot = OBEt(:,:,m);
  OBEt_plot(squeeze(hFacC(end,:,:))==0) = NaN;

  %%% Plot temperature section
  figure(12);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  pcolor(YY(:,:),-ZZ(:,:)/1000,OBEt_plot(:,:));
  shading interp;
  set(gca,'YDir','reverse');
  caxis([-2.2 1.2])
  if (m >= 10)
    xlabel('Latitude','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  set(gca,'XLim',[-71 -67]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('thermal',20));
    set(handle,'Position',[.93 .05 .01 .9]);
  end

  %%% Compute density gradient
  drhodz = zeros(Ny,Nr);
  for k = 1:Nr
    DENS = densjmd95(OBEs_plot(:,k),OBEt_plot(:,k),-zzf(k+1));
    DENS(squeeze(hFacC(end,:,k))==0) = NaN;
    if (k < Nr)     
      DENS_below = densjmd95(OBEs_plot(:,k+1),OBEt_plot(:,k+1),-zzf(k+1));       
      DENS_below(squeeze(hFacC(end,:,k+1))==0) = NaN;
      drhodz(:,k) = (DENS - DENS_below) / (zz(k) - zz(k+1));
    end      
  end

  %%% Plot buoyancy frequency
  figure(13);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  pcolor(YY(:,:),-ZZ(:,:)/1000,-gravity*drhodz/rho0);
  shading interp;
  set(gca,'YDir','reverse');
  caxis([-1 1]*1e-5)
  if (m >= 10)
    xlabel('Latitude','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  set(gca,'XLim',[-71 -67]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(redblue(20));
    set(handle,'Position',[.93 .05 .01 .9]);
  end




  jrange = 1:find(ymc>-68.5,1);
  krange = 1:find(zz<-1000,1);

  %%% T/S diagram
  figure(14);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  scatter(reshape(OBEs_plot(jrange,krange),1,[]),reshape(OBEt_plot(jrange,krange),1,[])); 
  if (m >= 10)
    xlabel('Salinity (g/kg)','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Temperature ($^\circ$C)','interpreter','latex','fontsize',fontsize);
  end
  axis([33 35 -2.2 1.2]);

end





%%% 
%%% switchLons
%%%
%%% Convenience function to switch grids from SOSE's [0,360] convention to
%%% our [-180,180] convention, generate indices for switching other
%%% matrices, and generate indices that restrict SOSE data to our model
%%% domain.
%%%
function [SOSElonC,SOSElatC] = switchLons (XC,YC,xmin,xmax,ymin,ymax)

  %%% Indices of "western" and "eastern" halves of SOSE domain
  idx_west = find(XC(:,1)>=180);
  idx_east = find(XC(:,1)<180);

  %%% Switch longitude and latitude grids around
  XC = [XC(idx_west,:)-360 ; XC(idx_east,:)];
  YC = [YC(idx_west,:) ; YC(idx_east,:)];

  %%% Indices defining the subset of SOSE that contains our model grid
  idx_OBN = find(XC(:,1)>xmin & XC(:,1)<xmax);
  idx_OBE = find(YC(1,:)>ymin & YC(1,:)<ymax);
  
  %%% Restrict grid matrices to our model domain
  XC = XC(idx_OBN,idx_OBE);
  YC = YC(idx_OBN,idx_OBE);
  
  %%% Extract 1D grid vectors
  SOSElonC = XC(:,1);
  SOSElatC= YC(1,:);

end


%%%
%%% interpBdyData2
%%%
%%% Convenience function to interpolate boundary data onto the model grid;
%%%
function data_interp = interpBdyData2 (data,XD,ZD,XI,ZI)
 
  %%% Input grid size  
  Nt = size(data,3);
  
  %%% Remove land points
  data(data==0) = NaN;
  
  %%% To store interpolated data
  data_interp = zeros(size(XI,2),size(XI,1),Nt);

  
  for n = 1:Nt
    data_interp(:,:,n) = interp2(XD,ZD,data(:,:,n)',XI,ZI,'linear')';
    data_interp(:,:,n) = inpaint_nans(data_interp(:,:,n),0); %%% Crude extrapolation %%% 11/2025 ALS: Changed from method 4 to method 0
  end

end


%%%
%%% interpBdyData1
%%%
%%% Convenience function to interpolate boundary data onto the model grid.
%%%
function data_interp = interpBdyData1 (data,xd,xi,landidx)
 
  %%% Input grid size  
  Nt = size(data,2);
  
  %%% Remove land points
  xd(landidx) = [];
  
  %%% To store interpolated data
  data_interp = zeros(length(xi),Nt);
  
  for n = 1:Nt  
    data_tmp = data(:,n);
    data_tmp(landidx) = [];
    data_interp(:,n) = interp1(xd,data_tmp,xi,'nearest','extrap')';
  end

end



%%% 
%%% freezingTemp
%%%
%%% Calculates freezing temperature
%%%
function Tfreeze = freezingTemp (salt,pressure)

  Pa1dbar = 1e4;
  Tfreeze = .0901 - .0575*salt - (7.61e-4 *(pressure/Pa1dbar));                
  
end

