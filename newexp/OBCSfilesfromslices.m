%%%setting initial.bin files for mitgcm parameters for momemtum/tracers





%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%% SET-UP %%%%%%%%%%%%%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%  

defineGrid

        
%%% Matlab utilities 

addpath ../newexp_utils
addpath ../utils/matlab


%%%% Setting min salinity at eastern boundary
min_salt_EB = true;
min_salt_NB = true;
% min_salt = 34;
min_salt = 34.15;

%%% Eastern boundary continental shelf properties
set_shelf_properties = 1;
% shelf_salt = 34;
shelf_salt = 34.15;
bathy_max = -400; %%% Limits of bathymetry at eastern boundary over which
bathy_min = -600; %%% to feather modification of shelf temperature
% bathy_max = -600; %%% Limits of bathymetry at eastern boundary over which
% bathy_min = -1200; %%% to feather modification of shelf temperature

%%% Set true to modify boundary sea ice concentrations
mod_bdry_iceflux = false;
mod_bdry_icethic = false;
mod_bdry_iceconc = false;

%%% Shifting pycnocline at eastern boundary
shift_pyc = 1;
shift_lat_min = -70;
shift_lat_max = -68.5;
shift_depth_max = 200;







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





%%% Modify eastern boundary stratification to set properties on continental
%%% shelf
if (set_shelf_properties)
  for i = 1:size(OBEt,3)

    for j=1:Ny
      for k = 1:Nr
        Pressure = -rho0*g*zz(1);
        Tf = freezingTemp(OBEs(j,k,i),Pressure);

        if (bathy_east(j) >= bathy_max)      
          OBEt(j,k,i) = Tf;
          OBEs(j,k,i) = shelf_salt;
        else
          if (bathy_east(j) >= bathy_min)
            OBEt(j,k,i) = ( Tf * (bathy_east(j)-bathy_min) + OBEt(j,k,i) * (bathy_max-bathy_east(j)) ) / (bathy_max - bathy_min);
            OBEs(j,k,i) = ( shelf_salt * (bathy_east(j)-bathy_min) + OBEs(j,k,i) * (bathy_max-bathy_east(j)) ) / (bathy_max - bathy_min);       
          end
        end
      end
    end    
  end
end

%%% Set minimum salinity and accompanying freezing temperature at eastern boundary
if (min_salt_EB)
  OBEs(OBEs<min_salt) = min_salt;
  OBEt(OBEs==min_salt) = freezingTemp(min_salt,0);
end

%%% Set minimum salinity and accompanying freezing temperature at eastern boundary
if (min_salt_NB)
  OBNs(OBNs<min_salt) = min_salt;
  OBNt(OBNs==min_salt) = freezingTemp(min_salt,0);
end


%%% Shift the stratification at the eastern boundary to deepen the
%%% pycnocline
if (shift_pyc)
  
  %%% Look over monthly climatologies
  for n = 1:size(OBEs,3)
   
    
    %%% Loop through latitudinal indices
    for j = 1:Ny    
      
      %%% Change in depth varies with latitude
      shift_depth = (shift_lat_max - ymc(j)) / (shift_lat_max - shift_lat_min) * shift_depth_max;
      shift_depth = max(min(shift_depth,shift_depth_max),0);
      
      %%% Find index of first gridpoint below the shift depth; everything
      %%% above this gridpoint is simply replaced with the surface T/S
      kmin = find(zz<-(shift_depth+dz(1)/2),1,'first');
      OBEs_tmp = OBEs(j,:,n);
      OBEt_tmp = OBEt(j,:,n);
      OBEs_tmp(1:kmin-1) = OBEs(j,1,n);
      OBEt_tmp(1:kmin-1) = OBEt(j,1,n);
      
      %%% Everything below the shift depth is linearly interpolated from
      %%% the water column properties above
      OBEs_tmp(kmin:end) = interp1(zz,OBEs(j,:,n),zz(kmin:end)+shift_depth,'linear');
      OBEt_tmp(kmin:end) = interp1(zz,OBEt(j,:,n),zz(kmin:end)+shift_depth,'linear');
      
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
    data_interp(:,:,n) = inpaint_nans(data_interp(:,:,n),4); %%% Crude extrapolation
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
