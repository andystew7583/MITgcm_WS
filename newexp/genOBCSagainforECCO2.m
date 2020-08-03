%%%%%%GenOBCSfiles%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% script to generate OBCS files


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Data format parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ieee='b';
prec='real*8';
realdigits = 8;
realfmt=['%.',num2str(realdigits),'e'];


%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%% SET-UP %%%%%%%%%%%%%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%% Matlab utilities 

addpath ../newexp_utils
addpath ../utils/matlab



%Defining grid that we want data to be interpolated onto
    
defineGrid

% name of SOSE data file

datadir = '/data3/MITgcm_WS/data/ECCO2';



%%% Eastern boundary continental shelf properties
shelf_salt = 34.5;
bathy_max = -400; %%% Limits of bathymetry at eastern boundary over which
bathy_min = -600; %%% to feather modification of shelf temperature







%%%%%
%%%%% Load ECCO2 grid
%%%%%

ref = fullfile(datadir,'SALT.1440x720x50.JAN.nc');

latitude = ncread(ref,'LATITUDE_T');
longitude = ncread(ref,'LONGITUDE_T');

[LAT,LON] = meshgrid(latitude,longitude);

idx1 = find(longitude(:,1)>=180);
idx2 = find(longitude(:,1)<180);

LON = [LON(idx1,:)-360 ; LON(idx2,:)];
LAT = [LAT(idx1,:) ; LAT(idx2,:)];


idx5 = find(LON(:,1)>=xmin & LON(:,1) <=21);
idx6 = find(LAT(1,:)>=ymin & LAT(1,:)<=-64);

lon2 = LON(idx5,idx6);
lat2 = LAT(idx5,idx6);

ECCOlon = lon2(:,1);
ECCOlon = double(ECCOlon);
ECCOlat = lat2(1,:);
ECCOlat = double(ECCOlat);
ECCO_RC = -ncread(ref,'DEPTH_T');
ECCO_RC = double(ECCO_RC)';

SX = size(lon2,1);
SY = size(lat2,2);

Ny_ecco = length(latitude);
Nx_ecco = length(longitude);
Nz_ecco = length(ECCO_RC);









%%% Load climatology on domain boundaries on ECCO grid

[Theta_north,Theta_east,Theta_west] = loadECCOrecs(datadir,'THETA',SX,SY,Nz_ecco,idx1,idx2,idx5,idx6);
[Salt_north,Salt_east,Salt_west] = loadECCOrecs(datadir,'SALT',SX,SY,Nz_ecco,idx1,idx2,idx5,idx6);
[Uvel_north,Uvel_east,Uvel_west] = loadECCOrecs(datadir,'UVEL',SX,SY,Nz_ecco,idx1,idx2,idx5,idx6);
[Vvel_north,Vvel_east,Vvel_west] = loadECCOrecs(datadir,'VVEL',SX,SY,Nz_ecco,idx1,idx2,idx5,idx6);






%%% Interpolate onto our model grid

[ECCO_XX,ECCO_ZZ] = meshgrid(ECCOlon,ECCO_RC);
[ECCO_YY,ECCO_ZY] = meshgrid(ECCOlat,ECCO_RC);



%%%%%%%%%%% OUR GRID %%%%%%%%%%%

[XM,ZM] = meshgrid(xmc,mrc);
[YM,GR] = meshgrid(ymc,mrc);


MPY = 12; %%%% months per year

ModX = Nx;
ModY = Ny;
ModZ = Nr;


%%%%%%%%%%%
%%%%%%------->
%%%% read in bathymetry 


fid = fopen(fullfile(inputconfigdir,bathyFile),'r','b');
h = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);

bathy_east = squeeze(h(end,:));











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Northern Boundary Interp %%%

OBNtEc = interpBdyData2(squeeze(Theta_north),ECCO_XX,ECCO_ZZ,XM,ZM);
OBNsEc = interpBdyData2(squeeze(Salt_north),ECCO_XX,ECCO_ZZ,XM,ZM);
OBNuEc = interpBdyData2(squeeze(Uvel_north),ECCO_XX,ECCO_ZZ,XM,ZM);
OBNvEc = interpBdyData2(squeeze(Vvel_north),ECCO_XX,ECCO_ZZ,XM,ZM);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Eastern Boundary Interp %%%

OBEtEc = interpBdyData2(squeeze(Theta_east),ECCO_YY,ECCO_ZY,YM,GR);
OBEsEc = interpBdyData2(squeeze(Salt_east),ECCO_YY,ECCO_ZY,YM,GR);
OBEuEc = interpBdyData2(squeeze(Uvel_east),ECCO_YY,ECCO_ZY,YM,GR);
OBEvEc = interpBdyData2(squeeze(Vvel_east),ECCO_YY,ECCO_ZY,YM,GR);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Western Boundary Interp %%%

% OBWtEc = interpBdyData2(squeeze(Theta_west),ECCO_YY,ZY,YM,GR);
% OBWsEc = interpBdyData2(squeeze(Salt_west),ECCO_YY,ZY,YM,GR);
% OBWuEc = interpBdyData2(squeeze(Uvel_west),ECCO_YY,ZY,YM,GR);
% OBWvEc = interpBdyData2(squeeze(Vvel_west),ECCO_YY,ZY,YM,GR);





%%% Modify eastern boundary stratification to set properties on continental
%%% shelf
for i = 1:size(OBEtEc,3)
    
  for j=1:Ny
    for k = 1:Nr
      Pressure = -rho0*g*zz(1);
      Tf = freezingTemp(OBEsEc(j,k,i),Pressure);
      
      if (bathy_east(j) >= bathy_max)      
        OBEtEc(j,k,i) = Tf;
        OBEsEc(j,k,i) = shelf_salt;
      else
        if (bathy_east(j) >= bathy_min)
          OBEtEc(j,k,i) = ( Tf * (bathy_east(j)-bathy_min) + OBEtEc(j,k,i) * (bathy_max-bathy_east(j)) ) / (bathy_max - bathy_min);
          OBEsEc(j,k,i) = ( shelf_salt * (bathy_east(j)-bathy_min) + OBEsEc(j,k,i) * (bathy_max-bathy_east(j)) ) / (bathy_max - bathy_min);       
        end
      end
    end
  end    
end








%%%
%%% Write data files
%%%



data = OBNtEc;
writeDataset(data,fullfile(inputconfigdir,OBNtFile),ieee,prec);
clear data



data = OBEtEc;
writeDataset(data,fullfile(inputconfigdir,OBEtFile),ieee,prec);
clear data




data = OBNsEc;
writeDataset(data,fullfile(inputconfigdir,OBNsFile),ieee,prec);
clear data




data = OBEsEc;
writeDataset(data,fullfile(inputconfigdir,OBEsFile),ieee,prec);
clear data




data = OBNuEc;
writeDataset(data,fullfile(inputconfigdir,OBNuFile),ieee,prec);
clear data





data = OBEuEc;
writeDataset(data,fullfile(inputconfigdir,OBEuFile),ieee,prec);
clear data



data = OBNvEc;
writeDataset(data,fullfile(inputconfigdir,OBNvFile),ieee,prec);
clear data




data = OBEvEc;
writeDataset(data,fullfile(inputconfigdir,OBEvFile),ieee,prec);
clear data















%%% 
%%% loadECCOrecs
%%%
%%% Convenience function to load records from SOSE output, switch
%%% longitudinal coordinates, and then retain only a subsection of the grid.
%%%
function [data_north,data_east,data_west] = loadECCOrecs(datadir,varname,Nx,Ny,Nr,lon_idx_west,lon_idx_east,region_lon_idx,region_lat_idx)

  %%% To store boundary data  
  data_north = zeros(Nx,1,Nr,12);
  data_east = zeros(1,Ny,Nr,12);
  data_west = zeros(1,Ny,Nr,12);
  
  %%% Naming convention for ECCO output files by month
  months = {'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'};
  
  %%% Loop through records
  for k = 1:length(months)
    
    %%% Load next record
    data = ncread(fullfile(datadir,[varname,'.1440x720x50.',months{k},'.nc']),varname);
    
    %%% Switch to correct longitude convention
    data = [data(lon_idx_west,:,:) ; data(lon_idx_east,:,:)];
    
    %%% Crop data to region of interest
    data = data(region_lon_idx,region_lat_idx,:);
    
    %%% Extract boundary data
    data_north(:,:,:,k) = data(:,end,:);
    data_east(:,:,:,k) = data(end,:,:);
    data_west(:,:,:,k) = data(1,:,:);
    
  end

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
  data(data<-9e22) = NaN;
  
  %%% To store interpolated data
  data_interp = zeros(size(XI,2),size(XI,1),Nt);
  
  for n = 1:Nt
    data_interp(:,:,n) = interp2(XD,ZD,data(:,:,n)',XI,ZI,'linear')';
    data_interp(:,:,n) = inpaint_nans(data_interp(:,:,n),4); %%% Crude extrapolation
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













