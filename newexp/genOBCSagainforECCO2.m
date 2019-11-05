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



%Defining grid that we want data to be interpolated onto
    
defineGrid

% name of SOSE data file

datadir = '/data1/MITgcm_WS/data/ECCO2';



% Switch to [-180,180] longitude range for SOSE Data

%%%%% find dimensions

ref = fullfile(datadir,'SALT.1440x720x50.JAN.nc');

latitude = ncread(ref,'LATITUDE_T');
longitude = ncread(ref,'LONGITUDE_T');

[LAT,LON] = meshgrid(latitude,longitude);

idx1 = find(longitude(:,1)>=180);
idx2 = find(longitude(:,1)<180);

LON = [LON(idx1,:)-360 ; LON(idx2,:)];
LAT = [LAT(idx1,:) ; LAT(idx2,:)];


idx5 = find(LON(:,1)>-83 & LON(:,1) <21);
idx6 = find(LAT(1,:)>-84 & LAT(1,:)<-64);

lon2 = LON(idx5,idx6);
lat2 = LAT(idx5,idx6);

ECCOlon = lon2(:,1);
ECCOlon = double(ECCOlon);
ECCOlat = lat2(1,:);
ECCOlat = double(ECCOlat);
RC = -ncread(ref,'DEPTH_T');
RC = double(RC)';

%%% Matlab utilities 

addpath ../newexp_utils
addpath ../utils/matlab



SX = size(lon2,1);
SY = size(lat2,2);



X_num = 84;
Y_num = 80;


Ny_ecco = 720;
Nx_ecco = size(longitude,1);
Nz = 50;


%loop through total number of records to get ocean state


%%%%%%% SALT %%%%%%%
months = {'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'};

Salt_month = NaN(Nx_ecco,Ny_ecco,Nz,12);
for k = 1:size(months,2)
        file = fullfile(datadir,['SALT.1440x720x50.',sprintf('%s',(months{k}),'.nc')]);
        Salt = ncread(file,'SALT');
        Salt_month(:,:,:,k) = Salt;
end


Uvel_month = NaN(Nx_ecco,Ny_ecco,Nz,12);
for k = 1:size(months,2)
        Uvel = ncread(fullfile(datadir,['UVEL.1440x720x50.',sprintf('%s',(months{k}),'.nc')]),'UVEL');
        Uvel_month(:,:,:,k) = Uvel;
end

Vvel_month = NaN(Nx_ecco,Ny_ecco,Nz,12);
for k = 1:size(months,2)
        Vvel = ncread(fullfile(datadir,['VVEL.1440x720x50.',sprintf('%s',(months{k}),'.nc')]),'VVEL');
        Vvel_month(:,:,:,k) = Vvel;
end
    
THETA_month = NaN(Nx_ecco,Ny_ecco,Nz,12);
for k = 1:size(months,2)
        THETA = ncread(fullfile(datadir,['THETA.1440x720x50.',sprintf('%s',(months{k}),'.nc')]),'THETA');
        THETA_month(:,:,:,k) = THETA;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Extracting Slices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
THETA_month = [THETA_month(idx1,:,:,:) ; THETA_month(idx2,:,:,:)];
Theta = THETA_month(idx5,idx6,:,:); 
Theta_east = Theta(end,:,:,:);
Theta_north = Theta(:,end,:,:);
Theta_west = Theta(1,:,:,:);
    


Salt_month = [Salt_month(idx1,:,:,:) ; Salt_month(idx2,:,:,:)];
Salt = Salt_month(idx5,idx6,:,:);
Salt_north = Salt(:,end,:,:);
Salt_east = Salt(end,:,:,:);
Salt_west = Salt(1,:,:,:);
    

Uvel_month = [Uvel_month(idx1,:,:,:) ; Uvel_month(idx2,:,:,:)];
Uvel = Uvel_month(idx5,idx6,:,:);
Uvel_north = Uvel(:,end,:,:);
Uvel_east = Uvel(end,:,:,:);
Uvel_west = Uvel(1,:,:,:);
    
    

Vvel_month = [Vvel_month(idx1,:,:,:) ; Vvel_month(idx2,:,:,:)];
Vvel = Vvel_month(idx5,idx6,:,:);
Vvel_north = Vvel(:,end,:,:);
Vvel_east = Vvel(end,:,:,:);
Vvel_west = Vvel(1,:,:,:);




%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Now the Open Boundary Interpolation should be pretty straight forward
%%%%%%%%%%%%%%%%%%%%%%%%%

[XX,ZZ] = meshgrid(ECCOlon,RC);
[YY,ZY] = meshgrid(ECCOlat,RC);



%%%%%%%%%%% OUR GRID %%%%%%%%%%%

[XM,ZM] = meshgrid(xmc,mrc);
[YM,GR] = meshgrid(ymc,mrc);


MPY = 12; %%%% months per year

ModX = Nx;
ModY = Ny;
ModZ = Nr;


%%%%%%%%%%%
%%%%%%------->
%%%% read in bathymetry and shelf ice


fid = fopen(fullfile(inputfolder,bathyFile),'r','b');
h = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);

fid = fopen(fullfile(inputfolder,SHELFICEtopoFile),'r','b');
icedraft = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);


hydrogTheta = zeros(Nx,Ny,Nr);
fid = fopen(fullfile(inputfolder,hydrogThetaFile),'r','b');
for k=1:Nr
     hydrogTheta(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);

hFacC = zeros(Nx,Ny,Nr);
fid = fopen(fullfile(inputfolder,'hFacC.bin'),'r','b');
for k=1:Nr
     hFacC(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);

hydrogSalt = zeros(Nx,Ny,Nr);
fid = fopen(fullfile(inputfolder,hydrogSaltFile),'r','b');
for k=1:Nr
      hydrogSalt(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);

salt_e = squeeze(hydrogSalt(end,:,:));
theta_e = squeeze(hydrogTheta(end,:,:));


icedraft_east = icedraft(end,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Northern Boundary Interp %%%

%%%%%%%%% Temperature

OBNtEc = NaN(ModX,ModZ,MPY);
Theta_north = squeeze(Theta_north);
Theta_north(Theta_north<-9e+22) = NaN; 
Theta_north(Theta_north==0)=NaN;
TN_temp = zeros(size(Theta_north,1),size(Theta_north,2));
for i = 1:size(Theta_north,3)
    TN_temp = Theta_north(:,:,i);
    OBNtEc(:,:,i) = interp2(XX,ZZ,TN_temp',XM',ZM');
    OBNtEc(:,:,i) = inpaint_nans(OBNtEc(:,:,i),4);
end



OBNsEc = NaN(ModX,ModZ,MPY);
Salt_north = squeeze(Salt_north);
Salt_north(Salt_north<-9e+22) = NaN; 
Salt_north(Salt_north==0)=NaN;
SN_temp = zeros(size(Salt_north,1),size(Salt_north,2));
for i = 1:size(Salt_north,3)
    SN_temp = Salt_north(:,:,i);
    OBNsEc(:,:,i) = interp2(XX,ZZ,SN_temp',XM',ZM');
    OBNsEc(:,:,i) = inpaint_nans(OBNsEc(:,:,i),4);
end



OBNuEc = NaN(ModX,ModZ,MPY);
Uvel_north = squeeze(Uvel_north);
Uvel_north(Uvel_north<-9e+22) = NaN; 
Uvel_north(Uvel_north==0)=NaN;
UN_temp = zeros(size(Uvel_north,1),size(Uvel_north,2));
for i = 1:size(Uvel_north,3)
    UN_temp = Uvel_north(:,:,i);
    OBNuEc(:,:,i) = interp2(XX,ZZ,UN_temp',XM',ZM');
    OBNuEc(:,:,i) = inpaint_nans(OBNuEc(:,:,i),4);
end



OBNvEc = NaN(ModX,ModZ,MPY);
Vvel_north = squeeze(Vvel_north);
Vvel_north(Vvel_north<-9e+22) = NaN;
Vvel_north(Vvel_north==0)=NaN;
VN_temp = zeros(size(Vvel_north,1),size(Vvel_north,2));
for i = 1:size(Vvel_north,3)
    VN_temp = Vvel_north(:,:,i);
    OBNvEc(:,:,i) = interp2(XX,ZZ,VN_temp',XM',ZM');
    OBNvEc(:,:,i) = inpaint_nans(OBNvEc(:,:,i),4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%eASterNNNNBoundArYYY%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Eastern Boundary Interp %%%

%%%%%%%%% Temperature

% OBEtEc = NaN(ModY,ModZ,MPY);
% Theta_east = squeeze(Theta_east);
% Theta_east(Theta_east<-9e+22) = NaN; 
% Theta_east(Theta_east==0)=NaN;
% TE_temp = zeros(size(Theta_east,1),size(Theta_east,2));
% for i = 1:size(Theta_east,3)
%     TE_temp = Theta_east(:,:,i);
%     OBEtEc(:,:,i) = interp2(YY,ZY,TE_temp',YM',GR');
%     OBEtEc(:,:,i) = inpaint_nans(OBEtEc(:,:,i),4);
% end
% 
% 
% 
% 
% OBEsEc = NaN(ModY,ModZ,MPY);
% Salt_east = squeeze(Salt_east);
% Salt_east(Salt_east<-9e+22) = NaN; 
% Salt_east(Salt_east==0)=NaN;
% SE_temp = zeros(size(Salt_east,1),size(Salt_east,2));
% for i = 1:size(Salt_east,3)
%     SE_temp = Salt_east(:,:,i);
%     OBEsEc(:,:,i) = interp2(YY,ZY,SE_temp',YM',GR');
%     OBEsEc(:,:,i) = inpaint_nans(OBEsEc(:,:,i),4);
% end


OBEtEc = NaN(Ny,Nr,MPY);
Theta_east = squeeze(Theta_east);
hFacC_e = squeeze(hFacC(end,:,:));
Theta_east(Theta_east<-9e+22)=NaN;
TE_temp = zeros(size(Theta_east,1),size(Theta_east,2));
for i = 1:size(Theta_east,3)
    TE_temp = Theta_east(:,:,i);
    OBEtEc(:,:,i) = interp2(YY,ZY,TE_temp',YM',GR');
    
        for j=1:Ny
            for k=1:Nr
        
                 if zz(k)<icedraft_east(j)&&(sum(hFacC_e(j,k+1:end))) > 0
                    OBEtEc(j,k,i) = OBEtEc(j,k,i);
                 else
                  Pressure = -rho0*g*zz(k);
                  OBEtEc(j,k,i) = .0901 - .0575*salt_e(j,k) - (7.61e-4 *(Pressure/Pa1dbar));
                 end
            end
        end
    OBEtEc(:,:,i) = inpaint_nans(OBEtEc(:,:,i),4);
end








OBEsEc = NaN(Ny,Nr,MPY);
Salt_east = squeeze(Salt_east);
hFacC_e = squeeze(hFacC(end,:,:));
Salt_east(Salt_east<-9e+22)=NaN;
SE_temp = zeros(size(Salt_east,1),size(Salt_east,2));
for i = 1:size(Salt_east,3)
    SE_temp = Salt_east(:,:,i);
    OBEsEc(:,:,i) = interp2(YY,ZY,SE_temp',YM',GR');
        for j=1:Ny
            for k=1:Nr
        
                 if zz(k)<icedraft_east(j)&&(sum(hFacC_e(j,k+1:end))) > 0
                    OBEsEc(j,k,i) = OBEsEc(j,k,i);
                 else
                    Pressure = -rho0*g*zz(k);
                    OBEsEc(j,k,i) = 34;
                 end
            end
        end
    OBEsEc(:,:,i) = inpaint_nans(OBEsEc(:,:,i),4);
    
    
end




OBEuEc = NaN(ModY,ModZ,MPY);
Uvel_east = squeeze(Uvel_east);
Uvel_east(Uvel_east<-9e+22) = NaN; 
Uvel_east(Uvel_east==0)=NaN;
UE_temp = zeros(size(Uvel_east,1),size(Uvel_east,2));
for i = 1:size(Uvel_east,3)
    UE_temp = Uvel_east(:,:,i);
    OBEuEc(:,:,i) = interp2(YY,ZY,UE_temp',YM',GR');
    OBEuEc(:,:,i) = inpaint_nans(OBEuEc(:,:,i),4);
end




OBEvEc = NaN(ModY,ModZ,MPY);
Vvel_east = squeeze(Vvel_east);
Vvel_east(Vvel_east<-9e+22) = NaN; 
Vvel_east(Vvel_east==0)=NaN;
VE_temp = zeros(size(Vvel_east,1),size(Vvel_east,2));
for i = 1:size(Vvel_east,3)
    VE_temp = Vvel_east(:,:,i);
    OBEvEc(:,:,i) = interp2(YY,ZY,VE_temp',YM',GR');
    OBEvEc(:,:,i) = inpaint_nans(OBEvEc(:,:,i),4);
end



data = OBNtEc;
writeDataset(data,fullfile(inputfolder,OBNtFile),ieee,prec);
clear data



data = OBEtEc;
writeDataset(data,fullfile(inputfolder,OBEtFile),ieee,prec);
clear data



data = OBNsEc;
writeDataset(data,fullfile(inputfolder,OBNsFile),ieee,prec);
clear data




data = OBEsEc;
writeDataset(data,fullfile(inputfolder,OBEsFile),ieee,prec);
clear data




data = OBNuEc;
writeDataset(data,fullfile(inputfolder,OBNuFile),ieee,prec);
clear data





data = OBEuEc;
writeDataset(data,fullfile(inputfolder,OBEuFile),ieee,prec);
clear data



data = OBNvEc;
writeDataset(data,fullfile(inputfolder,OBNvFile),ieee,prec);
clear data




data = OBEvEc;
writeDataset(data,fullfile(inputfolder,OBEvFile),ieee,prec);
clear data



















