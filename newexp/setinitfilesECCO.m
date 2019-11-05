%%%%%%%%%%%%%%
%%% setinitfiles for ECCO2
%%%%
%%%% initial conditions for ecco2

%%%setting initial.bin files for mitgcm parameters for momemtum/tracers

%%%10/13: trying to interpolate each grid point to local freezing point


%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%% SET-UP %%%%%%%%%%%%%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%  
  

%Defining grid that we want data to be interpolated onto

defineGrid

% name of SOSE data file

datadir = '/data3/MITgcm_WS/data/ECCO2';



%%% Matlab utilities 

addpath ../newexp_utils
addpath ../utils/matlab

%define record numbeR

recnum = 12;

%%%% Run Topography File

fid = fopen(fullfile(inputfolder,bathyFile),'r','b');
h = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);

fid = fopen(fullfile(inputfolder,SHELFICEtopoFile),'r','b');
icedraft = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);



%%% Limit values of dependent variables
Theta_max = 3;
Theta_min = -3;
Salt_max = 37;
Salt_min = 33;
Uvel_max = 0.1;
Uvel_min = -0.1;
Vvel_max = 0.1;
Vvel_min = -0.1;

%%% Salinity of ice shelf water - used to set initial salinity under shelf
%%% ice
S_ISW = 34.9;
UVEL_ISW = 0.0;
VVEL_ISW = 0.0;




%%%%%%%%%%%%%%%%%%%%%%
%%%% LOAD ECCO2 DATA   
%%%%%%%%%%%%%%%%%%%%%%
ECCOlatC  = (-90+1/8) : (1/4) :  90;
ECCOlonC =     (1/8) : (1/4) : 360;
dz_ecco = [10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01, ...
 10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04 , 19.82, 24.85,... 
 31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18, ...
 93.96, 96.58, 98.25, 99.25,100.01,101.33,104.56,111.33,122.83,...
 139.09,158.94,180.83,203.55,226.50,249.50,272.50,295.50,318.50,...
 341.50,364.50,387.50,410.50,433.50,456.50];

zz_ecco = -cumsum((dz_ecco+[0 dz_ecco(1:end-1)])/2);


Ny_ecco = 720;
Nx_ecco = 1440;
NZ_ECCO = 50;

%%%%%%% data %%%%%%%
months = {'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'};

Salt_month = NaN(Nx_ecco,Ny_ecco,NZ_ECCO,12);
for k = 1:size(12)
        file = fullfile(datadir,['SALT.1440x720x50.',sprintf('%s',(months{k}),'.nc')]);
        Salt_month = ncread(file,'SALT');
        Salt_month(:,:,:,k) = Salt_month;
end
Salt_month(Salt_month<-9e+22) = NaN;

Salt_month = (nanmean(Salt_month,4));


Uvel_month = NaN(Nx_ecco,Ny_ecco,NZ_ECCO,12);
for k = 1:12
        Uvel_month = ncread(fullfile(datadir,['UVEL.1440x720x50.',sprintf('%s',(months{k}),'.nc')]),'UVEL');
        Uvel_month(:,:,:,k) = Uvel_month;
end
Uvel_month = nanmean(Uvel_month,4);
Uvel_month(Uvel_month<-9e+22) = NaN;

Vvel_month = NaN(Nx_ecco,Ny_ecco,NZ_ECCO,12);
for k = 1:12
        Vvel_month = ncread(fullfile(datadir,['VVEL.1440x720x50.',sprintf('%s',(months{k}),'.nc')]),'VVEL');
        Vvel_month(:,:,:,k) = Vvel_month;
end
Vvel_month = nanmean(Vvel_month,4);
Vvel_month(Vvel_month<-9e+22)= NaN;

    
Theta_month = NaN(Nx_ecco,Ny_ecco,NZ_ECCO,12);
for k = 1:12
        THETA = ncread(fullfile(datadir,['THETA.1440x720x50.',sprintf('%s',(months{k}),'.nc')]),'THETA');
        Theta_month(:,:,:,k) = THETA;
end


Theta_month(Theta_month<-9e+22)=NaN;
Theta_month = nanmean(Theta_month,4);



% Switch to [-180,180] longitude range from ECCO2 Data
[XC,YC] = meshgrid(ECCOlonC,ECCOlatC);
YC = YC';
XC = XC';
idx1 = find(XC(:,1)>=180);
idx2 = find(XC(:,1)<180);

XC = [XC(idx1,:)-360 ; XC(idx2,:)];
YC = [YC(idx1,:) ; YC(idx2,:)];




%%%%% In define Grid
% %%%%%% for tracers ---> 
ECCOlon = XC(:,1);
ECCOlat= YC(1,:);



Theta_month = [Theta_month(idx1,:,:) ; Theta_month(idx2,:,:)];
Salt_month = [Salt_month(idx1,:,:) ; Salt_month(idx2,:,:)];
Uvel_month = [Uvel_month(idx1,:,:) ; Uvel_month(idx2,:,:)];
Vvel_month = [Vvel_month(idx1,:,:) ; Vvel_month(idx2,:,:)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% TWIST FUNCTION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Remove topography



%%%%%%------->


TwistedTheta = Tweddell(ECCOlon,ECCOlat,Theta_month);
TwistedSalt = Tweddell(ECCOlon,ECCOlat,Salt_month);
TwistedUvel = Tweddell(ECCOlon,ECCOlat,Uvel_month);
TwistedVvel = Tweddell(ECCOlon,ECCOlat,Vvel_month);


% %%%%%% meshgrids %%%%%%%%%%%%%%%%%%%%%%%
[XSC,YSC,RSC]=meshgrid(ECCOlon,ECCOlat,zz_ecco);


%%%%%%%%%%%%% This means we have to flip y,x indices now.
 
 
TwistedTheta = transpose3D(TwistedTheta);
TwistedSalt = transpose3D(TwistedSalt);
TwistedUvel = transpose3D(TwistedUvel);
TwistedVvel = transpose3D(TwistedVvel);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ITERPOLATE DATA ONTO OUR GRID %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%% Twisted Salinity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TwistedSaltPretzal = interp3(XSC,YSC,RSC,(TwistedSalt),Xm3C,Ym3C,RM3C,'linear'); 
for k=1:Nz
%   salt_slice = TwistedSaltPretzal(:,:,k);
%   salt_slice(icedraft<0) = S_ISW;
%   TwistedSaltPretzal(:,:,k) = salt_slice;
end
for k = 1:size(TwistedSaltPretzal,3)
  TwistedSaltPretzal(:,:,k)=inpaint_nans(TwistedSaltPretzal(:,:,k),4);
end
TwistedSaltPretzal(TwistedSaltPretzal==0) = NaN;
for i=1:Nx
  Salt_slice = inpaint_nans(squeeze(TwistedSaltPretzal(:,i,:)),4);
  TwistedSaltPretzal(:,i,:) = reshape(Salt_slice,[Ny 1 Nz]); 
end
TwistedSaltPretzal(TwistedSaltPretzal==0) = NaN;
for k = 1:size(TwistedSaltPretzal,3)
  TwistedSaltPretzal(:,:,k)=inpaint_nans(TwistedSaltPretzal(:,:,k),4);
end

TwistedSaltPretzal=transpose3D(TwistedSaltPretzal);
% TwistedSaltPretzal(isnan(TwistedSaltPretzal)) = 0;
% TwistedSaltPretzal(TwistedSaltPretzal == 0) = Salt_max;
% 
% TwistedSaltPretzal(TwistedSaltPretzal>Salt_max) = Salt_max;
% TwistedSaltPretzal(TwistedSaltPretzal<Salt_min) = Salt_min;


%%%%%%%%%%%%%%%% UVELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Please use XG,YC,RC as SOSE original grid (west point) 

%%%%%%%%%%%%%%%%%%% Twisted UVEL %%%%%%%%%%%%%%%%%%%%%%%%%%


MTwistUvel = interp3(XSC,YSC,RSC,(TwistedUvel),Xm3C,Ym3C,RM3C,'linear'); 
for k=1:Nz
%   UVel_slice = MTwistUvel(:,:,k);
%   UVel_slice(icedraft<0) = UVEL_ISW;
%   MTwistUvel(:,:,k) = UVel_slice;
end
for k = 1:size(MTwistUvel,3)
  MTwistUvel(:,:,k)=inpaint_nans(MTwistUvel(:,:,k),4);
end
MTwistUvel(MTwistUvel==0) = NaN;
for i=1:Nx
  UVel_slice = inpaint_nans(squeeze(MTwistUvel(:,i,:)),4);
  MTwistUvel(:,i,:) = reshape(UVel_slice,[Ny 1 Nz]); 
end
MTwistUvel(MTwistUvel==0) = NaN;
for k = 1:size(MTwistUvel,3)
  MTwistUvel(:,:,k)=inpaint_nans(MTwistUvel(:,:,k),4);
end

MTwistUvel=transpose3D(MTwistUvel);
% MTwistUvel(isnan(MTwistUvel)) = 0; %%% NaNs will blow up the model
% MTwistUvel(MTwistUvel>Uvel_max) = Uvel_max;
% MTwistUvel(MTwistUvel<Uvel_min) = Uvel_min;

%%%%%%%%%%%%%%% VVELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PlZ use XC,YG,RC as original SOSE grid (on South Point) %%

%%%%%%%%%%%%%%%%%%% Twisted VVEL %%%%%%%%%%%%%%%%%%%%%%%%%%


MTwistVvel = interp3(XSC,YSC,RSC,(TwistedVvel),Xm3C,Ym3C,RM3C,'linear'); 
% for k=1:Nz
%   VVel_slice = MTwistVvel(:,:,k);
%   VVel_slice(icedraft<0) = VVEL_ISW;
%   MTwistVvel(:,:,k) = VVel_slice;
% end
for k = 1:size(MTwistVvel,3)
  MTwistVvel(:,:,k)=inpaint_nans(MTwistVvel(:,:,k),4);
end
MTwistVvel(MTwistVvel==0) = NaN;
for i=1:Nx
  VVel_slice = inpaint_nans(squeeze(MTwistVvel(:,i,:)),4);
  MTwistVvel(:,i,:) = reshape(VVel_slice,[Ny 1 Nz]); 
end
MTwistVvel(MTwistVvel==0) = NaN;
for k = 1:size(MTwistVvel,3)
  MTwistVvel(:,:,k)=inpaint_nans(MTwistVvel(:,:,k),4);
end

MTwistVvel=transpose3D(MTwistVvel);
% MTwistVvel(isnan(MTwistVvel)) = 0; %%% NaNs will blow up the model
% MTwistVvel(MTwistVvel>Vvel_max) = Vvel_max;
% MTwistVvel(MTwistVvel<Vvel_min) = Vvel_min;

%%%%%%%%%%% VELOCITY VECTOR TWIST %%%%%%%%%%%%%%%%%%%%%%%%%%


% MTwistVelvec = interp3(XM3C,YM3C,R3C,(MTwistedVel_Vec),XUmG,YUmG,RM3C); 

% MTwistVelvec = reshape(MTwistVelvec,([Nx,Ny,Nz]));

%%%%%%%%
%%%%%%%%%%



% Surface Pressure Anomaly---------------->>>>>>>>>>> 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pressure


% ModPressure(isnan(ModPressure)) = 0;
Modelsurfacepressure = zeros(size(TwistedSaltPretzal,1),size(TwistedSaltPretzal,2),1);

Twisted_SSH_anom = (Modelsurfacepressure/g);      

%%% Write initial model surface height anomaly to filename.  Try making it
%%% zero....Initiaize with zero deviations from free surface
Modelsurfaceheight = zeros(size(TwistedSaltPretzal,1),size(TwistedSaltPretzal,2));



%%%%%%%%%%% accounting for supercooled water under the ice shelves
%%%%%%%%%%% TEMPERATURE %%%%%%%%%%%%%%%%%

%%%%%%%%%% Twisted TEMPERATURE  %%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% MTwistTemp = NaN(Ny,Nx,NZ_SOSE);
% thetatwist_temp=zeros(size(TwistedTheta,1),size(TwistedTheta,2));
% for k = 1:size(TwistedTheta,3);
%       thetatwist_temp = TwistedTheta(:,:,k);
%       MTwistTemp(:,:,k)=interp2(XC',YC',(thetatwist_temp)',XMC,YMC,'nearest');           
%       MTwistTemp(:,:,k) = inpaint_nans(MTwistTemp(:,:,k),4);
% end


ModelTwistedTheta = interp3(XSC,YSC,RSC,(TwistedTheta),Xm3C,Ym3C,RM3C,'linear');
% Theta_slice = ModelTwistedTheta(:,:,k);
% for i = 1:Nx
%     for j = 1:Ny
%      SP = (rho0*Modelsurfacepressure(i,j))/Pa1dbar;
% 
%       if (icedraft(i,j)<0)
%             for k = 1:Nz
%                 Pressure = (-rho0*g*zz(k)/Pa1dbar);
%                 ModelTwistedTheta(:,:,k) = .0901 - .0575*TwistedSaltPretzal(i,j,k) - (7.61e-4 *(SP+Pressure));
%             end
%       else
%           ModelTwistedTheta(:,:,k) = Theta_slice;
% 
%       end
%     end
% end
    
        
for k = 1:size(ModelTwistedTheta,3)
  ModelTwistedTheta(:,:,k)=inpaint_nans(ModelTwistedTheta(:,:,k),4);
end
ModelTwistedTheta(ModelTwistedTheta==0) = NaN;
for i=1:Nx
  Theta_slice = inpaint_nans(squeeze(ModelTwistedTheta(:,i,:)),4);
  ModelTwistedTheta(:,i,:) = reshape(Theta_slice,[Ny 1 Nz]); 
end
ModelTwistedTheta(ModelTwistedTheta==0) = NaN;
for k = 1:size(ModelTwistedTheta,3)
  ModelTwistedTheta(:,:,k)=inpaint_nans(ModelTwistedTheta(:,:,k),4);
end

ModelTwistedTheta=transpose3D(ModelTwistedTheta);



ModelTwistedTheta(ModelTwistedTheta>Theta_max) = Theta_max;
ModelTwistedTheta(ModelTwistedTheta<Theta_min) = Theta_min;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ice Shelf Load Anomaly 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Equation is Ptop-g*Σρ0*delz       where ρ0 ~ 1027.5kg/m^3
%%% Ptop ~ Pa+ g*Σρ**delz


%%% These are token values, and need not correspond to what is used in the
%%% model


%%%%%% Trying to substract a purely depth dependent contribution  -g(rho0)dz
%%%%%% with a constant reference density from total pressure.  Anomaly is
%%%%%% zero where there is no shelf ice, because there is no rhoShelfIce

MIce_Anom = zeros(Nx,Ny);
for i=1:Nx
  for j=1:Ny

      MIce_Anom(i,j) = rho0*Modelsurfacepressure(i,j);

    for k=1:Nr
      if (zz(k) < icedraft(i,j))
          break
      end
%       if icedraft(i,j) == 0
%           MIce_Anom(i,j) = MIce_Anom(i,j) ;
%       else
%           rhoShelfIce = densmdjwf(TwistedSaltPretzal(i,j,k),ModelTwistedTheta(i,j,k),(rho0*ModPressure(i,j,k)));
      Pressure = -rho0*g*zz(k);
      rhoShelfIce = densmdjwf(TwistedSaltPretzal(i,j),ModelTwistedTheta(i,j),Pressure/Pa1dbar);
%       MIce_Anom(i,j) = ModPressure(i,j,1) + (g*(rhoShelfIce-rho0)*dz(k))/rho0;
      MIce_Anom(i,j) = MIce_Anom(i,j) + (g*(rhoShelfIce-rho0)*dz(k));                

      %%% N.B. this does not account for wet/dry grid cells properly
%           MIce_Anom(i,j) = MIce_Anom(i,j) + g*(rhoShelfIce-rho0)*dz(k);
     end
   end
end







  


data = MIce_Anom;
writeDataset(data,fullfile(inputfolder,SHELFICEloadAnomalyFile),ieee,prec);
clear data











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WRITING FILES TO INPUT PATH FOR INTITIALIZATION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = ModelTwistedTheta;
writeDataset(data,fullfile(inputfolder,hydrogThetaFile),ieee,prec);
clear data


data = TwistedSaltPretzal;
writeDataset(data,fullfile(inputfolder,hydrogSaltFile),ieee,prec);
clear data


data = MTwistUvel;
writeDataset(data,fullfile(inputfolder,uVelInitFile),ieee,prec);
clear data


data = MTwistVvel;
writeDataset(data,fullfile(inputfolder,vVelInitFile),ieee,prec);
clear data


data = Modelsurfaceheight;
writeDataset(data,fullfile(inputfolder,pSurfInitFile),ieee,prec);
clear data






