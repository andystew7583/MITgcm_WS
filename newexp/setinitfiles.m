%%%
%%% setinitfiles.m
%%%
%%% Creates initial .bin files for mitgcm parameters for momemtum/tracers
%%%
%%%
%%% 10/13: trying to interpolate each grid point to local freezing point
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%% SET-UP %%%%%%%%%%%%%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%  
  

%%%Defining grid that we want data to be interpolated onto
defineGrid

%%% Grid for SOSE variables
load sosegrid.mat;

%%% name of SOSE data file
datadir = '/data3/MITgcm_WS/data/SOSEdata/08-10';

%%% Matlab utilities 
addpath ../newexp_utils
addpath ../utils/matlab

%define record numbeR

recnum = 12;

%%%% Run Topography File

fid = fopen(fullfile(inputconfigdir,bathyFile),'r','b');
h = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);


fid = fopen(fullfile(inputconfigdir,SHELFICEtopoFile),'r','b');
icedraft = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);


%%% Limit values of dependent variables
Theta_max = 3;
Theta_min = -4;
Salt_max = 37;
Salt_min = 32;
Uvel_max = 0.1;
Uvel_min = -0.2;
Vvel_max = 0.2;
Vvel_min = -0.2;




%%%%%%%%%%%%%%%%%%%%%%
%%%% LOAD SOSE DATA   
%%%%%%%%%%%%%%%%%%%%%%


Theta = rdmds(fullfile(datadir,'THETA'),'rec',recnum);
Salt = rdmds(fullfile(datadir,'SALT'),'rec',recnum);
Uvel = rdmds(fullfile(datadir,'UVEL'),'rec',recnum);
Vvel = rdmds(fullfile(datadir,'VVEL'),'rec',recnum);
PHIHYD = rdmds(fullfile(datadir,'PHIHYD'),'rec',recnum);


% Switch from [-180,180] longitude range for SOSE Data

idx1 = find(XC(:,1)>=180);
idx2 = find(XC(:,1)<180);
idx3 = find(XG(:,1)>=180);
idx4 = find(XG(:,1)<180);
XC = [XC(idx1,:)-360 ; XC(idx2,:)];
YC = [YC(idx1,:) ; YC(idx2,:)];
XG = [XG(idx3,:)-360 ; XG(idx4,:)];
YG = [YG(idx3,:) ; YG(idx4,:)];




%Change longitude change of respective variable field

Theta = [Theta(idx1,:,:) ; Theta(idx2,:,:)];
Salt = [Salt(idx1,:,:) ; Salt(idx2,:,:)];
Uvel = [Uvel(idx1,:,:) ; Uvel(idx2,:,:)];
Vvel = [Vvel(idx1,:,:) ; Vvel(idx2,:,:)];
PHIHYD = [PHIHYD(idx1,:,:) ; PHIHYD(idx2,:,:)];




%%%%% Change lat/lon of partial cell fields

hFacC = [hFacC(idx1,:,:) ; hFacC(idx2,:,:)];
hFacW = [hFacW(idx1,:,:) ; hFacW(idx2,:,:)];
hFacS = [hFacS(idx1,:,:) ; hFacS(idx2,:,:)];






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% TWIST FUNCTION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lat/lon variable for Tweddell function for SOSE

%%%%% In define Grid
% %%%%%% for tracers ---> 
SOSElonC = XC(:,1);
SOSElatC= YC(1,:);


%%%%%% for velocity points ----> 
SOSElonG= XG(:,1);
SOSElatG= YG(1,:);

NZ_SOSE = size(Theta,3);
%%% Remove topography

PHIHYD(Salt==0) = NaN;
Theta(Theta==0) = NaN; 
Salt(Salt==0) = NaN;
Uvel(Uvel==0) = NaN;
Vvel(Vvel==0) = NaN;


%%%%%------->


TwistedTheta = Tweddell(SOSElonC,SOSElatC,Theta);
TwistedSalt = Tweddell(SOSElonC,SOSElatC,Salt);
TwistedUvel = Tweddell(SOSElonG,SOSElatC,Uvel);
TwistedVvel = Tweddell(SOSElonC,SOSElatG,Vvel);


% %%%%%% Straight up 3D interp!!!!!!!!! %%%%%%%%%%%%%%%%%%%%%%%
[XSC,YSC,RSC]=meshgrid(SOSElonC,SOSElatC,RC);
%%%% Uvel
[XSU,YSU,RSU]=meshgrid(SOSElonG,SOSElatC,RC);
%%%%% Vvel
[XSV,YSV,RSV]=meshgrid(SOSElonC,SOSElatG,RC);



%%%%%%%%%%%% This means we have to flip y,x indices now.
 
 
TwistedTheta = transpose3D(TwistedTheta);
TwistedSalt = transpose3D(TwistedSalt);
TwistedUvel = transpose3D(TwistedUvel);
TwistedVvel = transpose3D(TwistedVvel);









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INTERPOLATE DATA ONTO OUR GRID %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Twisted Salinity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


TwistedSaltPretzal = interp3(XSC,YSC,RSC,(TwistedSalt),Xm3C,Ym3C,RM3C,'linear'); 
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

TwistedSaltPretzal(isnan(TwistedSaltPretzal)) = 0;
TwistedSaltPretzal(TwistedSaltPretzal == 0) = Salt_max;

TwistedSaltPretzal(TwistedSaltPretzal>Salt_max) = Salt_max;
TwistedSaltPretzal(TwistedSaltPretzal<Salt_min) = Salt_min;






%%%%%%%%%%%%%%%% UVELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Please use XG,YC,RC as SOSE original grid (west point) 

%%%%%%%%%%%%%%%%%%% Twisted UVEL %%%%%%%%%%%%%%%%%%%%%%%%%%


MTwistUvel = interp3(XSC,YSC,RSC,(TwistedUvel),Xm3C,Ym3C,RM3C,'linear'); 
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

MTwistUvel(isnan(MTwistUvel)) = 0; %%% NaNs will blow up the model
MTwistUvel(MTwistUvel>Uvel_max) = Uvel_max;
MTwistUvel(MTwistUvel<Uvel_min) = Uvel_min;







%%%%%%%%%%%%%% VVELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PlZ use XC,YG,RC as original SOSE grid (on South Point) %%

%%%%%%%%%%%%%%%%%% Twisted VVEL %%%%%%%%%%%%%%%%%%%%%%%%%%


MTwistVvel = interp3(XSC,YSC,RSC,(TwistedVvel),Xm3C,Ym3C,RM3C,'linear'); 
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
MTwistVvel(isnan(MTwistVvel)) = 0; %%% NaNs will blow up the model
MTwistVvel(MTwistVvel>Vvel_max) = Vvel_max;
MTwistVvel(MTwistVvel<Vvel_min) = Vvel_min;








% Surface Pressure Anomaly---------------->>>>>>>>>>> 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pressure

pressure_temp = PHIHYD;
Twisted_Pressure = Tweddell(SOSElonC,SOSElatC,pressure_temp);
Twisted_Pressure=transpose3D(Twisted_Pressure);


ModPressure = interp3(XSC,YSC,RSC,(Twisted_Pressure),Xm3C,Ym3C,RM3C,'linear'); 
for k = 1:size(ModPressure,3)
  ModPressure(:,:,k)=inpaint_nans(ModPressure(:,:,k),4);
end
ModPressure(ModPressure==0) = NaN;
for i=1:Nx
  Pres_slice = inpaint_nans(squeeze(ModPressure(:,i,:)),4);
  ModPressure(:,i,:) = reshape(Pres_slice,[Ny 1 Nz]); 
end
ModPressure(ModPressure==0) = NaN;
for k = 1:size(ModPressure,3)
  ModPressure(:,:,k)=inpaint_nans(ModPressure(:,:,k),4);
end

ModPressure=transpose3D(ModPressure);
% ModPressure(isnan(ModPressure)) = 0;
Modelsurfacepressure = ModPressure(:,:,1);

Twisted_SSH_anom = (Modelsurfacepressure/g);      

%%% Write initial model surface height anomaly to filename.  Try making it
%%% zero....Initiaize with zero deviations from free surface
Modelsurfaceheight = zeros(size(TwistedSaltPretzal,1),size(TwistedSaltPretzal,2));








%%%%%%%%%%% TEMPERATURE %%%%%%%%%%%%%%%%%


ModelTwistedTheta = interp3(XSC,YSC,RSC,(TwistedTheta),Xm3C,Ym3C,RM3C,'linear');
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











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading files for Sea Ice Model%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Have decided today to make snow cover .3 throughout


SIThick = rdmds(fullfile(datadir,'SIheff'),'rec',recnum);
SIArea = rdmds(fullfile(datadir,'IceConc'),'rec',recnum);
SSS = rdmds(fullfile(datadir,'SSSdaily'),'rec',recnum);
SIUvel = rdmds(fullfile(datadir,'SIuice'),'rec',recnum);
SIVvel = rdmds(fullfile(datadir,'SIvice'),'rec',recnum);
SIsalt = rdmds(fullfile(datadir,'oceSflux'),'rec',recnum);


%Change longitude change of respective variable field
SIThick = [SIThick(idx1,:,:) ; SIThick(idx2,:,:)];
SIArea = [SIArea(idx1,:,:) ; SIArea(idx2,:,:)];
SSS = [SSS(idx1,:,:) ; SSS(idx2,:,:)];
SIUvel = [SIUvel(idx1,:,:) ; SIUvel(idx2,:,:)];
SIVvel = [SIVvel(idx1,:,:) ; SIVvel(idx2,:,:)];
SIsalt = [SIsalt(idx1,:,:) ; SIsalt(idx2,:,:)];


%%% Remove topography
SIThick(hFacC(:,:,1)==0) = NaN; 
SIArea(hFacC(:,:,1)==0) = NaN;
SSS(hFacC(:,:,1)==0) = NaN;
SIUvel(hFacC(:,:,1)==0) = NaN;
SIVvel(hFacC(:,:,1)==0) = NaN;
SIsalt(hFacC(:,:,1)==0) = NaN;


%%%%%%%%% Twisting Sea Ice //\\\\///\\\//\\//\\ %%%%%%%%%

Twisted_SIThick = Tweddell(SOSElonC,SOSElatC,SIThick);
Twisted_SIArea = Tweddell(SOSElonC,SOSElatC,SIArea);
Twisted_SSS = Tweddell(SOSElonC,SOSElatC,SSS);
Twisted_SIUvel = Tweddell(SOSElonC,SOSElatC,SIUvel);
Twisted_SIVvel = Tweddell(SOSElonC,SOSElatC,SIVvel);
Twisted_SIsalt = Tweddell(SOSElonC,SOSElatC,SIsalt);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating respective sea ice fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Twisted versions %%%%%%%%%%%%%%%%%%%%%%%%%%%


Twisted_MSIThick = interp2(XC',YC',(Twisted_SIThick)',XMC,YMC,'linear');
Twisted_MSIThick = inpaint_nans(Twisted_MSIThick,4);
Twisted_MSIThick = Twisted_MSIThick';



Twisted_MSIArea = interp2(XC',YC',(Twisted_SIArea)',XMC,YMC,'linear');
Twisted_MSIArea = inpaint_nans(Twisted_MSIArea,4);
Twisted_MSIArea = Twisted_MSIArea';

%SSS
Twisted_MSSS = interp2(XC',YC',(Twisted_SSS)',XMC,YMC,'linear');
Twisted_MSSS = inpaint_nans(Twisted_MSSS,4);
Twisted_MSSS = Twisted_MSSS';


Twisted_MSIUvel = interp2(XC',YC',Twisted_SIUvel',XMC,YMC,'linear');
Twisted_MSIUvel = inpaint_nans(Twisted_MSIUvel,4);
Twisted_MSIUvel = Twisted_MSIUvel';



Twisted_MSIVvel = interp2(XC',YC',Twisted_SIVvel',XMC,YMC,'linear');
Twisted_MSIVvel = inpaint_nans(Twisted_MSIVvel,4);
Twisted_MSIVvel = Twisted_MSIVvel';


Twisted_SIsalt = (Twisted_SIsalt)*86400;
Twisted_MSIsalt = interp2(XC',YC',Twisted_SIsalt',XMC,YMC,'linear');
Twisted_MSIsalt = inpaint_nans(Twisted_MSIsalt,4);
Twisted_MSIsalt = Twisted_MSIsalt'; 


%%%%%%%%%%% Prescribing %%%%%%%%%%%%%%
%%%%%%%%%%%.3m of snow everywhere  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Twisted_Snow = ones(size(Twisted_MSIThick,1),size(Twisted_MSIThick,2))*.3;


%%%%%%%%%%%%%%% Initializing SI salinity at 0 m/s to try


Twisted_MSIsalt = zeros(size(Twisted_MSIsalt,1),size(Twisted_MSIsalt,2));










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ice Shelf Load Anomaly 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



salt_ref = 34;
temp_ref = -1.9;
MIce_Anom = zeros(Nx,Ny);
for i=1:Nx
  for j=1:Ny

      MIce_Anom(i,j) = rho0*Modelsurfaceheight(i,j);

    for k=1:Nr
      if (zz(k) < icedraft(i,j))
          break
      end

      Pressure = -rho0*g*zz(k);
%       rhoShelfIce = densmdjwf(hydrogSalt(i,j),hydrogTheta(i,j),Pressure/Pa1dbar);
      rhoShelfIce = densmdjwf(salt_ref,temp_ref,Pressure/Pa1dbar);
      MIce_Anom(i,j) = MIce_Anom(i,j) + (g*(rhoShelfIce-rho0)*dz(k));                

     end
   end
end

data = MIce_Anom;
writeDataset(data,fullfile(inputconfigdir,SHELFICEloadAnomalyFile),ieee,prec);
clear data











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WRITING FILES TO INPUT PATH FOR INTITIALIZATION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = ModelTwistedTheta;
writeDataset(data,fullfile(inputconfigdir,hydrogThetaFile),ieee,prec);
clear data


data = TwistedSaltPretzal;
writeDataset(data,fullfile(inputconfigdir,hydrogSaltFile),ieee,prec);
clear data


data = MTwistUvel;
writeDataset(data,fullfile(inputconfigdir,uVelInitFile),ieee,prec);
clear data


data = MTwistVvel;
writeDataset(data,fullfile(inputconfigdir,vVelInitFile),ieee,prec);
clear data


data = Modelsurfaceheight;
writeDataset(data,fullfile(inputconfigdir,pSurfInitFile),ieee,prec);
clear data







%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%SEA ICE%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%




data = Twisted_MSIThick;
writeDataset(data,fullfile(inputconfigdir,HeffFile),ieee,prec);
clear data



data = Twisted_MSIArea;
writeDataset(data,fullfile(inputconfigdir,AreaFile),ieee,prec);
clear data



data = Twisted_MSIsalt;
writeDataset(data,fullfile(inputconfigdir,HsaltFile),ieee,prec);
clear data



data = Twisted_Snow;
writeDataset(data,fullfile(inputconfigdir,HsnowFile),ieee,prec);
clear data



data = Twisted_MSIUvel;
writeDataset(data,fullfile(inputconfigdir,uIceFile),ieee,prec);
clear data



data = Twisted_MSIVvel;
writeDataset(data,fullfile(inputconfigdir,vIceFile),ieee,prec);
clear data



