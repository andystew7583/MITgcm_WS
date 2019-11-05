%%%
%%% HigherRes.m
%%%
%%% Takes an MITgcm experiment and doubles its resolution, producing input
%%% files for a new experiment at the higher resolution.
%%%
%%% The resolution-doubling process used in this file is like the one for
%%% initial conditions.

%%%%%% Set experiment name in analysis folder correctly before doing this

%%% For file I/O
addpath ../newexp_utils/
addpath ../utils/matlab

%%% Load experiment
run ../analysis/setExpname
run ../analysis/loadexp
expiter = 1166400;       
Nx = 304;
Ny = 232;
% Nr = 150;
Nr = 139;

%%% Limit values of dependent variables
Theta_max = 2;
Theta_min = -2;
Salt_max = 37;
Salt_min = 33;
Uvel_max = 0.3;
Uvel_min = -0.3;
Vvel_max = 0.3;
Vvel_min = -0.3;


[YSC,XSC,RSC]=meshgrid(yy,xx,zz);

%%% Pull out u,v,t,s from pickup file
A = rdmds(fullfile(expdir,expname,'results/pickup'),expiter);
uvtse1 = A(:,:,[1:4*Nr 4*Nr+1]);

OldUvel = uvtse1(:,:,1:Nr);
OldVvel = uvtse1(:,:,Nr+1:2*Nr);
OldTheta = uvtse1(:,:,2*Nr+1:3*Nr);
OldSalt = uvtse1(:,:,3*Nr+1:4*Nr);
OldpSurf = uvtse1(:,:,4*Nr+1);

OldUvel(OldUvel==0)=NaN;
OldVvel(OldVvel==0)=NaN;
OldTheta(OldTheta==0)=NaN;
OldSalt(OldSalt==0)=NaN;


Seaice = rdmds(fullfile(expdir,expname,'/results/pickup_seaice'),expiter);


OldSItices = Seaice(:,:,1);
OldSIArea = Seaice(:,:,8);
OldSIHeff= Seaice(:,:,9);
OldSIsnow =Seaice(:,:,10);
OldSIsalt =Seaice(:,:,11);
OldSIUice = Seaice(:,:,12);
OldSIVice = Seaice(:,:,13);

OldSIArea(OldSIArea==0)=NaN;
OldSIHeff(OldSIHeff==0)=NaN;
OldSIsnow (OldSIsnow==0)=NaN;
OldSIsalt(OldSIsalt==0)=NaN;
OldSIUice(OldSIUice==0)=NaN;
OldSIVice(OldSIVice==0)=NaN;

%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%% Do interpolation similar to setinitial files



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ITERPOLATE DATA ONTO OUR GRID %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defineGrid
 
 
fid = fopen(fullfile(inputfolder,bathyFile),'r','b');
h = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);

fid = fopen(fullfile(inputfolder,SHELFICEtopoFile),'r','b');
icedraft = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);


%%%%%%%%%%% TEMPERATURE %%%%%%%%%%%%%%%%%

% ModelTwistedTheta = interp3(YSC,XSC,RSC,(OldTheta),Ym3C,Xm3C,RM3C,'nearest');           
      




ModelTwistedTheta = interp3(YSC,XSC,RSC,(OldTheta),Ym3C,Xm3C,RM3C,'linear');
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
% 
ModelTwistedTheta(ModelTwistedTheta>Theta_max) = Theta_max;
ModelTwistedTheta(ModelTwistedTheta<Theta_min) = Theta_min;

%%%%%%%%%%%%%%%%%%%  Salinity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TwistedSaltPretzal = NaN(Ny,Nx,Nr);
% salt_temp=zeros(size(OldSalt,1),size(OldSalt,2));
% for k = 1:size(OldSalt,3)
%       salt_temp = OldSalt(:,:,k);
%       TwistedSaltPretzal(:,:,k)=interp2(XC',YC',(salt_temp)',XMC,YMC,'nearest');           
%       TwistedSaltPretzal(:,:,k) = inpaint_nans(TwistedSaltPretzal(:,:,k),4);
% end



TwistedSaltPretzal = interp3(YSC,XSC,RSC,(OldSalt),Ym3C,Xm3C,RM3C,'linear'); 
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
% 
TwistedSaltPretzal(TwistedSaltPretzal>Salt_max) = Salt_max;
TwistedSaltPretzal(TwistedSaltPretzal<Salt_min) = Salt_min;


%%%%%%%%%%%%%%%% UVELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Please use XG,YC,RC as SOSE original grid (west point) 

%%%%%%%%%%%%%%%%%%% Twisted UVEL %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MTwistUvel = NaN(Ny,Nx,Nr);
% Uvel_temp=zeros(size(OldUvel,1),size(OldUvel,2));
% for k = 1:size(OldUvel,3)
%       Uvel_temp = OldUvel(:,:,k);
%       MTwistUvel(:,:,k)=interp2(XC',YC',(Uvel_temp)',XMC,YMC,'nearest');           
%       MTwistUvel(:,:,k) = inpaint_nans(MTwistUvel(:,:,k),4);
% end


MTwistUvel = interp3(YSC,XSC,RSC,(OldUvel),Ym3C,Xm3C,RM3C,'linear'); 
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

%%%%%%%%%%%%%%% VVELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PlZ use XC,YG,RC as original SOSE grid (on South Point) %%

%%%%%%%%%%%%%%%%%%% Twisted VVEL %%%%%%%%%%%%%%%%%%%%%%%%%%

% MTwistVvel = NaN(Ny,Nx,Nr);
% Vvel_temp=zeros(size(OldVvel,1),size(OldVvel,2));
% for k = 1:size(OldVvel,3)
%       Vvel_temp = OldVvel(:,:,k);
%       MTwistVvel(:,:,k)=interp2(XC',YC',(Vvel_temp)',XMC,YMC,'nearest');           
%       MTwistVvel(:,:,k) = inpaint_nans(MTwistVvel(:,:,k),4);
% end

MTwistVvel = interp3(YSC,XSC,RSC,(OldVvel),Ym3C,Xm3C,RM3C,'linear'); 
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

%%%%%%%%%%% VELOCITY VECTOR TWIST %%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pressure

% pressure_temp = PHIHYD;
% Twisted_Pressure = Tweddell(SOSElonC,SOSElatC,pressure_temp);
% Twisted_Pressure=transpose3D(Twisted_Pressure);
% 
% 
% ModPressure = interp3(XSC,YSC,RSC,(Twisted_Pressure),Xm3C,Ym3C,RM3C,'linear'); 
% for k = 1:size(ModPressure,3)
%   ModPressure(:,:,k)=inpaint_nans(ModPressure(:,:,k),4);
% end
% ModPressure(ModPressure==0) = NaN;
% for i=1:Nx
%   Pres_slice = inpaint_nans(squeeze(ModPressure(:,i,:)),4);
%   ModPressure(:,i,:) = reshape(Pres_slice,[Ny 1 Nz]); 
% end
% ModPressure(ModPressure==0) = NaN;
% for k = 1:size(ModPressure,3)
%   ModPressure(:,:,k)=inpaint_nans(ModPressure(:,:,k),4);
% end
% 
% ModPressure=transpose3D(ModPressure);
% % ModPressure(isnan(ModPressure)) = 0;
% Modelsurfacepressure = ModPressure(:,:,1);
% 
% Twisted_SSH_anom = (Modelsurfacepressure/g);      

%%% Write initial model surface height anomaly to filename.  Try making it
%%% zero....Initiaize with zero deviations from free surface
Modelsurfaceheight = zeros(size(ModelTwistedTheta,1),size(ModelTwistedTheta,2));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading files for Sea Ice Model%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating respective sea ice fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Twisted versions %%%%%%%%%%%%%%%%%%%%%%%%%%%
max_SIThick = 10;

Twisted_MSIThick = interp2(XC',YC',(OldSIHeff)',XMC,YMC,'linear');
Twisted_MSIThick = inpaint_nans(Twisted_MSIThick,4);
Twisted_MSIThick = Twisted_MSIThick';
Twisted_MSIThick(Twisted_MSIThick>max_SIThick)=max_SIThick;


Twisted_MSIArea = interp2(XC',YC',(OldSIArea)',XMC,YMC,'linear');
Twisted_MSIArea = inpaint_nans(Twisted_MSIArea,4);
Twisted_MSIArea = Twisted_MSIArea';



Twisted_MSIUvel = interp2(XC',YC',OldSIUice',XMC,YMC,'linear');
Twisted_MSIUvel = inpaint_nans(Twisted_MSIUvel,4);
Twisted_MSIUvel = Twisted_MSIUvel';



Twisted_MSIVvel = interp2(XC',YC',OldSIVice',XMC,YMC,'linear');
Twisted_MSIVvel = inpaint_nans(Twisted_MSIVvel,4);
Twisted_MSIVvel = Twisted_MSIVvel';



Twisted_MSIsalt = interp2(XC',YC',OldSIsalt',XMC,YMC,'linear');
Twisted_MSIsalt = inpaint_nans(Twisted_MSIsalt,4);
Twisted_MSIsalt = Twisted_MSIsalt'; 


%%%%%%%%%%% Prescribing %%%%%%%%%%%%%%
%%%%%%%%%%%.3m of snow everywhere  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Twisted_Snow = interp2(XC',YC',OldSIsnow',XMC,YMC,'linear');
Twisted_Snow = inpaint_nans(Twisted_Snow,4);
Twisted_Snow = Twisted_Snow'; 







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
Modelsurfacepressure=zeros(Nx,Ny,1);
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







%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%SEA ICE%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%




data = Twisted_MSIThick;
writeDataset(data,fullfile(inputfolder,HeffFile),ieee,prec);
clear data



data = Twisted_MSIArea;
writeDataset(data,fullfile(inputfolder,AreaFile),ieee,prec);
clear data



data = Twisted_MSIsalt;
writeDataset(data,fullfile(inputfolder,HsaltFile),ieee,prec);
clear data



data = Twisted_Snow;
writeDataset(data,fullfile(inputfolder,HsnowFile),ieee,prec);
clear data



data = Twisted_MSIUvel;
writeDataset(data,fullfile(inputfolder,uIceFile),ieee,prec);
clear data



data = Twisted_MSIVvel;
writeDataset(data,fullfile(inputfolder,vIceFile),ieee,prec);
clear data









