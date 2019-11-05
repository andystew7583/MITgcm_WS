%%%setting initial.bin files for mitgcm parameters for momemtum/tracers





%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%% SET-UP %%%%%%%%%%%%%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%  

defineGrid

        
%%% Matlab utilities 

addpath ../newexp_utils
addpath ../utils/matlab


%%%% setting min salinity 
min_salt_true =1;
min_salt = 34.15;

%%%%% Input file path for OBCS SOSE generated files

input = '/data1/MITgcm_WS/newexp/OBCS';


%%%%%% Changing XC,YC of SOSE to match OBCS files


idx1 = find(XC(:,1)>=180);
idx2 = find(XC(:,1)<180);
idx3 = find(XG(:,1)>=180);
idx4 = find(XG(:,1)<180);
XC = [XC(idx1,:)-360 ; XC(idx2,:)];
YC = [YC(idx1,:) ; YC(idx2,:)];
XG = [XG(idx3,:)-360 ; XG(idx4,:)];
YG = [YG(idx3,:) ; YG(idx4,:)];

% idx5 = find(XC(:,1)>-83 & XC(:,1) <-3);
idx5 = find(XC(:,1)>-83 & XC(:,1) <21);
idx6 = find(YC(1,:)>-84 & YC(1,:)<-64);

XC = XC(idx5,idx6);
YC = YC(idx5,idx6);

XG = XG(idx5,idx6);
YG = YG(idx5,idx6);

hFacC = [hFacC(idx1,:,:) ; hFacC(idx2,:,:)];
hFacW = [hFacW(idx1,:,:) ; hFacW(idx2,:,:)];
hFacS = [hFacS(idx1,:,:) ; hFacS(idx2,:,:)];

hFacC = hFacC(idx5,idx6,:); 
hFacW = hFacW(idx5,idx6,:);
hFacS = hFacS(idx5,idx6,:); 


% X_num = 738;
X_num = 624;
Y_num = 84;

%%% Load respective OBCS files

fid = fopen(fullfile(input,'SIThicktwenty.bin'),'r','b');
SIThick = fread(fid,[X_num*Y_num,216],'real*8'); 
fclose(fid);

SIThick = reshape(SIThick,X_num,Y_num,216);



fid = fopen(fullfile(input,'Areatwenty.bin'),'r','b');
SIArea = fread(fid,[X_num*Y_num,1080],'real*8'); 
fclose(fid);

SIArea = reshape(SIArea,X_num,Y_num,1080);



fid = fopen(fullfile(input,'uIcetwenty.bin'),'r','b');
SIUvel = fread(fid,[X_num*Y_num,360],'real*8'); 
fclose(fid);

SIUvel = reshape(SIUvel,X_num,Y_num,360);



fid = fopen(fullfile(input,'vIcetwenty.bin'),'r','b');
SIVvel = fread(fid,[X_num*Y_num,360],'real*8'); 
fclose(fid);

SIVvel = reshape(SIVvel,X_num,Y_num,360);



fid = fopen(fullfile(input,'Theta.bin'),'r','b');
Theta = fread(fid,[600*108,42*216],'real*8'); 
fclose(fid);

Theta = reshape(Theta,600,108,42,216);


fid = fopen(fullfile(input,'Theta_northtwenty.bin'),'r','b');
Theta_north = fread(fid,[X_num*1,42*216],'real*8'); 
fclose(fid);

Theta_north = reshape(Theta_north,X_num,1,42,216);


fid = fopen(fullfile(input,'Theta_easttwenty.bin'),'r','b');
Theta_east = fread(fid,[Y_num*1,42*216],'real*8'); 
fclose(fid);

Theta_east = reshape(Theta_east,1,Y_num,42,216);


fid = fopen(fullfile(input,'Theta_westtwenty.bin'),'r','b');
Theta_west = fread(fid,[Y_num*1,42*216],'real*8'); 
fclose(fid);

Theta_west = reshape(Theta_west,1,Y_num,42,216);



fid = fopen(fullfile(input,'Salt_northtwenty.bin'),'r','b');
Salt_north = fread(fid,[X_num*1,42*216],'real*8'); 
fclose(fid);

Salt_north = reshape(Salt_north,X_num,1,42,216);


fid = fopen(fullfile(input,'Salt_easttwenty.bin'),'r','b');
Salt_east = fread(fid,[Y_num*1,42*216],'real*8'); 
fclose(fid);

Salt_east = reshape(Salt_east,1,Y_num,42,216);



fid = fopen(fullfile(input,'Uvel_northtwenty.bin'),'r','b');
Uvel_north = fread(fid,[X_num*1,42*216],'real*8'); 
fclose(fid);

Uvel_north = reshape(Uvel_north,X_num,1,42,216);


fid = fopen(fullfile(input,'Uvel_easttwenty.bin'),'r','b');
Uvel_east = fread(fid,[Y_num*1,42*216],'real*8'); 
fclose(fid);

Uvel_east = reshape(Uvel_east,1,Y_num,42,216);




fid = fopen(fullfile(input,'Vvel_northtwenty.bin'),'r','b');
Vvel_north = fread(fid,[X_num*1,42*216],'real*8'); 
fclose(fid);

Vvel_north = reshape(Vvel_north,X_num,1,42,216);


fid = fopen(fullfile(input,'Vvel_easttwenty.bin'),'r','b');
Vvel_east = fread(fid,[Y_num*1,42*216],'real*8'); 
fclose(fid);

Vvel_east = reshape(Vvel_east,Y_num,1,42,216);



fid = fopen(fullfile(input,'SIThick_northtwenty.bin'),'r','b');
SIThick_north = fread(fid,[X_num*1,216],'real*8'); 
fclose(fid);

SIThick_north = reshape(SIThick_north,X_num,1,216);

fid = fopen(fullfile(input,'SIThick_easttwenty.bin'),'r','b');
SIThick_east = fread(fid,[Y_num*1,216],'real*8'); 
fclose(fid);

SIThick_east = reshape(SIThick_east,1,Y_num,216);



fid = fopen(fullfile(input,'Area_northtwenty.bin'),'r','b');
SIArea_north = fread(fid,[X_num*1,1080],'real*8'); 
fclose(fid);

SIArea_north = reshape(SIArea_north,X_num,1,1080);

fid = fopen(fullfile(input,'Area_easttwenty.bin'),'r','b');
SIArea_east = fread(fid,[Y_num*1,1080],'real*8'); 
fclose(fid);

SIArea_east = reshape(SIArea_east,1,Y_num,1080);



fid = fopen(fullfile(input,'uIce_northtwenty.bin'),'r','b');
SIUvel_north = fread(fid,[X_num*1,360],'real*8'); 
fclose(fid);

SIUvel_north = reshape(SIUvel_north,X_num,1,360);

fid = fopen(fullfile(input,'uIce_easttwenty.bin'),'r','b');
SIUvel_east = fread(fid,[Y_num*1,360],'real*8'); 
fclose(fid);

SIUvel_east = reshape(SIUvel_east,1,Y_num,360);




fid = fopen(fullfile(input,'vIce_northtwenty.bin'),'r','b');
SIVvel_north = fread(fid,[X_num*1,360],'real*8'); 
fclose(fid);

SIVvel_north = reshape(SIVvel_north,X_num,1,360);

fid = fopen(fullfile(input,'vIce_easttwenty.bin'),'r','b');
SIVvel_east = fread(fid,[Y_num*1,360],'real*8'); 
fclose(fid);

SIVvel_east = reshape(SIVvel_east,1,Y_num,360);


fid = fopen(fullfile(input,'vIce_westtwenty.bin'),'r','b');
SIVvel_west = fread(fid,[1*Y_num,360],'real*8'); 
fclose(fid);

SIVvel_west = reshape(SIVvel_west,1,Y_num,360);



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

clear hydrogSalt
clear hydrogTheta


icedraft_east = icedraft(end,:);

bathy_east = h(end,:);

clear h icedraft





%%%%%%%%% SOSE POINTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOSElonC = XC(:,1);
SOSElatC= YC(1,:);


%%%%%% for velocity points ----> 
SOSElonG= XG(:,1);
SOSElatG= YG(1,:);

Theta_north = squeeze(Theta_north);
Salt_north = squeeze(Salt_north);
Uvel_north = squeeze(Uvel_north);
Vvel_north = squeeze(Vvel_north);
SIThick_north = squeeze(SIThick_north);
SIArea_north = squeeze(SIArea_north);
SIUvel_north = squeeze(SIUvel_north);
SIVvel_north = squeeze(SIVvel_north);







% % % % % % % % % meshgrids for 2D interp
[XS,YS]=meshgrid(SOSElonC,SOSElatC);
[WW,KK]=meshgrid(xmc,ymc);

[XX,ZZ] = meshgrid(SOSElonC,RC);
[YY,ZY] = meshgrid(SOSElatC,RC);



%%%%%%%%%%% OUR GRID %%%%%%%%%%%

[XM,ZM] = meshgrid(xmc,mrc);
[YM,GR] = meshgrid(ymc,mrc);


NumRec = 216;
%%%%Getting total number of months
years = NumRec/72;
months = NumRec/18;
DPY = 365;
DPM = 6;
YDPM = 18;
DM = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Northern Boundary Interp %%%

%%%%%%%%% Temperature

OBNt = NaN(Nx,Nz,NumRec);
Theta_north(Theta_north==0)=NaN;
TN_temp = zeros(size(Theta_north,1),size(Theta_north,2));
for i = 1:size(Theta_north,3)
    TN_temp = Theta_north(:,:,i);
    OBNt(:,:,i) = interp2(XX,ZZ,TN_temp',XM',ZM');
    OBNt(:,:,i) = inpaint_nans(OBNt(:,:,i),4);
end

OBNt = reshape(OBNt,size(OBNt,1),size(OBNt,2),DPM,months,years);
OBNt = squeeze(nanmean(nanmean(OBNt,3),5));

OBNs = NaN(Nx,Nz,NumRec);
Salt_north(Salt_north==0)=NaN;
SN_temp = zeros(size(Salt_north,1),size(Salt_north,2));
for i = 1:size(Salt_north,3)
    SN_temp = Salt_north(:,:,i);
    OBNs(:,:,i) = interp2(XX,ZZ,SN_temp',XM',ZM');
    OBNs(:,:,i) = inpaint_nans(OBNs(:,:,i),4);
end

OBNs = reshape(OBNs,size(OBNs,1),size(OBNs,2),DPM,months,years);
OBNs = squeeze(nanmean(nanmean(OBNs,3),5));


if min_salt_true ==1

    OBNs(OBNs<min_salt) = min_salt;
end



OBNu = NaN(Nx,Nz,NumRec);
Uvel_north(Uvel_north==0)=NaN;
UN_temp = zeros(size(Uvel_north,1),size(Uvel_north,2));
for i = 1:size(Uvel_north,3)
    UN_temp = Uvel_north(:,:,i);
    OBNu(:,:,i) = interp2(XX,ZZ,UN_temp',XM',ZM');
    OBNu(:,:,i) = inpaint_nans(OBNu(:,:,i),4);
end

OBNu = reshape(OBNu,size(OBNu,1),size(OBNu,2),DPM,months,years);
OBNu = squeeze(nanmean(nanmean(OBNu,3),5));

OBNv = NaN(Nx,Nz,NumRec);
Vvel_north(Vvel_north==0)=NaN;
VN_temp = zeros(size(Vvel_north,1),size(Vvel_north,2));
for i = 1:size(Vvel_north,3)
    VN_temp = Vvel_north(:,:,i);
    OBNv(:,:,i) = interp2(XX,ZZ,VN_temp',XM',ZM');
    OBNv(:,:,i) = inpaint_nans(OBNv(:,:,i),4);
end

OBNv = reshape(OBNv,size(OBNv,1),size(OBNv,2),DPM,months,years);
OBNv = squeeze(nanmean(nanmean(OBNv,3),5));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Eastern Boundary Interp %%%

%%%%%%%%%
%%%%% Try 2 - Use isobath to delineate depth of stratification
OBEt = NaN(Ny,Nz,NumRec);
Theta_east = squeeze(Theta_east);
Theta_east(Theta_east==0)=NaN;
TE_temp = zeros(size(Theta_east,1),size(Theta_east,2));
for i = 1:size(Theta_east,3)
    TE_temp = Theta_east(:,:,i);
    OBEt(:,:,i) = interp2(YY,ZY,TE_temp',YM',GR');
 
        for j=1:Ny
            for k = 1:Nr
                 Pressure = -rho0*g*zz(1);

                 if bathy_east(j)>=-400
                     OBEt(j,k,i) = .0901 - .0575*salt_e(j,1) - (7.61e-4 *(Pressure/Pa1dbar));

                 else

                      OBEt(j,k,i) = OBEt(j,k,i);
                     
                 end
            end
        end
    OBEt(:,:,i) = inpaint_nans(OBEt(:,:,i),4);
end

OBEt = reshape(OBEt,size(OBEt,1),size(OBEt,2),DPM,months,years);
OBEt = squeeze(nanmean(nanmean(OBEt,3),5));



OBEs = NaN(Ny,Nz,NumRec);
Salt_east = squeeze(Salt_east);
hFacC_e = squeeze(hFacC(end,:,:));
Salt_east(Salt_east==0)=NaN;
SE_temp = zeros(size(Salt_east,1),size(Salt_east,2));
for i = 1:size(Salt_east,3)
    SE_temp = Salt_east(:,:,i);
    OBEs(:,:,i) = interp2(YY,ZY,SE_temp',YM',GR');
        for j=1:Ny

           for k = 1:Nr
                 if  bathy_east(j) >= -400
                        OBEs(j,k,i) = 34.5;
                     else
                        OBEs(j,k,i) = OBEs(j,k,i);
                      
                 end
           end

        end
    OBEs(:,:,i) = inpaint_nans(OBEs(:,:,i),4);
%     
    
end

OBEs = reshape(OBEs,size(OBEs,1),size(OBEs,2),DPM,months,years);
OBEs = squeeze(nanmean(nanmean(OBEs,3),5));


if min_salt_true ==1

    OBEs(OBEs<min_salt) = min_salt;
end



%%%%% Try 4: set temperature-dependent salinity at the
%%%%% boundaries...salinities of 34, temperatures will be sfc freezing
%%%%% temperature

    for j = 1:Ny
        for k = 1:Nr
            for i = 1:12
                if OBEs(j,k,i)==min_salt
                    OBEt(j,k,i) = .0901 - .0575*34.15 - (7.61e-4 *(Pressure/Pa1dbar));
                end
            end
        end
    end
    

    for j = 1:Nx
        for k = 1:Nr
            for i = 1:12
                if OBNs(j,k,i)==min_salt
                    OBNt(j,k,i) = .0901 - .0575*34.15 - (7.61e-4 *(Pressure/Pa1dbar));
                end
            end
        end
    end



OBEu = NaN(Ny,Nz,NumRec);
Uvel_east = squeeze(Uvel_east);
Uvel_east(Uvel_east==0)=NaN;
UE_temp = zeros(size(Uvel_east,1),size(Uvel_east,2));
for i = 1:size(Uvel_east,3)
    UE_temp = Uvel_east(:,:,i);
    OBEu(:,:,i) = interp2(YY,ZY,UE_temp',YM',GR');
    OBEu(:,:,i) = inpaint_nans(OBEu(:,:,i),4);
end

OBEu = reshape(OBEu,size(OBEu,1),size(OBEu,2),DPM,months,years);
OBEu = squeeze(nanmean(nanmean(OBEu,3),5));


OBEv = NaN(Ny,Nz,NumRec);
Vvel_east = squeeze(Vvel_east);
Vvel_east(Vvel_east==0)=NaN;
VE_temp = zeros(size(Vvel_east,1),size(Vvel_east,2));
for i = 1:size(Vvel_east,3)
    VE_temp = Vvel_east(:,:,i);
    OBEv(:,:,i) = interp2(YY,ZY,VE_temp',YM',GR');
    OBEv(:,:,i) = inpaint_nans(OBEv(:,:,i),4);
end

OBEv = reshape(OBEv,size(OBEv,1),size(OBEv,2),DPM,months,years);
OBEv = squeeze(nanmean(nanmean(OBEv,3),5));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Western Boundary Interp %%%

%%%%%%%%% Temperature

% OBWt = NaN(Ny,Nz,NumRec);
% Theta_west = squeeze(Theta_west);
% Theta_west(Theta_west==0)=NaN;
% TW_temp = zeros(size(Theta_west,1),size(Theta_west,2));
% for i = 1:size(Theta_west,3)
%     TW_temp = Theta_west(:,:,i);
%     OBWt(:,:,i) = interp2(YY,ZY,TW_temp',YM',GR');
%     OBWt(:,:,i) = inpaint_nans(OBWt(:,:,i),4);
% end
% 
% 
% OBWt = reshape(OBWt,size(OBWt,1),size(OBWt,2),DPM,months,years);
% OBWt = squeeze(nanmean(nanmean(OBWt,3),5));
% 
% 
% OBWs = NaN(Ny,Nz,NumRec);
% Salt_west = squeeze(Salt_west);
% Salt_west(Salt_west==0)=NaN;
% SW_temp = zeros(size(Salt_west,1),size(Salt_west,2));
% for i = 1:size(Salt_west,3)
%     SW_temp = Salt_west(:,:,i);
%     OBWs(:,:,i) = interp2(YY,ZY,SW_temp',YM',GR');
%     OBWs(:,:,i) = inpaint_nans(OBWs(:,:,i),4);
% end
% 
% OBWs = reshape(OBWs,size(OBWs,1),size(OBWs,2),DPM,months,years);
% OBWs = squeeze(nanmean(nanmean(OBWs,3),5));
% 
% 
% OBWu = NaN(Ny,Nz,NumRec);
% Uvel_west = squeeze(Uvel_west);
% Uvel_west(Uvel_west==0)=NaN;
% UW_temp = zeros(size(Uvel_west,1),size(Uvel_west,2));
% for i = 1:size(Uvel_west,3)
%     UW_temp = Uvel_west(:,:,i);
%     OBWu(:,:,i) = interp2(YY,ZY,UW_temp',YM',GR');
%     OBWu(:,:,i) = inpaint_nans(OBWu(:,:,i),4);
% end
% 
% OBWu = reshape(OBWu,size(OBWu,1),size(OBWu,2),DPM,months,years);
% OBWu = squeeze(nanmean(nanmean(OBWu,3),5));
% 
% OBWv = NaN(Ny,Nz,NumRec);
% Vvel_west = squeeze(Vvel_west);
% Vvel_west(Vvel_west==0)=NaN;
% VW_temp = zeros(size(Vvel_west,1),size(Vvel_west,2));
% for i = 1:size(Vvel_west,3)
%     VW_temp = Vvel_west(:,:,i);
%     OBWv(:,:,i) = interp2(YY,ZY,VW_temp',YM',GR');
%     OBWv(:,:,i) = inpaint_nans(OBWv(:,:,i),4);
% end
% 
% OBWv = reshape(OBWv,size(OBWv,1),size(OBWv,2),DPM,months,years);
% OBWv = squeeze(nanmean(nanmean(OBWv,3),5));
% 
% 
% OBWh = NaN(Ny,NumRec);
% SIThick_west = squeeze(SIThick_west);
% % ThickW_Temp = zeros(size(SIThick_west,1));
% for i = 1:size(SIThick_west,2);
%     ThickW_Temp = SIThick_west(:,i);
%     OBWh(:,i) = interp1(SOSElatC,ThickW_Temp,ymc);
%     OBWh(:,i) = inpaint_nans(OBWh(:,i),4);
% end
% 
% OBWh = reshape(OBWh,size(OBWh,1),DPM,months,years);
% OBWh = squeeze(nanmean(nanmean(OBWh,2),4));
% 
% 
% OBWa = NaN(Ny,NumRec);
% SIArea_west = squeeze(SIArea_west);
% % AreaW_Temp = zeros(size(SIArea_west,1));
% for i = 1:size(SIArea_west,2)
%     AreaW_Temp = SIArea_west(:,i);
%     OBWa(:,i) = interp1(SOSElatC,AreaW_Temp,ymc);
%     OBWa(:,i) = inpaint_nans(OBWa(:,i),4);
% end
% 
% OBWa = reshape(OBWa,size(OBWa,1),YDPM,months,5);
% OBWa = squeeze(nanmean(nanmean(OBWa,2),4));
% 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Alternate SEA ICE WRITING%%%%%%%%%%

SeaThick = NaN(Nx,Ny,NumRec);
SITemp = zeros(size(SIThick,1),size(SIThick,2));
for i = 1:size(SIThick,3)
    SITemp = SIThick(:,:,i);
    SeaThick(:,:,i) = interp2(XS,YS,SITemp',WW',KK');
    SeaThick(:,:,i) = inpaint_nans(SeaThick(:,:,i),4);
end

SIEast = squeeze(SeaThick(end,:,:));
OBEh = reshape(SIEast,size(SIEast,1),DPM,months,years);
OBEh = squeeze(nanmean(nanmean(OBEh,2),4));

SINorth = squeeze(SeaThick(:,end,:));
OBNh = reshape(SINorth,size(SINorth,1),DPM,months,years);
OBNh = squeeze(nanmean(nanmean(OBNh,2),4));





SeaArea = NaN(Nx,Ny,NumRec);
AreaTemp = zeros(size(SIArea,1),size(SIArea,2));
for i = 1:size(SIArea,3)
    AreaTemp = SIArea(:,:,i);
    SeaArea(:,:,i) = interp2(XS,YS,AreaTemp',WW',KK');
    SeaArea(:,:,i) = inpaint_nans(SeaArea(:,:,i),4);
end

AreaEast = squeeze(SeaArea(end,:,:));
OBEa = reshape(AreaEast,size(AreaEast,1),DM,months,years);
OBEa = squeeze(nanmean(nanmean(OBEa,2),4));

AreaNorth = squeeze(SeaArea(:,end,:));
OBNa = reshape(AreaNorth,size(AreaNorth,1),DM,months,years);
OBNa = squeeze(nanmean(nanmean(OBNa,2),4));

years = 1;

iceu = NaN(Nx,Ny,NumRec);
iceu_temp = zeros(size(SIUvel,1),size(SIUvel,2));
for i = 1:size(SIUvel,3)
    iceu_temp = SIUvel(:,:,i);
    iceu(:,:,i) = interp2(XS,YS,iceu_temp',WW',KK');
    iceu(:,:,i) = inpaint_nans(iceu(:,:,i),4);
end

SI_UVELEast = squeeze(iceu(end,:,:));
OBEuice = reshape(SI_UVELEast,size(SI_UVELEast,1),DM,months,years);
OBEuice = squeeze(nanmean(nanmean(OBEuice,2),4));

SI_UVELNORTH = squeeze(iceu(:,end,:));
OBNuice = reshape(SI_UVELNORTH,size(SI_UVELNORTH,1),DM,months,years);
OBNuice = squeeze(nanmean(nanmean(OBNuice,2),4));


% OBWuice = NaN(Ny,NumRec);
% SIUvel_west = squeeze(SIUvel_west);
% % SIUvelW_Temp = zeros(size(SIUvel_west,1));
% for i = 1:size(SIUvel_west,2);
%     SIUvelW_Temp = SIUvel_west(:,i);
%     OBWuice(:,i) = interp1(SOSElatC,SIUvelW_Temp,ymc);
%     OBWuice(:,i) = inpaint_nans(OBWuice(:,i),4);
% end
% 
% OBWuice = reshape(OBWuice,size(OBWuice,1),DM,months,years);
% OBWuice = squeeze(nanmean(nanmean(OBWuice,2),4));
% 
% 
% OBWvice = NaN(Ny,NumRec);
% SIVvel_west = squeeze(SIVvel_west);
% % SIVvelW_Temp = zeros(size(SIVvel_west,1));
% for i = 1:size(SIVvel_west,2);
%     SIVvelW_Temp = SIVvel_west(:,i);
%     OBWvice(:,i) = interp1(SOSElatC,SIVvelW_Temp,ymc);
%     OBWvice(:,i) = inpaint_nans(OBWvice(:,i),4);
% end
% 
% OBWvice = reshape(OBWvice,size(OBWvice,1),DM,months,years);
% OBWvice = squeeze(nanmean(nanmean(OBWvice,2),4));

icev = NaN(Nx,Ny,NumRec);
icev_temp = zeros(size(SIVvel,1),size(SIVvel,2));
for i = 1:size(SIVvel,3)
    icev_temp = SIVvel(:,:,i);
    icev(:,:,i) = interp2(XS,YS,icev_temp',WW',KK');
    icev(:,:,i) = inpaint_nans(icev(:,:,i),4);
end

SI_VVELEast = squeeze(icev(end,:,:));
OBEvice = reshape(SI_VVELEast,size(SI_VVELEast,1),DM,months,years);
OBEvice = squeeze(nanmean(nanmean(OBEvice,2),4));

SI_VVELNORTH = squeeze(icev(:,end,:));
OBNvice = reshape(SI_VVELNORTH,size(SI_VVELNORTH,1),DM,months,years);
OBNvice = squeeze(nanmean(nanmean(OBNvice,2),4));


OBNsl = zeros(size(OBNvice,1),size(OBNvice,2),size(OBNvice,3));
OBNsn = .3*ones(size(OBNvice,1),size(OBNvice,2),size(OBNvice,3));


OBEsl = zeros(size(OBEvice,1),size(OBEvice,2));
OBEsn = .3*ones(size(OBEvice,1),size(OBEvice,2));

% 
% OBWsl = zeros(size(OBWvice,1),size(OBWvice,2),size(OBWvice,3));
% OBWsn = .3*ones(size(OBWvice,1),size(OBWvice,2),size(OBWvice,3));
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    WRITING FILES %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




data = OBNt;
writeDataset(data,fullfile(inputfolder,OBNtFile),ieee,prec);
clear data



data = OBEt;
writeDataset(data,fullfile(inputfolder,OBEtFile),ieee,prec);
clear data


data = OBNs;
writeDataset(data,fullfile(inputfolder,OBNsFile),ieee,prec);
clear data




data = OBEs;
writeDataset(data,fullfile(inputfolder,OBEsFile),ieee,prec);
clear data



data = OBNu;
writeDataset(data,fullfile(inputfolder,OBNuFile),ieee,prec);
clear data



data = OBEu;
writeDataset(data,fullfile(inputfolder,OBEuFile),ieee,prec);
clear data


data = OBNv;
writeDataset(data,fullfile(inputfolder,OBNvFile),ieee,prec);
clear data



data = OBEv;
writeDataset(data,fullfile(inputfolder,OBEvFile),ieee,prec);
clear data




data = OBNa;
writeDataset(data,fullfile(inputfolder,OBNaFile),ieee,prec);
clear data


data = OBEa;
writeDataset(data,fullfile(inputfolder,OBEaFile),ieee,prec);
clear data


data = OBNh;
writeDataset(data,fullfile(inputfolder,OBNhFile),ieee,prec);
clear data



data = OBEh;
writeDataset(data,fullfile(inputfolder,OBEhFile),ieee,prec);
clear data


data = OBNsn;
writeDataset(data,fullfile(inputfolder,OBNsnFile),ieee,prec);
clear data




data = OBEsn;
writeDataset(data,fullfile(inputfolder,OBEsnFile),ieee,prec);
clear data


data = OBNuice;
writeDataset(data,fullfile(inputfolder,OBNuiceFile),ieee,prec);
clear data


data = OBEuice;
writeDataset(data,fullfile(inputfolder,OBEuiceFile),ieee,prec);
clear data


data = OBNvice;
writeDataset(data,fullfile(inputfolder,OBNviceFile),ieee,prec);
clear data



data = OBEvice;
writeDataset(data,fullfile(inputfolder,OBEviceFile),ieee,prec);
clear data


data = OBNsl;
writeDataset(data,fullfile(inputfolder,OBNslFile),ieee,prec);
clear data



data = OBEsl;
writeDataset(data,fullfile(inputfolder,OBEslFile),ieee,prec);
clear data


