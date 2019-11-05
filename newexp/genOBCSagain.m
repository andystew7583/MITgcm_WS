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

datadir = '/Volumes/JHAZE/MITgcm_WS/data/SOSEdata/08-10';

input = '/Volumes/JHAZE/MITgcm_WS/newexp/OBCS';

% Switch to [-180,180] longitude range for SOSE Data

idx1 = find(XC(:,1)>=180);
idx2 = find(XC(:,1)<180);
idx3 = find(XG(:,1)>=180);
idx4 = find(XG(:,1)<180);
XC = [XC(idx1,:)-360 ; XC(idx2,:)];
YC = [YC(idx1,:) ; YC(idx2,:)];
XG = [XG(idx3,:)-360 ; XG(idx4,:)];
YG = [YG(idx3,:) ; YG(idx4,:)];

idx5 = find(XC(:,1)>-83 & XC(:,1) <40);
idx6 = find(YC(1,:)>-84 & YC(1,:)<-64);

XC = XC(idx5,idx6);
YC = YC(idx5,idx6);

XG = XG(idx5,idx6);
YG = YG(idx5,idx6);

%%%%%% Change lat/lon of partial cell fields

hFacC = [hFacC(idx1,:,:) ; hFacC(idx2,:,:)];
hFacW = [hFacW(idx1,:,:) ; hFacW(idx2,:,:)];
hFacS = [hFacS(idx1,:,:) ; hFacS(idx2,:,:)];

hFacC = hFacC(idx5,idx6,:);
hFacW = hFacW(idx5,idx6,:);
hFacS = hFacS(idx5,idx6,:);


Nrec = 360;
%%% Matlab utilities 

addpath ../newexp_utils
addpath ../utils/matlab



SX = size(XC,1);
SY = size(XC,2);

SOSElonC = XC(:,1);
SOSElatC= XC(1,:);

% X_num = 624;
X_num = 738;
Y_num = 84;



%loop through total number of records to get ocean state


% theta_obcs_north=zeros(SX,1,42,Nrec);
% theta_obcs_east = zeros(1,SY,42,Nrec);
% theta_obcs_west = zeros(1,SY,42,Nrec);
% 
% 
% salt_obcs_north=zeros(SX,1,42,Nrec);
% salt_obcs_east = zeros(1,SY,42,Nrec);
% salt_obcs_west = zeros(1,SY,42,Nrec);
% 
% uvel_obcs_north=zeros(SX,1,42,Nrec);
% uvel_obcs_east = zeros(1,SY,42,Nrec);
% uvel_obcs_west = zeros(1,SY,42,Nrec);
% 
% vvel_obcs_north=zeros(SX,1,42,Nrec);
% vvel_obcs_east = zeros(1,SY,42,Nrec);
% vvel_obcs_west = zeros(1,SY,42,Nrec);
% 
% PHIHYD_obcs_north=zeros(SX,1,42,Nrec);
% PHIHYD_obcs_east = zeros(1,SY,42,Nrec);
% PHIHYD_obcs_west = zeros(1,SY,42,Nrec);
% % 
% 
% SIThick_obcs = NaN(SX,SY,Nrec);
% SIArea_obcs = NaN(SX,SY,Nrec);
SIUvel_obcs = NaN(SX,SY,Nrec);
SIVvel_obcs = NaN(SX,SY,Nrec);

% SIThick_obcs_north= zeros(SX,1,Nrec);
% SIThick_obcs_east = zeros(1,SY,Nrec);
% SIThick_obcs_west = zeros(1,SY,Nrec);

% SIArea_obcs_north = zeros(SX,1,Nrec);
% SIArea_obcs_east = zeros(1,SY,Nrec);
% SIArea_obcs_west = zeros(1,SY,Nrec);

% SIUvel_obcs_north=zeros(SX,1,Nrec);
% SIUvel_obcs_east = zeros(1,SY,Nrec);
% SIUvel_obcs_west = zeros(1,SY,Nrec);
% 
% SIVvel_obcs_north=zeros(SX,1,Nrec);
% SIVvel_obcs_east = zeros(1,SY,Nrec);
% SIVvel_obcs_west = zeros(1,SY,Nrec);

for i = 1:(Nrec)
%     Theta = rdmds(fullfile(datadir,'THETA'),'rec',i);
%     Salt = rdmds(fullfile(datadir,'SALT'),'rec',i);
%     Uvel = rdmds(fullfile(datadir,'UVEL'),'rec',i);
%     Vvel = rdmds(fullfile(datadir,'VVEL'),'rec',i);
%     PHIHYD = rdmds(fullfile(datadir,'PHIHYD'),'rec',i);
    
    
    
%Change longitude change of respective variable field

%     Theta = [Theta(idx1,:,:,:) ; Theta(idx2,:,:,:)];
%     Salt = [Salt(idx1,:,:,:) ; Salt(idx2,:,:,:)];
%     Uvel = [Uvel(idx1,:,:,:) ; Uvel(idx2,:,:,:)];
%     Vvel = [Vvel(idx1,:,:,:) ; Vvel(idx2,:,:,:)];
%     PHIHYD = [PHIHYD(idx1,:,:,:) ; PHIHYD(idx2,:,:,:)];
    

    
    
%     SIThick = rdmds(fullfile(datadir,'SIheff'),'rec',i);
%     SIArea = rdmds(fullfile(datadir,'IceConc'),'rec',i);
    SIUvel = rdmds(fullfile(datadir,'SIuice'),'rec',i);
    SIVvel = rdmds(fullfile(datadir,'SIvice'),'rec',i);

    
%     SIThick = [SIThick(idx1,:,:) ; SIThick(idx2,:,:)];
%     SIArea = [SIArea(idx1,:,:) ; SIArea(idx2,:,:)];
    SIUvel = [SIUvel(idx1,:,:) ; SIUvel(idx2,:,:)];
    SIVvel = [SIVvel(idx1,:,:) ; SIVvel(idx2,:,:)];
    

    
    
%     Theta = Theta(idx5,idx6,:,:); 
%     Theta_north = Theta(:,end,:,:);
%     Theta_east = Theta(end,:,:,:);
%     Theta_west = Theta(1,:,:,:);
%     
%     Salt = Salt(idx5,idx6,:,:);
%     Salt_north = Salt(:,end,:,:);
%     Salt_east = Salt(end,:,:,:);
%     Salt_west = Salt(1,:,:,:);
%     
%     Uvel = Uvel(idx5,idx6,:,:);
%     Uvel_north = Uvel(:,end,:,:);
%     Uvel_east = Uvel(end,:,:,:);
%     Uvel_west = Uvel(1,:,:,:);
%     
%     
%     Vvel = Vvel(idx5,idx6,:,:);
%     Vvel_north = Vvel(:,end,:,:);
%     Vvel_east = Vvel(end,:,:,:);
%     Vvel_west = Vvel(1,:,:,:);
%     
%     
%     PHIHYD = PHIHYD(idx5,idx6,:,:);
%     PHIHYD_north = PHIHYD(:,end,:,:);
%     PHIHYD_east = PHIHYD(end,:,:,:);
%     PHIHYD_west = PHIHYD(1,:,:,:);
%     
%     
%     theta_obcs_north(:,:,:,i)= Theta_north;
%     theta_obcs_east(:,:,:,i) = Theta_east;
%     theta_obcs_west(:,:,:,i) = Theta_west;
%     
%     salt_obcs_north(:,:,:,i)= Salt_north;
%     salt_obcs_east(:,:,:,i) = Salt_east;
%     salt_obcs_west(:,:,:,i) = Salt_west;
%     
%     uvel_obcs_north(:,:,:,i)= Uvel_north;
%     uvel_obcs_east(:,:,:,i) = Uvel_east;
%     uvel_obcs_west(:,:,:,i) = Uvel_west;
%     
%     vvel_obcs_north(:,:,:,i)= Vvel_north;
%     vvel_obcs_east(:,:,:,i) = Vvel_east;
%     vvel_obcs_west(:,:,:,i) = Vvel_west;
%     
%     PHIHYD_obcs_north(:,:,:,i) = PHIHYD_north;
%     PHIHYD_obcs_east(:,:,:,i) = PHIHYD_east;
%     PHIHYD_obcs_west(:,:,:,i) = PHIHYD_west;
    
 
    

    
    

%    Change longitude change of respective variable field
%     SIThick = SIThick(idx5,idx6,:);
%     SIThick_north = SIThick(:,end,:);
%     SIThick_east = SIThick(end,:,:);
%     SIThick_west = SIThick(1,:,:);
%     
%     SIArea = SIArea(idx5,idx6,:);
%     SIArea_north = SIArea(:,end,:);
%     SIArea_east = SIArea(end,:,:);
%     SIArea_west = SIArea(1,:,:);
%     
    
    SIUvel = SIUvel(idx5,idx6,:);
    SIUvel_north = SIUvel(:,end,:);
    SIUvel_east = SIUvel(end,:,:);
    SIUvel_west = SIUvel(1,:,:);
   
    
    SIVvel = SIVvel(idx5,idx6,:);
    SIVvel_north = SIVvel(:,end,:);
    SIVvel_east = SIVvel(end,:,:);
    SIVvel_west = SIVvel(1,:,:);
    
%     SIThick_obcs(:,:,i) = SIThick;
 
%     SIArea_obcs(:,:,i) = SIArea;
    
    SIUvel_obcs(:,:,i) = SIUvel;
   
    SIVvel_obcs(:,:,i) = SIVvel;
%     
%     SIThick_obcs_north(:,:,i) = SIThick_north;
%     SIThick_obcs_west(:,:,i) = SIThick_west;
%     SIThick_obcs_east(:,:,i) = SIThick_east;

%     SIArea_obcs_north(:,:,i) = SIArea_north;
%     SIArea_obcs_west(:,:,i) = SIArea_west;
%     SIArea_obcs_east(:,:,i) = SIArea_east;
    
%     SIUvel_obcs_north(:,:,i) = SIUvel_north;
%     SIUvel_obcs_west(:,:,i) = SIUvel_west;
%     SIUvel_obcs_east(:,:,i) = SIUvel_east;
%    
%     SIVvel_obcs_north(:,:,i) = SIVvel_north;
%     SIVvel_obcs_west(:,:,i) = SIVvel_west;
%     SIVvel_obcs_east(:,:,i) = SIVvel_east;
    
 
end

%%%%%%%%%% Reshaping to write in a binary file b/c too large for .mat
%%%%%%%%%% export

% theta_obcs_north = reshape(theta_obcs_north,X_num*1,42*216);
% theta_obcs_east = reshape(theta_obcs_east,Y_num*1,42*216);
% theta_obcs_west = reshape(theta_obcs_west,Y_num*1,42*216);
% 
% salt_obcs_north = reshape(salt_obcs_north,X_num*1,42*216);
% salt_obcs_east = reshape(salt_obcs_east,Y_num*1,42*216);
% salt_obcs_west = reshape(salt_obcs_west,Y_num*1,42*216);
% 
% Uvel_obcs_north = reshape(uvel_obcs_north,X_num*1,42*216);
% Uvel_obcs_east = reshape(uvel_obcs_east,Y_num*1,42*216);
% Uvel_obcs_west = reshape(uvel_obcs_west,Y_num*1,42*216);
% 
% Vvel_obcs_north = reshape(vvel_obcs_north,X_num*1,42*216);
% Vvel_obcs_east = reshape(vvel_obcs_east,Y_num*1,42*216);
% Vvel_obcs_west = reshape(vvel_obcs_west,Y_num*1,42*216);
% 
% PHIHYD_obcs_north = reshape(PHIHYD_obcs_north,X_num*1,42*216);
% PHIHYD_obcs_east = reshape(PHIHYD_obcs_east,Y_num*1,42*216);
% PHIHYD_obcs_west = reshape(PHIHYD_obcs_west,Y_num*1,42*216);
% 
% 
% 
% SIThick_obcs = reshape(SIThick_obcs,X_num*Y_num,216);
% SIArea_obcs = reshape(SIArea_obcs,738*84,1080);
SIUvel_obcs = reshape(SIUvel_obcs,X_num*Y_num,360);
SIVvel_obcs = reshape(SIVvel_obcs,X_num*Y_num,360);


% SIThick_obcs_north = reshape(SIThick_obcs_north,X_num*1,216);
% SIThick_obcs_east = reshape(SIThick_obcs_east,Y_num*1,216);
% SIThick_obcs_west = reshape(SIThick_obcs_west,Y_num*1,216);

% SIArea_obcs_north = reshape(SIArea_obcs_north,738*1,1080);
% SIArea_obcs_east = reshape(SIArea_obcs_east,84*1,1080);
% SIArea_obcs_west = reshape(SIArea_obcs_west,84*1,1080);

% SIUvel_obcs_north = reshape(SIUvel_obcs_north,X_num*1,360);
% SIUvel_obcs_east = reshape(SIUvel_obcs_east,Y_num*1,360);
% SIUvel_obcs_west = reshape(SIUvel_obcs_west,Y_num*1,360);
% 
% SIVvel_obcs_north = reshape(SIVvel_obcs_north,X_num*1,360);
% SIVvel_obcs_east = reshape(SIVvel_obcs_east,Y_num*1,360);
% SIVvel_obcs_west = reshape(SIVvel_obcs_west,Y_num*1,360);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WRITING FILES TO INPUT PATH                     %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% data = theta_obcs_north;
% writeDataset(data,fullfile(input,'Theta_northforty.bin'),ieee,prec);
% clear data
% 
% data = theta_obcs_west;
% writeDataset(data,fullfile(input,'Theta_westforty.bin'),ieee,prec);
% clear data
% 
% data = theta_obcs_east;
% writeDataset(data,fullfile(input,'Theta_eastforty.bin'),ieee,prec);
% clear data
% 
% 
% data = salt_obcs_north;
% writeDataset(data,fullfile(input,'Salt_northforty.bin'),ieee,prec);
% clear data
% 
% data = salt_obcs_east;
% writeDataset(data,fullfile(input,'Salt_eastforty.bin'),ieee,prec);
% clear data
% 
% data = salt_obcs_west;
% writeDataset(data,fullfile(input,'Salt_westforty.bin'),ieee,prec);
% clear data
% 
% 
% 
% 
% 
% data = Vvel_obcs_north;
% writeDataset(data,fullfile(input,'Vvel_northforty.bin'),ieee,prec);
% clear data
% 
% data = Vvel_obcs_east;
% writeDataset(data,fullfile(input,'Vvel_eastforty.bin'),ieee,prec);
% clear data
% 
% data = Vvel_obcs_west;
% writeDataset(data,fullfile(input,'Vvel_westforty.bin'),ieee,prec);
% clear data
% 
% 
% 
% 
% 
% data = Uvel_obcs_east;
% writeDataset(data,fullfile(input,'Uvel_eastforty.bin'),ieee,prec);
% clear data
% 
% data = Uvel_obcs_west;
% writeDataset(data,fullfile(input,'Uvel_westforty.bin'),ieee,prec);
% clear data
% 
% data = Uvel_obcs_north;
% writeDataset(data,fullfile(input,'Uvel_northforty.bin'),ieee,prec);
% clear data
% 
% 
% 
% 
% 
% data = PHIHYD_obcs_north;
% writeDataset(data,fullfile(input,'PHIHYD_northforty.bin'),ieee,prec);
% clear data
% 
% data = PHIHYD_obcs_east;
% writeDataset(data,fullfile(input,'PHIHYD_eastforty.bin'),ieee,prec);
% clear data
% 
% data = PHIHYD_obcs_west;
% writeDataset(data,fullfile(input,'PHIHYD_westforty.bin'),ieee,prec);
% clear data



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%SEA ICE%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


% data = SIThick_obcs;
% writeDataset(data,fullfile(input,'SIThickforty.bin'),ieee,prec);
% clear data
% % 
% 

% data = SIArea_obcs;
% writeDataset(data,fullfile(input,'Areaforty.bin'),ieee,prec);
% clear data




data = SIUvel_obcs;
writeDataset(data,fullfile(input,'uIceforty.bin'),ieee,prec);
clear data



data = SIVvel_obcs;
writeDataset(data,fullfile(input,'vIceforty.bin'),ieee,prec);
clear data


% data = SIThick_obcs_north;
% writeDataset(data,fullfile(input,'SIThick_northforty.bin'),ieee,prec);
% clear data
% 
% 
% data = SIThick_obcs_east;
% writeDataset(data,fullfile(input,'SIThick_eastforty.bin'),ieee,prec);
% clear data
% 
% data = SIThick_obcs_west;
% writeDataset(data,fullfile(input,'SIThick_westforty.bin'),ieee,prec);
% clear data



% data = SIArea_obcs_north;
% writeDataset(data,fullfile(input,'Area_northforty.bin'),ieee,prec);
% clear data
% 
% 
% data = SIArea_obcs_west;
% writeDataset(data,fullfile(input,'Area_westforty.bin'),ieee,prec);
% clear data
% 
% 
% data = SIArea_obcs_east;
% writeDataset(data,fullfile(input,'Area_eastforty.bin'),ieee,prec);
% clear data




% data = SIUvel_obcs_north;
% writeDataset(data,fullfile(input,'uIce_northforty.bin'),ieee,prec);
% clear data
% 
% data = SIUvel_obcs_east;
% writeDataset(data,fullfile(input,'uIce_eastforty.bin'),ieee,prec);
% clear data
% 
% data = SIUvel_obcs_west;
% writeDataset(data,fullfile(input,'uIce_westforty.bin'),ieee,prec);
% clear data
% 
% 
% 
% data = SIVvel_obcs_north;
% writeDataset(data,fullfile(input,'vIce_northforty.bin'),ieee,prec);
% clear data
% 
% data = SIVvel_obcs_west;
% writeDataset(data,fullfile(input,'vIce_westforty.bin'),ieee,prec);
% clear data
% 
% data = SIVvel_obcs_east;
% writeDataset(data,fullfile(input,'vIce_eastforty.bin'),ieee,prec);
% clear data
% 
% 
% 
% 
