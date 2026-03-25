%%%%%%%
%%% Plot instananeous Temp


run  '../newexp/defineGrid.m';
%%% Read experiment data
loadexp;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(7));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

tt = zeros(1,84);


%%%%Chose an arbitrary time to output 
Theta_tot = NaN(Nx,Ny,Nr,nDumps);
SIarea_tot = NaN(Nx,Ny,nDumps);
UVEL_tot = NaN(Nx,Ny,Nr,nDumps);
VVEL_tot = NaN(Nx,Ny,Nr,nDumps);
Salt_tot = NaN(Nx,Ny,Nr,nDumps);

%%% n = 1:nDumps (but experiment hasn't finished yet %%%
for n=1:84  
  
  Salt = rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n));
  Salt_tot(:,:,:,n) = Salt; 
  
  UVEL = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),dumpIters(n));
  UVEL_tot(:,:,:,n) = UVEL;
  
  VVEL = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));
  VVEL_tot(:,:,:,n) = VVEL;
  
  Theta = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));
  Theta_tot(:,:,:,n) = Theta;
    
  SIarea = rdmdsWrapper(fullfile(exppath,'/results/SIarea_inst'),dumpIters(n));
  SIarea_tot(:,:,n) = SIarea;
  
  tt(n) =  (dumpIters(n)*deltaT)/86400;
  %tt(n) is the day of the run
  
end

Nyears = 84/12;

%%%%%%%%%%%% Now we want to get seasons
SeasTemp = zeros(size(Theta_tot,1),size(Theta_tot,2),size(Theta_tot,3),Nyears,4);
SeasSalt = zeros(size(Salt_tot,1),size(Salt_tot,2),size(Salt_tot,3),Nyears,4);
SeasUVEL = zeros(size(UVEL_tot,1),size(UVEL_tot,2),size(UVEL_tot,3),Nyears,4);
SeasVVEL = zeros(size(VVEL_tot,1),size(VVEL_tot,2),size(VVEL_tot,3),Nyears,4);
SeasSIarea = zeros(size(SIarea_tot,1),size(SIarea_tot,2),Nyears,4);


for n=1:Nyears
  SeasTemp(:,:,:,n,1) = mean(Theta_tot(:,:,:,(n-1)*12+[12 1 2]),4);
  SeasTemp(:,:,:,n,2) = mean(Theta_tot(:,:,:,(n-1)*12+(3:5)),4);
  SeasTemp(:,:,:,n,3) = mean(Theta_tot(:,:,:,(n-1)*12+(6:8)),4);
  SeasTemp(:,:,:,n,4) = mean(Theta_tot(:,:,:,(n-1)*12+(9:11)),4);
 
  SeasSalt(:,:,:,n,1) = mean(Salt_tot(:,:,:,(n-1)*12+[12 1 2]),4);
  SeasSalt(:,:,:,n,2) = mean(Salt_tot(:,:,:,(n-1)*12+(3:5)),4);
  SeasSalt(:,:,:,n,3) = mean(Salt_tot(:,:,:,(n-1)*12+(6:8)),4);
  SeasSalt(:,:,:,n,4) = mean(Salt_tot(:,:,:,(n-1)*12+(9:11)),4);
  
  SeasUVEL(:,:,:,n,1) = mean(UVEL_tot(:,:,:,(n-1)*12+[12 1 2]),4);
  SeasUVEL(:,:,:,n,2) = mean(UVEL_tot(:,:,:,(n-1)*12+(3:5)),4);
  SeasUVEL(:,:,:,n,3) = mean(UVEL_tot(:,:,:,(n-1)*12+(6:8)),4);
  SeasUVEL(:,:,:,n,4) = mean(UVEL_tot(:,:,:,(n-1)*12+(9:11)),4);
  
  SeasVVEL(:,:,:,n,1) = mean(VVEL_tot(:,:,:,(n-1)*12+[12 1 2]),4);
  SeasVVEL(:,:,:,n,2) = mean(VVEL_tot(:,:,:,(n-1)*12+(3:5)),4);
  SeasVVEL(:,:,:,n,3) = mean(VVEL_tot(:,:,:,(n-1)*12+(6:8)),4);
  SeasVVEL(:,:,:,n,4) = mean(VVEL_tot(:,:,:,(n-1)*12+(9:11)),4);
    
  SeasSIarea(:,:,n,1) = mean(SIarea_tot(:,:,(n-1)*12+[12 1 2]),3);
  SeasSIarea(:,:,n,2) = mean(SIarea_tot(:,:,(n-1)*12+(3:5)),3);
  SeasSIarea(:,:,n,3) = mean(SIarea_tot(:,:,(n-1)*12+(6:8)),3);
  SeasSIarea(:,:,n,4) = mean(SIarea_tot(:,:,(n-1)*12+(9:11)),3);
  
end

SeasTemp = squeeze(nanmean(SeasTemp,4));
SeasSalt = squeeze(nanmean(SeasSalt,4));
SeasSIarea = squeeze(nanmean(SeasSIarea,3));
SeasUVEL = squeeze(nanmean(SeasUVEL,4));
SeasVVEL = squeeze(nanmean(SeasVVEL,4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make a plot of temperature


%%%%%% Choose a depth ~ Surface and 500m



ThetaJJA = squeeze(SeasTemp(244,:,:,3));
% hFacCsfc=hFacC(:,:,1);
% ThetaJJA(hFacCsfc==0)=NaN;
% ThetaJJA_EB = squeeze(SeasTemp(end,:,:,1));
ThetaJJA_NB = squeeze(SeasTemp(:,end,:,3));
ThetaDJF = squeeze(SeasTemp(244,:,:,1));
ThetaDJF_NB = squeeze(SeasTemp(:,end,:,1));

SaltDJF = squeeze(SeasSalt(244,:,:,1));
SaltDJF_NB = squeeze(SeasSalt(:,end,:,1));


SaltJJA = squeeze(SeasSalt(244,:,:,3));
SaltJJA_NB = squeeze(SeasSalt(:,end,:,3));

% hFacCsfc=hFacC(:,:,1);
% SaltDJF(hFacCsfc==0)=NaN;
% SaltDJF_EB = squeeze(SeasSalt(end,:,:,1));



% UVELDJF = SeasUVEL(1:244,:,1,1);
% UVELDJF(hFacCsfc==0)=NaN;
% 
% UVELJJA = SeasUVEL(244,:,1,3);
% UVELJJA(hFacCsfc==0)=NaN;
% 
% VVELDJF = SeasVVEL(244,:,1,1);
% VVELDJF(hFacCsfc==0)=NaN;
% 
% VVELJJA = SeasVVEL(244,:,1,3);
% VVELJJA(hFacCsfc==0)=NaN;
% 
% SIareaDJF = SeasSIarea(244,:,1);
% SIareaDJF(hFacCsfc==0)=NaN;
% 
% SIareaJJA = SeasSIarea(:,:,3);
% SIareaJJA(hFacCsfc==0)=NaN;
% 
% SIareaSON = SeasSIarea(244,:,4);
% SIareaSON(hFacCsfc==0)=NaN;
% 
% SIareaMAM = SeasSIarea(244,:,2);
% SIareaMAM(hFacCsfc==0)=NaN;
% 
% 
% 
% %%%%%%%%%% 500m
% Theta500DJF = SeasTemp(244,:,45,1);
% hFacC500=hFacC(244,:,45);
% Theta500DJF(hFacC500==0)=NaN;
% 
% Theta500JJA = SeasTemp(244,:,45,3);
% Theta500JJA(hFacC500==0)=NaN;
% 
% UVEL500DJF = SeasUVEL(244,:,45,1);
% UVEL500DJF(hFacC500==0)=NaN;
% 
% UVEL500JJA = SeasUVEL(244,:,45,3);
% UVEL500JJA(hFacC500==0)=NaN;
% 
% VVEL500DJF = SeasVVEL(244,:,45,1);
% VVEL500DJF(hFacC500==0)=NaN;
% 
% VVEL500JJA = SeasVVEL(244,:,45,3);
% VVEL500JJA(hFacC500==0)=NaN;


%%%%%%% PloTs

% lat vs depth -> EB
[YR,ZR]=meshgrid(ymc,zz);

% lon vs depth -> NB 
[XX,ZX]=meshgrid(xmc,zz);
%%%%%%%%%% Theta

figure(1)
pcolor(YR,ZR,ThetaJJA'),shading interp
colorbar
colormap jet(100)
title('Theta JJA A12 Transect Expanded');
print(figure(1),'-dpng',['/data1/MITgcm_WS/newexp/Expanded/ThetaJJA_A12.png']);

figure(41)
pcolor(XX,ZX,ThetaJJA_NB'),shading interp
colorbar
colormap jet(100)
title('Theta JJA SR4  Expanded');
print(figure(41),'-dpng',['/data1/MITgcm_WS/newexp/Expanded/ThetaJJA_SR4.png']);

figure(40)
pcolor(YR,ZR,ThetaDJF'),shading interp
colorbar
colormap jet(100)
title('Theta DJF A12  Expanded');
print(figure(1),'-dpng',['/data1/MITgcm_WS/newexp/Expanded/ThetaDJF_A12.png']);

figure(20)
pcolor(XX,ZX,ThetaDJF_NB'),shading interp
colorbar
colormap jet(100)
title('Theta DJF SR4');
xlabel('Latitude')
print(figure(20),'-dpng',['/data1/MITgcm_WS/newexp/Expanded/ThetaDJF_SR4.png']);

figure(21)
pcolor(YR,ZR,SaltDJF'),shading interp
colorbar
colormap jet(100)
caxis([33 35])
title('Salt DJF A12 Expanded');
xlabel('Longitude')
print(figure(21),'-dpng',['/data1/MITgcm_WS/newexp/Expanded/SaltDJF_A12.png']);

figure(22)
pcolor(XX,ZX,SaltDJF_NB'),shading interp
colorbar
colormap jet(100)
title('Theta DJF SR4');
xlabel('Latitude')
print(figure(22),'-dpng',['/data1/MITgcm_WS/newexp/Expanded/SaltDJF_SR4.png']);

figure(23)
pcolor(XX,ZX,SaltJJA_NB'),shading interp
colorbar
colormap jet(100)
title('Salt JJA SR4');
xlabel('Latitude')
print(figure(23),'-dpng',['/data1/MITgcm_WS/newexp/Expanded/SaltJJA_SR4.png']);


figure(2)
pcolor(YR,ZR,SaltJJA'),shading interp
colorbar
colormap jet(100)
caxis([33 35])
title('Salt JJA A12 Expanded');
print(figure(2),'-dpng',['/data1/MITgcm_WS/newexp/Expanded/SaltJJA_A12.png']);

% figure(3)
% pcolor(XC,YC,ThetaJJA),shading interp
% colorbar
% colormap jet(100)
% title('Theta JJA Surface');
% print(figure(3),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/ThetaJJA.png']);
% 
% figure(4)
% pcolor(XC,YC,Theta500JJA),shading interp
% colorbar
% colormap jet(100)
% title('Theta JJA 500m');
% print(figure(4),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/ThetaJJA500m.png']);
% 
% %%%%%%% UVEL
% 
% figure(5)
% pcolor(XC,YC,UVELDJF),shading interp
% colorbar
% colormap jet(100)
% title('UVEL DJF Surface');
% print(figure(5),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/UVELDJF.png']);
% 
% figure(6)
% pcolor(XC,YC,UVEL500DJF),shading interp
% colorbar
% colormap jet(100)
% title('UVEL DJF 500m');
% print(figure(6),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/UVELDJF500m.png']);
% 
% figure(7)
% pcolor(XC,YC,UVELJJA),shading interp
% colorbar
% colormap jet(100)
% title('UVEL JJA Surface');
% print(figure(7),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/UVELJJA.png']);
% 
% figure(8)
% pcolor(XC,YC,UVEL500JJA),shading interp
% colorbar
% colormap jet(100)
% title('UVEL JJA 500m');
% print(figure(8),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/UVELJJA500m.png']);
% 
% 
% %%%%%%%%% VVEL
% 
% figure(9)
% pcolor(XC,YC,VVELDJF),shading interp
% colorbar
% colormap jet(100)
% title('VVEL DJF Surface');
% print(figure(9),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/VVELDJF.png']);
% 
% figure(10)
% pcolor(XC,YC,VVEL500DJF),shading interp
% colorbar
% colormap jet(100)
% title('VVEL DJF 500m');
% print(figure(10),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/VVELDJF500m.png']);
% 
% figure(11)
% pcolor(XC,YC,VVELJJA),shading interp
% colorbar
% colormap jet(100)
% title('VVEL JJA Surface');
% print(figure(11),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/VVELJJA.png']);
% 
% figure(12)
% pcolor(XC,YC,VVEL500JJA),shading interp
% colorbar
% colormap jet(100)
% title('VVEL JJA 500m');
% print(figure(12),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/VVELJJA500m.png']);
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% SIArea
% 
% figure(13)
% pcolor(XC,YC,SIareaDJF),shading interp
% colorbar
% colormap jet
% title('SIarea DJF');
% print(figure(13),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/SIareaDJF.png']);
% 
% figure(14)
% pcolor(XC,YC,SIareaMAM),shading interp
% colorbar
% colormap jet
% title('SIarea MAM');
% print(figure(14),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/SIareaMAM.png']);
% 
% 
% figure(15)
% pcolor(XC,YC,SIareaJJA),shading interp
% colorbar
% colormap jet
% title('SIarea JJA');
% print(figure(15),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/SIareaJJA.png']);
% 
% figure(16)
% pcolor(XC,YC,SIareaSON),shading interp
% colorbar
% colormap jet
% title('SIarea SON');
% print(figure(16),'-dpng',['/data1/MITgcm_WS/newexp/30yearplots/SIareaSON.png']);
% 
% 



