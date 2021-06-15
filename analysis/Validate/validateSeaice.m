%%%
%%% validateSeaice.m
%%%
%%% Validates model output against NODC sea ice concentration data.
%%%


%%% Load experiment
expdir = '/data3/MITgcm_WS/experiments';
% expname = 'hires_seq_onethird_RTOPO2';
% tmin = 19.05;
% tmax = 27.05;
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% tmin = 10.05;
% tmax = 18.05;
expname = 'hires_seq_onetwelfth_notides_RTOPO2';
tmin = 1.05;
tmax = 9.05;
cd ..
loadexp;
cd Validate











%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD CLIMATOLOGY %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load NODC data
load(['NODCiceconc.mat']);

%%% Create monthly climatology
NODCiceconc = nanmean(NODCiceconc,4);

%%% Interpolate to model grid
LA = reshape(lat,790*830,1);
LO = reshape(lon,790*830,1);
interpolatedIce = zeros(Nx,Ny,12);
for months = 1:12
  V = NODCiceconc(:,:,months);
  V = reshape(V,790*830,1);
  F = scatteredInterpolant(double(LO),double(LA),double(V),'natural');
  interpolatedIce(:,:,months) = F(double(XC),double(YC));
end
            
   
%%% Seasonal climatology
SeasSIconc = zeros(Nx,Ny,4);
SeasSIconc(:,:,1) = nanmean(interpolatedIce(:,:,[12 1 2]),3);
SeasSIconc(:,:,2) = nanmean(interpolatedIce(:,:,(3:5)),3);
SeasSIconc(:,:,3) = nanmean(interpolatedIce(:,:,(6:8)),3);
SeasSIconc(:,:,4) = nanmean(interpolatedIce(:,:,(9:11)),3);

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRE-DETERMINE ITERATION NUMBERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(7));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Determine iteration numbers to process
itersToRead = [];
times = [];
for n=1:length(dumpIters)
 
  tyears = dumpIters(n)*deltaT/86400/365;
 
  if ((tyears >= tmin) && (tyears <= tmax))    
    itersToRead = [itersToRead dumpIters(n)];
    times = [times dumpIters(n)*deltaT];
  end
  
end
Ntime = length(itersToRead);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD MODEL SEA ICE CLIMATOLOGY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load monthly concentrations
SIconc = NaN(Nx,Ny,Ntime);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))]      
  
  %%% Read velocity field
  SIconc(:,:,n) = rdmdsWrapper(fullfile(exppath,'/results/SIarea'),itersToRead(n));
  
end

%%% Seasonal climatology
ModelSeasSIconc = zeros(Nx,Ny,Ntime/12,4);
for n=1:Ntime/12
  ModelSeasSIconc(:,:,n,1) = mean(SIconc(:,:,(n-1)*12+[12 1 2]),3);   
  ModelSeasSIconc(:,:,n,2) = mean(SIconc(:,:,(n-1)*12+(3:5)),3);
  ModelSeasSIconc(:,:,n,3) = mean(SIconc(:,:,(n-1)*12+(6:8)),3);
  ModelSeasSIconc(:,:,n,4) = mean(SIconc(:,:,(n-1)*12+(9:11)),3);
end
ModelSeasSIconc = ModelSeasSIconc*100;
Mod_Ice_seas = squeeze(nanmean(ModelSeasSIconc,3));






%%%% comparing JJA
figure(3)
Ice_seas_JJA=(SeasSIconc(:,:,3));
Mod_Ice_seas_JJA=(Mod_Ice_seas(:,:,3));

Mod_Ice_seas_JJA(hFacC(:,:,1)<1)=NaN;
Ice_seas_JJA(hFacC(:,:,1)<1)=NaN;
Ice_seas_JJA(Ice_seas_JJA==0)=NaN;

Ice_seas_DJF=(SeasSIconc(:,:,1));
Mod_Ice_seas_DJF=(Mod_Ice_seas(:,:,1));

Mod_Ice_seas_DJF(hFacC(:,:,1)<1)=NaN;
Ice_seas_DJF(hFacC(:,:,1)<1)=NaN;
Ice_seas_DJF(Ice_seas_DJF==0)=NaN;




%%% Plotting options
latMin = -78.5;
latMax = YC(1,end-spongethickness);
lonMin = min(min(XC));
lonMax = XC(end-spongethickness,1);






figure(2)
clf
set(gcf,'Position',[464    68   993   914]);
ww1=subplot(2,1,1)


axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)
pcolorm(YC,XC,(Mod_Ice_seas_DJF)),shading interp;
set(gca,'Position',[0.2 0.5 .6 .45]);

e =colorbar;
set(e,'Position',[0.93 0.5 0.02 .45])
caxis([0 100])
h=title('Model DJF Sea Ice Concentration (\%cover) ','interpreter','latex','fontsize',18);
set(h,'Position',[0 -.93 0])
colormap(ww1,brewermap(50,'PuBu'))
xlabel('Longitude','interpreter','latex','FontSize',12);
ylabel('Latitude','interpreter','latex','FontSize',12);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom+.02 ax_width+.02 ax_height-.1];
hold on

ww=subplot(2,1,2)

axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)
pcolorm(YC,XC,(Mod_Ice_seas_DJF-Ice_seas_DJF)),shading interp;
caxis([-50 50]);
set(gca,'Position',[0.2 0.0 .6 .45]);
h=title('DJF Sea Ice Concentration Anomaly (Model-Observed) (\%cover) ','interpreter','latex','fontsize',18);
set(h,'Position',[0 -.94 0])
h = colorbar;
set(h,'Position',[0.93 0.01 0.02 .45])
colormap(ww,brewermap(50,'PuOr'))
xlabel('Longitude','interpreter','latex','FontSize',12);
ylabel('Latitude','interpreter','latex','FontSize',12);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom+.02 ax_width+.02 ax_height-.1];











figure(3)
clf
set(gcf,'Position',[464    68   993   914]);

ww1=subplot(2,1,1)
axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)
pcolorm(YC,XC,(Mod_Ice_seas_JJA)),shading interp;
set(gca,'Position',[0.2 0.5 .6 .45]);

e =colorbar;
set(e,'Position',[0.93 0.5 0.02 .45])
caxis([0 100])
h=title('Model JJA Sea Ice Concentration (\%cover) ','interpreter','latex','fontsize',18);
set(h,'Position',[0 -.93 0])
colormap(ww1,brewermap(50,'PuBu'))
xlabel('Longitude','interpreter','latex','FontSize',12);
ylabel('Latitude','interpreter','latex','FontSize',12);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom+.02 ax_width+.02 ax_height-.1];
hold on

ww=subplot(2,1,2)

axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin latMax], ...
  'MapLonLimit',[lonMin lonMax], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)
pcolorm(YC,XC,(Mod_Ice_seas_JJA-Ice_seas_JJA));shading interp;
caxis([-50 50]);
set(gca,'Position',[0.2 0.0 .6 .45]);
h=title('JJA Sea Ice Concentration Anomaly (Model-Observed) (\%cover) ','interpreter','latex','fontsize',18);
set(h,'Position',[0 -.94 0])
h = colorbar;
set(h,'Position',[0.93 0.01 0.02 .45])
colormap(ww,brewermap(50,'PuOr'))
xlabel('Longitude','interpreter','latex','FontSize',12);
ylabel('Latitude','interpreter','latex','FontSize',12);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom+.02 ax_width+.02 ax_height-.1];