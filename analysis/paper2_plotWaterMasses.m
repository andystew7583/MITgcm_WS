% %%%
% %%% paper2_plotWaterMasses.m
% %%%
% %%% Plots simulated water mass distribution and compares with observations
% %%%
% 
% %%%%%%%%%%%%%%%%
% %%%%% DATA %%%%%
% %%%%%%%%%%%%%%%%
% 
% %%% Load experiment data
% expdir = '../experiments';
% % expname = 'hires_seq_onethird_RTOPO2';
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% % expname = 'hires_seq_onetwelfth_notides_RTOPO2';
% loadexp;
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%
% %%%%% OPTIONS %%%%%
% %%%%%%%%%%%%%%%%%%%
% 
% %%% T/S grids
% dT = 0.025;
% dS = 0.005;
% Tmin = -3;
% Tmax = 2;
% Smin = 33.7;
% Smax = 35;
% SS = Smin:dS:Smax;
% TT = Tmin:dT:Tmax;
% NS = length(SS);
% NT = length(TT);
% 
% %%% Quasi-latitude horizontal grid
% ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,false);
% ETA = repmat(ETA,[1 1 Nr]);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% LOAD REFERENCE STATE %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Load data file
% outfname = [expname,'_TSfluxes'];
% outfname = [outfname,'.mat'];
% pt_ref = load(fullfile('products',outfname),'theta_tavg');
% pt_ref = pt_ref.theta_tavg;
% ss_ref = load(fullfile('products',outfname),'salt_tavg');
% ss_ref = ss_ref.salt_tavg;
% 
% %%% Pressure is just Boussinesq hydrostatic reference pressure
% pp_ref = -gravity*rhoConst*repmat(RC,[Nx Ny 1])/1e4;
% 
% %%% Remove topography
% pt_ref(hFacC==0) = NaN;
% ss_ref(hFacC==0) = NaN;
% 
% %%% Potential density
% pd_ref = densjmd95(ss_ref,pt_ref,pp_ref);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% COMPUTATIONS %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Volumes in T/S bins
% TSvol_full = binByVolume(ss_ref,pt_ref,[], ...
%                   Smin,Smax,dS,Tmin,Tmax,dT, ...
%                   RAC,DRF,hFacC);               
%                                              
% %%% Calculate density
% sigma = densmdjwf(ss_ref,pt_ref,500*ones(size(pt_ref)));
% 
% %%% Load CTD data
% load_NARE85_IC92;
                
                


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE THE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
fontsize = 14;
framepos = [417    34   681   926];
              

figure(203);                
clf;
set(gcf,'Position',framepos);            
                
                



%%%
%%% T/S DIAGRAM
%%%

[TTT,SSS] = meshgrid(TT,SS);
DDD = densjmd95(SSS,TTT,-RC(1)*gravity*rhoConst/1e4*ones(NS,NT)) - 1000;

subplot(3,2,1);
pcolor(SSS,TTT,log10(TSvol_full));
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
% pcolor(SSS,TTT,(TSvol));
shading flat;
colorbar;
% colormap(gca,flip(haxby));
colormap(gca,cmocean('amp'));
caxis([10 14]);





%%%
%%% MODEL, FILCHNER OUTFLOW
%%%

%%% Meshgrid for plots
[ZZ,XX]=meshgrid(zz,XC(:,1));
jidx = find(YC(1,:)>-74.61,1,'first');
dmax = 650;

%%% Make plots
subplot(3,2,3);
colormap jet;
pcolor(XX,-ZZ,squeeze(pt_ref(:,jidx,:)))
shading interp
caxis([-2.2 -1.6])
colormap redblue
axis([-37.16 -30.99 0 dmax])
colorbar
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);







%%% 
%%% OBS, FILCHNER OUTTTFLOW
%%%

subplot(3,2,4);
pcolor(NARE85_HH,NARE85_ZZ,NARE85_pt)
shading interp;
colorbar;
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);
axis([-37.16 -30.99 0 dmax]);





%%% 
%%% MODEL, WESTERN WEDDELL
%%%

%%% Meshgrid for plots
[ZZ,XX]=meshgrid(zz,XC(:,1));
jidx = find(YC(1,:)>-67.62,1,'first');
dmax = 3000;

%%% Make plots
subplot(3,2,5);
colormap jet;
pcolor(XX,-ZZ,squeeze(pt_ref(:,jidx,:)))
shading interp
% caxis([-2.2 -1.6])
colormap redblue
axis([-57.52 -52.07 0 dmax])
colorbar
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);






subplot(3,2,6);
pcolor(IC92_HH,IC92_ZZ,IC92_pt)
shading interp;
colorbar;
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);
axis([-57.52 -52.07 0 dmax]);

