%%%
%%% plotTSdiagram.m
%%%
%%% Plots a volumetric T/S diagram of our model domain.
%%%

%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
% expname = 'hires_seq_onethird_RTOPO2';
expname = 'hires_seq_onesixth_notides_RTOPO2';
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
loadexp;

% %%% Time frame over which to average thermodynamic variables to create
% %%% climatology
% % tmin = 18.05*86400*365;
% % tmax = 27.05*86400*365;
% tmin = 9.05*86400*365;
% tmax = 18.05*86400*365;




%%%%%%%%%%%%%%%%%%%
%%%%% OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%

%%% T/S grids
dT = 0.025;
dS = 0.005;
Tmin = -3;
Tmax = 2;
Smin = 33.7;
Smax = 35;
SS = Smin:dS:Smax;
TT = Tmin:dT:Tmax;
NS = length(SS);
NT = length(TT);

%%% Quasi-latitude horizontal grid
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,false);
ETA = repmat(ETA,[1 1 Nr]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD REFERENCE STATE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% Frequency of diagnostic output - should match that specified in
% %%% data.diagnostics.
% dumpFreq = abs(diag_frequency(1));
% nDumps = round(endTime/dumpFreq);
% dumpIters = round((1:nDumps)*dumpFreq/deltaT);
% dumpIters = dumpIters(dumpIters > nIter0);
% nDumps = length(dumpIters);
% 
% %%% Time-average temperature and salinity
% pt_ref = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
% ss_ref = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);

%%% Load data file
outfname = [expname,'_TSfluxes'];
outfname = [outfname,'.mat'];
pt_ref = load(fullfile('products',outfname),'theta_tavg');
pt_ref = pt_ref.theta_tavg;
ss_ref = load(fullfile('products',outfname),'salt_tavg');
ss_ref = ss_ref.salt_tavg;

%%% Pressure is just Boussinesq hydrostatic reference pressure
pp_ref = -gravity*rhoConst*repmat(RC,[Nx Ny 1])/1e4;

%%% Remove topography
pt_ref(hFacC==0) = NaN;
ss_ref(hFacC==0) = NaN;

%%% Potential density
pd_ref = densjmd95(ss_ref,pt_ref,pp_ref);















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% COMPUTE VOLUMES IN T/S BINS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Do the binning
TSvol_full = binByVolume(ss_ref,pt_ref,[], ...
                  Smin,Smax,dS,Tmin,Tmax,dT, ...
                  RAC,DRF,hFacC);               
                
%%% Version for southern WS continental shelf and cavity only     
ss_ref_shelf = ss_ref;
pt_ref_shelf = pt_ref;
icedraft3D = repmat(SHELFICEtopo,[1 1 Nr]);
bathy3D = repmat(bathy,[1 1 Nr]);
ss_ref_shelf((icedraft3D<0) | (ETA>5) | ((ETA>3) & (bathy3D<-1500))) = NaN;
pt_ref_shelf((icedraft3D<0) | (ETA>5) | ((ETA>3) & (bathy3D<-1500))) = NaN;
TSvol_shelf = binByVolume(ss_ref_shelf,pt_ref_shelf,[], ...
                  Smin,Smax,dS,Tmin,Tmax,dT, ...
                  RAC,DRF,hFacC);
                
%%% Version for northern restoring regions only
ss_ref_sponge = ss_ref;
pt_ref_sponge = pt_ref;
XC3D = repmat(XC,[1 1 Nr]);
YC3D = repmat(YC,[1 1 Nr]);
ss_ref_sponge((YC3D<YC(1,end-spongethickness+1) & (XC3D<XC(end-spongethickness+1,1)))) = NaN;
pt_ref_sponge((YC3D<YC(1,end-spongethickness+1) & (XC3D<XC(end-spongethickness+1,1)))) = NaN;
TSvol_sponge = binByVolume(ss_ref_sponge,pt_ref_sponge,[], ...
                  Smin,Smax,dS,Tmin,Tmax,dT, ...
                  RAC,DRF,hFacC);

%%% NaN out any areas with zero volume
TSvol_full(TSvol_full==0) = NaN;
TSvol_shelf(TSvol_shelf==0) = NaN;
TSvol_sponge(TSvol_sponge==0) = NaN;

%%% Bin with weighting by eta
TSvol_ETA_full = binByVolume(ss_ref,pt_ref,ETA, ...
                  Smin,Smax,dS,Tmin,Tmax,dT, ...
                  RAC,DRF,hFacC);
TSvol_ETA_full = TSvol_ETA_full ./ TSvol_full;








%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%

[TTT,SSS] = meshgrid(TT,SS);
DDD = densjmd95(SSS,TTT,-RC(1)*gravity*rhoConst/1e4*ones(NS,NT)) - 1000;

figure(61);
pcolor(SSS,TTT,log10(TSvol_full));
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
% pcolor(SSS,TTT,(TSvol));
shading flat;
colorbar;
colormap(flip(haxby));
caxis([10 14]);


figure(62);
pcolor(SSS,TTT,log10(TSvol_shelf));
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
% pcolor(SSS,TTT,(TSvol));
shading flat;
colorbar;
colormap(flip(haxby));
caxis([10 14]);

figure(63);
pcolor(SSS,TTT,log10(TSvol_sponge));
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
% pcolor(SSS,TTT,(TSvol));
shading flat;
colorbar;
colormap(flip(haxby));
caxis([10 14]);

figure(64);
pcolor(SSS,TTT,TSvol_ETA_full);
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
% pcolor(SSS,TTT,(TSvol));
shading flat;
colorbar;
colormap(flip(haxby));
% caxis([10 14]);

figure(65)
ss_bdry = [squeeze(ss_ref(:,Ny,:)) ; squeeze(ss_ref(Nx,:,:))];
pt_bdry = [squeeze(pt_ref(:,Ny,:)) ; squeeze(pt_ref(Nx,:,:))];
scatter(ss_bdry(:),pt_bdry(:));
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);





