%%%
%%% plotWScrosssection.m
%%%
%%% Plots a cross-section of time-mean properties through the Weddell Sea.
%%%



%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
expname = 'hires_seq_onethird_RTOPO2';
loadexp;

%%% Time frame over which to average thermodynamic variables to create
%%% climatology
tmin = 18.05*86400*365;
tmax = 27.05*86400*365;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% COMPUTE REFERENCE STATE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Time-average temperature and salinity
pt_tavg = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
ss_tavg = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);

%%% Pressure is just Boussinesq hydrostatic reference pressure
pp_ref = -gravity*rhoConst*repmat(RC,[Nx Ny 1])/1e4;

%%% Remove topography
pt_tavg(hFacC==0) = NaN;
ss_tavg(hFacC==0) = NaN;

%%% Potential densities
pd0 = densjmd95(ss_tavg,pt_tavg,0*pp_ref) - 1000;
pd5 = densjmd95(ss_tavg,pt_tavg,5000*ones(size(pp_ref))) - 1000;

%%% Neutral density
%%% TODO NEED TO CHANGE TO MATCH EXPERIMENT
load('hires_seq_onethird_ND1');







%%%%%%%%%%%%%%%%%
%%%%% PLOTS %%%%%
%%%%%%%%%%%%%%%%%

framepos = [386         211        1024         652];
fontsize = 14;
plotloc = [0.1215    0.1100    0.7245    0.8150];

%%% Pot temp
idx = find(XC(:,1)<-38,1,'last');
[ZZ,YY]=meshgrid(RC,YC(1,:));
pt_min = -3;
pt_max = 1;
pt_plot = pt_tavg;
pt_plot(hFacC==0)= NaN;
pt_plot = squeeze(pt_plot(idx,:,:));
pt_plot(pt_plot<pt_min) = pt_min;
pt_plot(pt_plot>pt_max) = pt_max;
figure(6);
clf;
set(gcf,'Position',framepos);
contourf(YY,ZZ,pt_plot,[pt_min:.1:pt_max],'EdgeColor','None');
shading interp;
colorbar;
colormap(cmocean('thermal'));
set(gca,'FontSize',fontsize);
set(gca,'Position',plotloc);
xlabel('Latitude');
ylabel('Depth (m)');
title('Potential temperatue (^oC) at 37W');
print('-dpng',['Figures/',expname,'_theta_section.png']);


%%% Salinity
idx = find(XC(:,1)<-38,1,'last');
[ZZ,YY]=meshgrid(RC,YC(1,:));
ss_min = 34.2;
ss_max = 34.9;
ss_plot = ss_tavg;
ss_plot(hFacC==0)= NaN;
ss_plot = squeeze(ss_plot(idx,:,:));
ss_plot(ss_plot<ss_min) = ss_min;
ss_plot(ss_plot>ss_max) = ss_max;
figure(7);
clf;
set(gcf,'Position',framepos);
contourf(YY,ZZ,ss_plot,[ss_min:.05:ss_max],'EdgeColor','None');
shading interp;
colorbar;
colormap(cmocean('haline'));
set(gca,'FontSize',fontsize);
set(gca,'Position',plotloc);
xlabel('Latitude');
ylabel('Depth (m)');
title('Salinity (g/kg) at 37W');
print('-dpng',['Figures/',expname,'_salt_section.png']);

%%% Sigma0
idx = find(XC(:,1)<-38,1,'last');
[ZZ,YY]=meshgrid(RC,YC(1,:));
dens_min = 27.5;
dens_max = 28.1;
dens_plot = pd0;
dens_plot(hFacC==0)= NaN;
dens_plot = squeeze(dens_plot(idx,:,:));
dens_plot(dens_plot<dens_min) = dens_min;
dens_plot(dens_plot>dens_max) = dens_max;
figure(8);
clf;
set(gcf,'Position',framepos);
contourf(YY,ZZ,dens_plot,[dens_min:.05:dens_max],'EdgeColor','None');
shading interp;
colorbar;
colormap(cmocean('dense'));
set(gca,'FontSize',fontsize);
set(gca,'Position',plotloc);
xlabel('Latitude');
ylabel('Depth (m)');
title('\sigma_0 (kg/m^3) at 37W');
print('-dpng',['Figures/',expname,'_sigma0_section.png']);

%%% Sigma5
idx = find(XC(:,1)<-38,1,'last');
[ZZ,YY]=meshgrid(RC,YC(1,:));
dens_min = 50.2;
dens_max = 50.8;
dens_plot = pd5;
dens_plot(hFacC==0)= NaN;
dens_plot = squeeze(dens_plot(idx,:,:));
dens_plot(dens_plot<dens_min) = dens_min;
dens_plot(dens_plot>dens_max) = dens_max;
figure(9);
clf;
set(gcf,'Position',framepos);
contourf(YY,ZZ,dens_plot,[dens_min:.05:dens_max],'EdgeColor','None');
shading interp;
colorbar;
colormap(cmocean('dense'));
set(gca,'FontSize',fontsize);
set(gca,'Position',plotloc);
xlabel('Latitude');
ylabel('Depth (m)');
title('\sigma_5 (kg/m^3) at 37W');
print('-dpng',['Figures/',expname,'_sigma5_section.png']);


%%% gamma
idx = find(XC(:,1)<-38,1,'last');
[ZZ,YY]=meshgrid(RC,YC(1,:));
dens_min = 27.6;
dens_max = 28.7;
dens_plot = gg_ref;
dens_plot(hFacC==0)= NaN;
dens_plot = squeeze(dens_plot(idx,:,:));
dens_plot(dens_plot<dens_min) = dens_min;
dens_plot(dens_plot>dens_max) = dens_max;
figure(10);
clf;
set(gcf,'Position',framepos);
contourf(YY,ZZ,dens_plot,[dens_min:.05:dens_max],'EdgeColor','None');
shading interp;
colorbar;
colormap(cmocean('dense'));
set(gca,'FontSize',fontsize);
set(gca,'Position',plotloc);
xlabel('Latitude');
ylabel('Depth (m)');
title('\gamma (kg/m^3) at 37W');
print('-dpng',['Figures/',expname,'_gamma_section.png']);