%%%
%%% plotFilchnerOutflow.m
%%%
%%% Plots T and S sections along the Filchner outflow
%%%


%%% Read experiment data
loadexp;

%%% Set time range
% tmin = (18-0.05)*86400*365;
% tmax = (18+0.05)*86400*365;
tmin = (5.5-0.05)*86400*365;
tmax = (5.5+0.05)*86400*365;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(4));
nDumps = round(nTimeSteps*(deltaT/dumpFreq));
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Load mean T and S
exppath = fullfile(expdir,expname);
thetawc = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
saltwc = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
vvelwc = readIters(exppath,'VVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
saltwc(hFacC==0) = NaN;
thetawc(hFacC==0) = NaN;
vvelwc(hFacS==0) = NaN;

%%% Calculate density
denswc = densmdjwf(saltwc,thetawc,500*ones(size(thetawc)));

%%% Meshgrid for plots
[ZZ,XX]=meshgrid(zz,XC(:,1));
jidx = 558;
YC(1,jidx)
dmax = 1500;

%%% Plot bathy for reference
figure(9);
pcolor(XC,YC,bathy);
shading interp;
colorbar;
colormap haxby(20);
caxis([-2000 0]);
hold on;
plot(XC(:,1),YC(1,jidx)*ones(Nx,1),'k--');
hold off;

%%% Make plots
figure(10);
colormap jet;
pcolor(XX,-ZZ,squeeze(thetawc(:,jidx,:)))
shading interp
caxis([-2.2 -1.6])
colormap redblue
axis([-64 -15 0 dmax])
colorbar
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);

figure(11)
pcolor(XX,-ZZ,squeeze(saltwc(:,jidx,:)))
shading interp
axis([-64 -15 0 dmax])
colormap jet
colorbar
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');

figure(12)
pcolor(XX,-ZZ,squeeze(denswc(:,jidx,:)))
shading interp
axis([-64 -15 0 dmax])
colormap jet
colorbar
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');

figure(13)
pcolor(XX,-ZZ,squeeze(vvelwc(:,jidx,:)))
shading interp
axis([-64 -15 0 dmax])
colormap redblue
colorbar
caxis([-.25 .25]);
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');

imin = find(XC(:,1)>-45,1);
imax = find(XC(:,1)>-27,1);
figure(14);
theta_f = squeeze(thetawc(imin:imax,jidx,:));
salt_f = squeeze(saltwc(imin:imax,jidx,:));
zz_f = squeeze(ZZ(imin:imax,:));
scatter(salt_f(:),theta_f(:),10,zz_f(:));
caxis([-1000 0]);
