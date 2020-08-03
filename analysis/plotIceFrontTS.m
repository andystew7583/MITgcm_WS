%%%
%%% plotIceFrontTS.m
%%%
%%% Plots T and S sections along the Filchner-Ronne ice front
%%%


%%% Read experiment data
setExpname;
expname = 'n_34452';
loadexp;

%%% Set time range
tmin = 9*86400*360;
tmax = 18*86400*360;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
% nIter0 = 0;
nTimeSteps = 1843855;
nDumps = round(nTimeSteps*(deltaT/dumpFreq));
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);
 

%%% Load mean T and S
exppath5 = '/data3/MITgcm_WS/experiments/n_34452';
thetawc = readIters(exppath5,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
saltwc = readIters(exppath5,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);

%%% Define section
saltwc(saltwc==0) = NaN;
thetawc(thetawc==0) = NaN;

figure(3);
[YY,XX] = meshgrid(yy,xx);
pcolor(XX,YY,saltwc(:,:,1));
colormap jet;


x1 = -60.59;
y1 = -74.37;
x2 = -39.38;
y2 = -78.34;
Npts = 201;
dx = (x2-x1)/(Npts-1);
dy = (y2-y1)/(Npts-1);
xxl = x1:dx:x2;
yyl = y1:dy:y2;

%%% Interpolate onto section
[XXi,YYi] = meshgrid(xx,yy);
TTl = NaN*ones(Npts,Nr);
SSl = NaN*ones(Npts,Nr);
for k=1:Nr
  TTl(:,k) = interp2(XXi,YYi,thetawc(:,:,k)',xxl,yyl,'linear');
  SSl(:,k) = interp2(XXi,YYi,saltwc(:,:,k)',xxl,yyl,'linear');
end

[ZZl,XXl] = meshgrid(zz,xxl);

figure(1);
pcolor(XXl,-ZZl,TTl);
shading interp;
axis([x1 x2 0 1000]);
colorbar;
colormap jet;
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');

figure(2);
pcolor(XXl,-ZZl,SSl);
shading interp;
axis([x1 x2 0 1000]);
colorbar;
colormap jet;
set(gca,'YDir','reverse');
xlabel('Longitude');
ylabel('Depth (m)');