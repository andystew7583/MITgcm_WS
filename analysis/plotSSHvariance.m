%%%
%%% calcSSHvariance.m
%%%
%%% Calculates variance of SSH using instantaneous model output.
%%%

% loadexp;

%%% Write to data file
ncfname = fullfile('products',[expname,'_SSHvar.nc']);
XC = ncread(ncfname,'XC');
YC = ncread(ncfname,'YC');
bathy = ncread(ncfname,'bathy');
ssh_var = ncread(ncfname,'ssh_var');
ish_var = ncread(ncfname,'ish_var');



%%% Make plots

bathy_plot = bathy;
bathy_plot(sum(hFacC,3)==0) = NaN;
figure(1);
clf;
pcolor(XC,YC,sqrt(abs(ssh_var)));
hold on;
[C,h] = contour(XC,YC,-bathy_plot,[0 250 500 1000 2000 3000 4000],'EdgeColor',[.3 .3 .3]); 
clabel(C,h);
shading interp;
colormap(cmocean('amp',16));
colorbar;
caxis([0 0.08]);
xlabel('Longitude');
ylabel('Latitude');
title('SSH standard deviation (m)');
set(gca,'FontSize',14);


figure(2);
pcolor(XC,YC,sqrt(abs(ish_var)));
shading interp;
colormap(cmocean('amp',50));
colorbar;
caxis([0 1]);
