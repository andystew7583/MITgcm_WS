%%%
%%% plotOverturning.m
%%%
%%% Plots the overturning circulation in density space.
%%%

load([expname,'_MOC_dens.mat']);


figure(2);
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,mean(psi_mean,3)/1e6);
% pcolor(EE,DD,mean(psi_mean+psi_eddy,3)/1e6);
% pcolor(EE,DD,mean(psi_eddy,3)/1e6);
shading interp;
caxis([-4 4]);
set(gca,'YDir','reverse');
% set(gca,'YLim',[36.8 37.5]);
set(gca,'YLim',[27.3 28.2]);
colormap redblue(32);
colorbar;