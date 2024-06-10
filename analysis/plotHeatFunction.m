%%%
%%% plotHeatFunction.m
%%%
%%% Plots the heat function in eta/z space.
%%%

addpath CDT/cdt;

%%% Load experiment
expdir = '../experiments';
% expname = 'hires_seq_onethird_RTOPO2';
expname = 'hires_seq_onesixth_RTOPO2';
% expname = 'hires_seq_onetwelfth_RTOPO2';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% loadexp;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to decompose eddy fluxes
calc_eddy_decomp = false;

%%% Parameters
rho0 = 1000;
Cp = 4000;
ylim = [0 2000];

%%% Construct output file name
outfname = [expname,'_HeatFunction'];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else
  if (deform_cavity)
    outfname = [outfname,'_deform'];
  end
end
outfname = [outfname,'.mat'];


%%% Load MOC data file
load(fullfile('products',outfname));

theta0 = -1.9;
% theta0 = -2.1;
psiT_plot = mean(psiT_tot-psi_tot*theta0,3) * rho0*Cp/1e12;
for j=1:Neta
  idx = find(psiT_plot(j,:)==psiT_plot(j,1));
  idx(end) = [];
  psiT_plot(j,idx) = NaN;
  idx = find(abs(psiT_plot(j,2:end))<1e-12,1,'first');
  psiT_plot(j,idx+1:end) = NaN;
end
figure(1);
clf;
set(gcf,'Position',[325         460        1023         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
% contourf(EE,DD,psi_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
pcolor(EE,ZZ,psiT_plot); 
shading interp;
% caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(cmocean('balance',40));
colorbar;
xlabel('\eta');
ylabel('Depth (m)');
title('Heat function (TW)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);
caxis([-5 5])

S0 = 34.72;
psiS_plot = mean(psiS_tot-psi_tot*S0,3) * rho0/1e9;
figure(2);
clf;
set(gcf,'Position',[425         460        1023         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
% contourf(EE,DD,psi_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
pcolor(EE,ZZ,psiS_plot); 
shading interp;
% caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(cmocean('balance',40));
colorbar;
xlabel('\eta');
ylabel('Depth (m)');
title('Salt function (Gg/s)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);
caxis([-1 1])


psi_plot = mean(psi_tot,3);
figure(3);
clf;
set(gcf,'Position',[525         460        1023         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
pcolor(EE,ZZ,psi_plot); 
shading interp;
% caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(cmocean('balance',30));
colorbar;
xlabel('\eta');
ylabel('Depth (m)');
title('Eulerian-mean streamfunction (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);
% caxis([-20 20])


thflux_tot = psiT_tot(:,1,:)*rho0*Cp/1e12;
figure(4);

plot(eta,mean(thflux_tot,3));
hold on;
plot(eta,std(thflux_tot,1,3));
hold off;


thflux_eddy = psiT_eddy(:,2,:)*rho0*Cp/1e12;
figure(5);
plot(eta,mean(thflux_eddy,3));
hold on;
plot(eta,std(thflux_eddy,1,3));
hold off;


psiT_mod = psiT_tot-psi_tot*theta0;
eidx_icefront = 80;
eidx_shelfbreak = 126;
zidx_icefront = 26;
hflux_icefront = squeeze(psiT_mod(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_shelfbreak = squeeze(psiT_mod(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

hflux_eddy_icefront = squeeze(psiT_eddy(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_eddy_shelfbreak = squeeze(psiT_eddy(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

figure(6);
plot(times/t1year,hflux_icefront);
hold on;
plot(times/t1year,hflux_shelfbreak);

plot(times/t1year,hflux_eddy_icefront);

plot(times/t1year,hflux_eddy_shelfbreak);
hold off