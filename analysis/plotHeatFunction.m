%%%
%%% plotHeatFunction.m
%%%
%%% Plots the heat function in eta/z space.
%%%

addpath CDT/cdt;

%%% Load experiment
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_RTOPO2';
% expname = 'hires_seq_onetwelfth_RTOPO2';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% loadexp;

%%% Reference surface freezing temperature
theta0 = -1.9;

%%% Reference salinity
salt0 = 34.6;
% salt0 = 34.72;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;


%%% Set true to use depth-averaged temperature as the coordinate system
use_meanT = false;

%%% Set true to decompose eddy fluxes
calc_eddy_decomp = false;

%%% Parameters
rho0 = 1000;
Cp = 4000;
ylim = [0 2000];
psimax = 4;
psistep = 0.2;
fontsize = 14;
xlim = [-2.6 0];
% xlim = [-9 4];

%%% Load HeatFunction data file
outfname = [expname,'_HeatFunction'];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else 
  if (use_meanT)
    outfname = [outfname,'_meanT'];
  else 
    if (deform_cavity)
      outfname = [outfname,'_deform'];
    end
  end
end
load(fullfile('products',outfname));
Neta = length(eta);

%%% Some quantities computed separately for higher-resolution experiments
if (strcmp(expname,'hires_seq_onetwentyfourth_notides_RTOPO2'))

  %%% Load eddy-induced components of heatfunction/transport
  outfname = [expname,'_HeatFunctionEddyDecomp'];
  if (use_PsiBT)
    outfname = [outfname,'_PsiBT'];
  else
    if (use_meanT)
      outfname = [outfname,'_meanT'];
    else 
      if (deform_cavity)
        outfname = [outfname,'_deform'];
      end
    end
  end
  outfname = [outfname,'.mat'];
  load(fullfile('products',outfname));

  psiT_eddy_adv = psiT_eddy_adv - psi_eddy*theta0;
  psiT_eddy_stir = psiT_eddy - psiT_eddy_adv;
  
  %%% Load positive and negative heatfunction components
  outfname = [expname,'_PosNegHeatFunction'];
  if (use_PsiBT)
    outfname = [outfname,'_PsiBT'];
  else
    if (use_meanT)
      outfname = [outfname,'_meanT'];
    else 
      if (deform_cavity)
        outfname = [outfname,'_deform'];
      end
    end
  end
  outfname = [outfname,'.mat'];
  load(fullfile('products',outfname));

end

%%% Mask for ice/land
psiT_tot_mean = mean(psiT_tot-psi_tot*theta0,3)* rho0*Cp/1e12;
msk = ones(size(psiT_tot_mean));
for j=1:Neta  
  idx = find(psiT_tot_mean(j,:)==psiT_tot_mean(j,1));
  idx(end) = [];
  msk(j,idx) = NaN;
  idx = find(abs(psiT_tot_mean(j,:))<1e-12,1,'first');
  msk(j,idx+1:end) = NaN;
end

psiT_plot = -mean(psiT_tot-psi_tot*theta0,3) * rho0*Cp/1e12;
psiT_plot = psiT_plot.*msk;
figure(81);
clf;
set(gcf,'Position',[325         460        1023         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
[C,h] = contourf(EE,ZZ,psiT_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(cmocean('balance',40));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Heat function (TW)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);
caxis([-psimax psimax]);

% axpos = [0.0887    0.1300    0.8283    0.7950];
axpos = [0.1487    0.1700    0.7083    0.7250];

psiT_plot = -mean(psiT_mean-psi_tot*theta0,3) * rho0*Cp/1e12;
psiT_plot = psiT_plot.*msk;
figure(82);
clf;
set(gcf,'Position',[325         460        550         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
[C,h] = contourf(EE,ZZ,psiT_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(cmocean('balance',40));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Heat function, mean component (TW)');
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);




psiT_plot = -mean(psiT_pos_mean,3) * rho0*Cp/1e12;
psiT_plot = psiT_plot.*msk;
figure(84);
clf;
set(gcf,'Position',[325         460        550         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
[C,h] = contourf(EE,ZZ,psiT_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(cmocean('balance',40));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Heat function, positive component (TW)');
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);


psiT_plot = -mean(psiT_neg_mean,3) * rho0*Cp/1e12;
psiT_plot = psiT_plot.*msk;
figure(85);
clf;
set(gcf,'Position',[325         460        550         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
[C,h] = contourf(EE,ZZ,psiT_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(cmocean('balance',40));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Heat function, negative component (TW)');
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);

psiT_plot = -mean(psiT_eddy,3) * rho0*Cp/1e12;
psiT_plot = psiT_plot.*msk;
figure(86);
clf;
set(gcf,'Position',[325         460        550         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
[C,h] = contourf(EE,ZZ,psiT_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
colormap(cmocean('balance',40));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Heat function, eddy component (TW)');
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);


psiT_plot = -mean(psiT_eddy_adv,3) * rho0*Cp/1e12;
psiT_plot = psiT_plot.*msk;
figure(87);
clf;
set(gcf,'Position',[325         460        550         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
[C,h] = contourf(EE,ZZ,psiT_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
colormap(cmocean('balance',40));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Heat function, eddy advection component (TW)');
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);

psiT_plot = -mean(psiT_eddy_stir,3) * rho0*Cp/1e12;
psiT_plot = psiT_plot.*msk;
figure(88);
clf;
set(gcf,'Position',[325         460        550         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
[C,h] = contourf(EE,ZZ,psiT_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
colormap(cmocean('balance',40));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Heat function, eddy stirring component (TW)');
set(gca,'Position',axpos);
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);



psimax = 1;
psistep = 0.05;

psiS_plot = -mean(psiS_tot-psi_tot*salt0,3) * rho0/1e9;
psiS_plot = psiS_plot.*msk;
figure(89);
clf;
set(gcf,'Position',[425         460        1023         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
contourf(EE,ZZ,psiS_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
set(gca,'XLim',xlim);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7]);
colormap(cmocean('balance',40));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Salt function (Gg/s)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);


psiS_plot = -mean(psiS_mean-psi_tot*salt0,3) * rho0/1e9;
psiS_plot = psiS_plot.*msk;
figure(90);
clf;
set(gcf,'Position',[425         460        1023         398]);
[ZZ,EE] = meshgrid(squeeze(-RF),eta);
contourf(EE,ZZ,psiS_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
set(gca,'XLim',xlim);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7]);
colormap(cmocean('balance',40));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Salt function, mean component (Gg/s)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);



psi_plot = mean(psi_tot,3);
figure(90);
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


psiT_mod = psiT_tot-psi_tot*theta0;
psiT_mean_mod = psiT_mean-psi_tot*theta0;
 
hflux_int = -mean(psiT_mod(:,1,:),3) * rho0*Cp/1e12;
hflux_mean_int = -mean(psiT_mean_mod(:,1,:),3)* rho0*Cp/1e12;
hflux_eddy_int = -mean(psiT_eddy(:,1,:),3)* rho0*Cp/1e12;
hflux_eddy_adv_int = -mean(psiT_eddy_adv(:,1,:),3)* rho0*Cp/1e12;
hflux_eddy_stir_int = -mean(psiT_eddy_stir(:,1,:),3)* rho0*Cp/1e12;

figure(5);
plot(eta,hflux_int,'LineWidth',1.5);
hold on;
plot(eta,hflux_mean_int,'LineWidth',1.5);
plot(eta,hflux_eddy_int,'LineWidth',1.5);
plot(eta,hflux_eddy_adv_int,'LineWidth',1.5);
plot(eta,hflux_eddy_stir_int,'LineWidth',1.5);
hold off;
set(gca,'XLim',xlim);
xlabel('MOC coordinate \eta');
ylabel('Heat flux (TW)');
set(gca,'FontSize',14);
legend('Total','Mean','Eddy','Eddy advection','Eddy stirring');

figure(6);
[C,h]=contourf(XC,YC,ETA,eta,'EdgeColor','k');
clabel(C,h);
colormap jet;
xlabel('Longitude');
ylabel('Latitude');
title('Heatfunction coordinate');