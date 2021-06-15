%%%
%%% plotOverturning.m
%%%
%%% Plots the overturning circulation in density space.
%%%

addpath CDT/cdt;

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';

%%% Options (see calcOverturning)
calc_psi_eddy = true;
deform_cavity = false;
use_layers = true;
densvar = 'PD0';
psimax = 4;
% psistep = 0.5;
% psimax = 2;
psistep = 0.05;
ylim = [27.3 28.2];
% ylim = [27.3 29];

%%% Construct output file name
outfname = [expname,'_MOC_',densvar];
if (calc_psi_eddy)
  if (use_layers)
    estr = '_layers';
  else
    estr = '_TRM';
  end
else
  estr = '_noeddy';
end
outfname = [outfname,estr];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];


%%% Load MOC data file
load(fullfile('products',outfname));

psi_plot = mean(psi_mean+psi_eddy,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(1);
clf;
set(gcf,'Position',[325         460        1023         398]);
[DD,EE] = meshgrid(dens_levs,eta);
contourf(EE,DD,psi_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
% pcolor(EE,DD,psi_plot);
% shading interp;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
% colormap redblue(32);
colormap(cmocean('balance',round(2*psimax/psistep)));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Overturning streamfunction (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);

psi_plot = mean(psi_mean,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(2);
clf;
set(gcf,'Position',[325         460        1023         398]);
[DD,EE] = meshgrid(dens_levs,eta);
contourf(EE,DD,psi_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
% pcolor(EE,DD,psi_plot);
% shading interp;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
% colormap redblue(32);
colormap(cmocean('balance',round(2*psimax/psistep)));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Overturning streamfunction (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);




psi_plot = mean(psi_eddy,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(3);
clf;
set(gcf,'Position',[325         460        1023         398]);
[DD,EE] = meshgrid(dens_levs,eta);
contourf(EE,DD,psi_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
% pcolor(EE,DD,psi_plot);
% shading interp;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
% set(gca,'YLim',[36.8 37.5]);
set(gca,'YLim',ylim);
% colormap redblue(32);
colormap(cmocean('balance',round(2*psimax/psistep)));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Overturning streamfunction (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);

psi_plot = mean(psi_mean_stand,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(4);
clf;
set(gcf,'Position',[325         460        1023         398]);
[DD,EE] = meshgrid(dens_levs,eta);
contourf(EE,DD,psi_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
% pcolor(EE,DD,psi_plot);
% shading interp;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
% colormap redblue(32);
colormap(cmocean('balance',round(2*psimax/psistep)));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Overturning streamfunction (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);

psi_plot = mean(psi_mean_fluc,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(5);
clf;
set(gcf,'Position',[325         460        1023         398]);
[DD,EE] = meshgrid(dens_levs,eta);
contourf(EE,DD,psi_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
% pcolor(EE,DD,psi_plot);
% shading interp;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
% colormap redblue(32);
colormap(cmocean('balance',round(2*psimax/psistep)));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Overturning streamfunction (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);


Nt = length(times);
Neta = length(eta);
Nd = length(dens_levs);
psi_tot = psi_mean+psi_eddy;
% psi_tot = psi_tot - repmat(squeeze(mean(reshape(psi_tot,[Neta Nd 12,size(psi_tot,3)/12]),4)),[1 1 size(psi_tot,3)/12]);
[eof_maps,pc,expvar] = eof((psi_tot)/1e6);


psi_plot = std((psi_mean+psi_eddy)/1e6,[],3);
% psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(6);
clf;
set(gcf,'Position',[325         460        1023         398]);
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
% caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
% colormap redblue(32);
colormap hot;
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Overturning streamfunction (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);

psi_plot = eof_maps(:,:,end)*std(pc(end,:));
% psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(7);
clf;
set(gcf,'Position',[325         460        1023         398]);
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
% caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
% colormap redblue(32);
colormap(cmocean('balance'));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Overturning streamfunction (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);

figure(8);
plot(times/t1year,pc(end,:)./std(pc(end,:)));