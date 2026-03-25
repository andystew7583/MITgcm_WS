%%%
%%% compareOverturning.m
%%%
%%% Compares tidal/tide-free overturning circulation in density space.
%%%

%%% Load experiment
expdir = '../experiments';
expname_tides = 'hires_seq_onethird_RTOPO2';
expname_notides = 'hires_seq_onethird_notides_RTOPO2';

%%% Options (see calcOverturning)
calc_psi_eddy = true;
deform_cavity = false;
use_layers = true;
densvar = 'PD0';
% psimax = 6;
% psistep = 0.5;
psimax = 1;
psistep = 0.05;

%%% Construct output file name
outfname_base = ['_MOC_',densvar];
if (calc_psi_eddy)
  if (use_layers)
    estr = '_layers';
  else
    estr = '_TRM';
  end
else
  estr = '_noeddy';
end
outfname_base = [outfname_base,estr];
if (deform_cavity)
  outfname_base = [outfname_base,'_deform'];
end
outfname_base = [outfname_base,'.mat'];
outfname_tides = [expname_tides,outfname_base];
outfname_notides = [expname_notides,outfname_base];



%%% Load MOC data file
load(fullfile('products',outfname_tides));
psi_mean_tides = psi_mean;
psi_eddy_tides = psi_eddy;
psi_total_tides = psi_mean_tides+psi_eddy_tides;
load(fullfile('products',outfname_notides));
psi_mean_notides = psi_mean;
psi_eddy_notides = psi_eddy;
psi_total_notides = psi_mean_notides+psi_eddy_notides;




%%% Difference in total MOC
psi_plot = mean(psi_total_tides-psi_total_notides,3)/1e6;
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
% set(gca,'YLim',[36.8 37.5]);
set(gca,'YLim',[27.3 28.2]);
% colormap redblue(32);
colormap(cmocean('balance',round(2*psimax/psistep)));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Total MOC, tides - noTides (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);


%%% Difference in mean MOC
psi_plot = mean(psi_mean_tides-psi_mean_notides,3)/1e6;
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
% set(gca,'YLim',[36.8 37.5]);
set(gca,'YLim',[27.3 28.2]);
% colormap redblue(32);
colormap(cmocean('balance',round(2*psimax/psistep)));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Mean MOC, tides - noTides (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);


%%% Difference in eddy MOC
psi_plot = mean(psi_eddy_tides-psi_eddy_notides,3)/1e6;
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
set(gca,'YLim',[27.3 28.2]);
% colormap redblue(32);
colormap(cmocean('balance',round(2*psimax/psistep)));
colorbar;
xlabel('\eta');
ylabel('Potential density (kg/m^3)');
title('Transient MOC, tides - noTides (Sv)');
set(gca,'Position',[0.0787    0.1300    0.8383    0.7950]);
set(gca,'FontSize',12);