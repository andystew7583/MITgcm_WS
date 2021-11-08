%%%
%%% paper2_plotMOCDecompDens.m
%%%
%%% Plots the m overturning circulation in density space.
%%%


%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%

% %%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_RTOPO2';
loadexp;
expname_notides = 'hires_seq_onetwelfth_notides_RTOPO2';

%%% Options (see calcOverturning)
calc_psi_eddy = true;
deform_cavity = false;
use_layers = true;
densvar = 'PD0';
psimax = 2;
psistep = 0.05;
% psimax = 2;
% psistep = 0.1;
ylim = [27.3 28.2];
% ylim = [27.3 29];

%%% Construct output file name
outfname = ['','_MOC_',densvar];
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
outfname_tides = [expname,outfname];
outfname_notides = [expname_notides,outfname];

%%% Load MOC data file
clear('psi_mean');
load(fullfile('products',outfname_tides));
psi_mean_tides = mean(psi_mean,3)/1e6;
psi_eddy_tides = mean(psi_eddy,3)/1e6;
psi_tot_tides = mean(psi_mean+psi_eddy,3)/1e6;
clear('psi_mean');
load(fullfile('products',outfname_notides));
psi_mean_plot = mean(psi_mean_stand,3)/1e6;
psi_fluc_plot = mean(psi_mean_fluc,3)/1e6;
psi_eddy_plot = mean(psi_eddy,3)/1e6;
psi_mean_notides = mean(psi_mean,3)/1e6;
psi_eddy_notides = mean(psi_eddy,3)/1e6;
psi_tot_notides = mean(psi_mean+psi_eddy,3)/1e6;
psi_tide_plot = psi_tot_tides-psi_tot_notides;





%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.06 0.54 .4 .41];
axpos(2,:) = [0.53 0.54 .4 .41];
axpos(3,:) = [0.06 0.05 .4 .41];
axpos(4,:) = [0.53 0.05 .4 .41];
cbpos = [0.95 0.05 0.015 .91];
axlabels = {'(a)','(b)','(c)','(d)'};








%%% Set up the figure
figure(209)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[382          55        1120         930]);




%%% Mean streamfunction
subplot('Position',axpos(1,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_mean_plot);
shading interp;
hold on;
[C,h] = contour(EE,DD,psi_mean_plot,[-1 -0.5 -0.25 -0.1 0.1 0.25 0.5 1],'EdgeColor','k');
clabel(C,h);
[C,h] = contour(EE,DD,psi_mean_plot,[-6 -5 -4 -3 -2],'EdgeColor','w');
clabel(C,h,'Color','w');
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
hold off;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 10.5]);
colormap(gca,cmocean('balance',round(2*psimax/(psistep))));%,'pivot',0));
% xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)')
set(gca,'FontSize',fontsize);
text(-6.5,27.35,'Mean component','FontSize',fontsize+2);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Top');
set(ax2,'YLim',get(ax1,'YLim'));
set(ax2,'YTick',[]);
box off;
set(ax2,'XLim',get(ax1,'XLim')-78.16);
set(ax2,'FontSize',fontsize);
set(get(ax2,'XLabel'),'String','Pseudo-Latitude');

%%% Fluctuating streamfunction
subplot('Position',axpos(2,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_fluc_plot);
shading interp;
hold on;
[C,h] = contour(EE,DD,psi_fluc_plot,[-6 -5 -4 -3 -2 -1 -0.5 -0.25 -0.1 0.1 0.25 0.5 1],'EdgeColor','k');
clabel(C,h);
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
hold off;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 10.5]);
colormap(gca,cmocean('balance',round(2*psimax/(psistep))));%,'pivot',0));
% xlabel('MOC coordinate, \eta');
% ylabel('Potential density (kg/m^3)')
set(gca,'FontSize',fontsize);
text(-6.5,27.35,'Seasonal/interannual component','FontSize',fontsize+2);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Top');
set(ax2,'YLim',get(ax1,'YLim'));
set(ax2,'YTick',[]);
box off;
set(ax2,'XLim',get(ax1,'XLim')-78.16);
set(ax2,'FontSize',fontsize);
set(get(ax2,'XLabel'),'String','Pseudo-Latitude');

%%% Eddy streamfunction
subplot('Position',axpos(3,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_eddy_plot);
shading interp;
hold on;
[C,h] = contour(EE,DD,psi_eddy_plot,[-6 -5 -4 -3 -2 -1 -0.5 -0.25 -0.1 0.1 0.25 0.5 1],'EdgeColor','k');
clabel(C,h);
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
hold off;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 10.5]);
colormap(gca,cmocean('balance',round(2*psimax/(psistep))));%,'pivot',0));
xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)')
set(gca,'FontSize',fontsize);
text(-6.5,27.35,'Eddy component','FontSize',fontsize+2);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Top');
set(ax2,'YLim',get(ax1,'YLim'));
set(ax2,'YTick',[]);
box off;
set(ax2,'XLim',get(ax1,'XLim')-78.16);
set(ax2,'FontSize',fontsize);
% set(get(ax2,'XLabel'),'String','Pseudo-Latitude');

%%% Tidal streamfunction
subplot('Position',axpos(4,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_tide_plot);
shading interp;
hold on;
[C,h] = contour(EE,DD,psi_tide_plot,[-6 -5 -4 -3 -2 -1 -0.5 -0.25 -0.1 0.1 0.25 0.5 1],'EdgeColor','k');
clabel(C,h);
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
hold off;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 10.5]);
colormap(gca,cmocean('balance',round(2*psimax/(psistep))));%,'pivot',0));
xlabel('MOC coordinate, \eta');
set(gca,'FontSize',fontsize);
text(-6.5,27.35,'Tides - No Tides','FontSize',fontsize+2);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);

%%% Add colorbar
cbhandle = colorbar;
set(cbhandle,'Position',cbpos);
title(cbhandle,'Sv');

%%% Finish last plot
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Top');
set(ax2,'YLim',get(ax1,'YLim'));
set(ax2,'YTick',[]);
box off;
set(ax2,'XLim',get(ax1,'XLim')-78.16);
set(ax2,'FontSize',fontsize);
% set(get(ax2,'XLabel'),'String','Pseudo-Latitude');



%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.04 axpos(4,2)-0.05 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


[j1,k1]=find(psi_eddy_plot==min(psi_eddy_plot(:)))
% j1 = 131
[k2] = find(psi_tot_notides(j1,:)==min(psi_tot_notides(j1,:)))
total_AABW_flux = -psi_tot_notides(j1,k2)
total_eddy_flux_below_psi_max = -psi_eddy_plot(j1,k2)
total_AABW_flux_below_eddy_max = -psi_tot_notides(j1,k1)
total_eddy_flux = -psi_eddy_plot(j1,k1)