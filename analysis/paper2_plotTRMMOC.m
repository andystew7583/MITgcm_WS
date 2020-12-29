%%%
%%% paper2_plotTRMMOC.m
%%%
%%% Compares eddy-induced overturning from LAYERS and TRM.
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_notides_RTOPO2';

%%% Options (see calcOverturning)
calc_psi_eddy = true;
deform_cavity = false;
densvar = 'PD0';
% psimax = 6;
% psistep = 0.5;
psimax = 2;
psistep = 0.1;


%%% Load MOC data file with LAYERS
use_layers = true;
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
load(fullfile('products',outfname));
psi_layers = mean(psi_eddy,3)/1e6;


%%% Load MOC data file with TRM
use_layers = false;
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
load(fullfile('products',outfname));
psi_trm = mean(psi_eddy,3)/1e6;





%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE THE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
fontsize = 14;
framepos = [325         460        923         398];
labelspacing = 200;
axpos = zeros(6,4);
axpos(1,:) = [0.08 0.12 0.39 0.76];
axpos(2,:) = [0.55 0.12 0.39 0.76];
cb1_pos = [0.95 0.12 0.015 0.76];
axlabels = {'(a)','(b)'};
ylim = [27.4 28.1];
xlim = [2.5 6];



figure(216);
clf;
set(gcf,'Position',framepos);   





subplot('Position',axpos(1,:));
[DD,EE] = meshgrid(dens_levs,eta);
contourf(EE,DD,psi_layers,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
colormap(cmocean('balance',round(2*psimax/psistep)));
xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)');
text(2.6,27.45,'Diagnosed','FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
box off;

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Top');
set(ax2,'YAxisLocation','Right');
set(ax2,'YLim',get(ax1,'YLim'));
set(ax2,'YTick',[]);
box off;
set(ax2,'XLim',get(ax1,'XLim')-78.16);
set(ax2,'FontSize',fontsize);
set(get(ax2,'XLabel'),'String','Pseudo-Latitude');

subplot('Position',axpos(2,:));
[DD,EE] = meshgrid(dens_levs,eta);
contourf(EE,DD,psi_trm,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
colormap(cmocean('balance',round(2*psimax/psistep)));
cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
title(cbhandle,'Sv');
xlabel('MOC coordinate, \eta');
% ylabel('Potential density (kg/m^3)');
set(gca,'FontSize',fontsize);
text(2.6,27.45,'NTRM','FontSize',fontsize+2);
box off;

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Top');
set(ax2,'YAxisLocation','Right');
set(ax2,'YLim',get(ax1,'YLim'));
set(ax2,'YTick',[]);
box off;
set(ax2,'XLim',get(ax1,'XLim')-78.16);
set(ax2,'FontSize',fontsize);
set(get(ax2,'XLabel'),'String','Pseudo-Latitude');

%%% Panel labels
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.07 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
