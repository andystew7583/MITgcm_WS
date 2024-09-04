%%%
%%% paper3_plotTSfluxes.m
%%%
%%% Plots heat/salt fluxes in quasi-latitude coordinates
%%%

%%% Load experiment
expdir = '../experiments';
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

%%% Mask for ice/land
psiT_tot_mean = mean(psiT_tot-psi_tot*theta0,3)* rho0*Cp/1e12;
msk = ones(size(psiT_tot_mean));
msk_ice = NaN*msk;
for j=1:Neta  
  idx = find(psiT_tot_mean(j,:)==psiT_tot_mean(j,1));
  idx(end) = [];
  msk(j,idx) = NaN;
  if (~isempty(idx))
    msk_ice(j,1:idx(end)) = 1;
  end
  idx = find(abs(psiT_tot_mean(j,:))<1e-12,1,'first');
  msk(j,idx+1:end) = NaN;
end


psiT_tot_plot = -mean(psiT_tot-psi_tot*theta0,3) * rho0*Cp/1e12 .* msk;
psiT_mean_plot = -mean(psiT_mean-psi_tot*theta0,3) * rho0*Cp/1e12 .* msk;
psiT_pos_plot = -mean(psiT_pos_mean,3) * rho0*Cp/1e12 .* msk;
psiT_neg_plot = -mean(psiT_neg_mean,3) * rho0*Cp/1e12 .* msk;



%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.05 0.55 .9 .4];
axpos(2,:) = [0.05 0.05 .4 .4];
axpos(3,:) = [0.55 0.05 .4 .4];
cbpos = [0.96 0.05 0.01 .85];
axlabels = {'(a)','(b)','(c)'};
rho0 = 1027;
Cp = 4000;
colororder = get(gca,'ColorOrder');
linewidth = 1.5;
ylim = [0 2000];
psimax = 4;
psistep = 0.2;
xlim = [-8.6 4];
icecolor = [186 242 239]/255;
psisteps = [-psimax:psistep:-psistep psistep:psistep:psimax];
[ZZ,EE] = meshgrid(squeeze(-RF),eta);

%%% Plotting range for salinity figure
latMin_b = min(min(YC));
latMax_b = YC(1,end-spongethickness);
lonMin_b = min(min(XC));
lonMax_b = XC(end-spongethickness,1);







%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(203)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34  1000  926]);





%%% Total heat function
axes('Position',axpos(1,:));
pcolor(EE,ZZ,psiT_tot_plot);
shading flat;
hold on;
[C,h] = contour(EE,ZZ,psiT_tot_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
hold off;
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(gca,cmocean('balance',length(psisteps)));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('Total heat function (TW)');
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
pcolor(EE,ZZ,msk_ice);
shading flat;
colormap(ax2,icecolor);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
set(ax2,'YLim',ylim);
set(ax2,'YDir','reverse');
set(ax2,'XLim',xlim);
set(ax2,'Color','None')

%%% Positive heat function
axes('Position',axpos(2,:));
pcolor(EE,ZZ,psiT_pos_plot);
shading flat;
hold on;
[C,h] = contour(EE,ZZ,psiT_pos_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
hold off;
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(gca,cmocean('balance',length(psisteps)));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('``Warm'''' component (TW)');
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
pcolor(EE,ZZ,msk_ice);
shading flat;
colormap(ax2,icecolor);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
set(ax2,'YLim',ylim);
set(ax2,'YDir','reverse');
set(ax2,'XLim',xlim);
set(ax2,'Color','None')

%%% Negative heat function
axes('Position',axpos(3,:));
pcolor(EE,ZZ,psiT_neg_plot);
shading flat;
hold on;
[C,h] = contour(EE,ZZ,psiT_neg_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
hold off;
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(gca,cmocean('balance',length(psisteps)));
colorbar;
xlabel('MOC coordinate \eta');
ylabel('Depth (m)');
title('``Cold'''' component (TW)');
set(gca,'FontSize',fontsize);
caxis([-psimax psimax]);

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
pcolor(EE,ZZ,msk_ice);
shading flat;
colormap(ax2,icecolor);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
set(ax2,'YLim',ylim);
set(ax2,'YDir','reverse');
set(ax2,'XLim',xlim);
set(ax2,'Color','None')


%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
