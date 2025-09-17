%%%
%%% paper3_plotTSfluxes.m
%%%
%%% Plots heat/salt fluxes in quasi-latitude coordinates
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;

%%% Reference surface freezing temperature
theta0 = -1.9;

%%% Reference salinity
salt0 = 34.6;
% salt0 = 34.72;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to use grounding line coordinate
gl_coord = true;

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
    elseif (gl_coord)
      outfname = [outfname,'_GLcoord'];
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
    elseif (gl_coord)
      outfname = [outfname,'_GLcoord'];
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

% msk = ones(size(psiT_neg_mean));
% msk_ice = NaN*msk;
% for j=1:Neta  
%   idx = find(psiT_neg_mean(j,:)==psiT_neg_mean(j,1));
%   idx(end) = [];
%   msk(j,idx) = NaN;
%   if (~isempty(idx))
%     msk_ice(j,1:idx(end)) = 1;
%   end
%   idx = find(abs(psiT_neg_mean(j,:))<1e-12,1,'first');
%   msk(j,idx+1:end) = NaN;
% end


psiT_pos_plot = -mean(psiT_pos_mean,3) * rho0*Cp/1e12 .* msk;
psiT_neg_plot = -mean(psiT_neg_mean,3) * rho0*Cp/1e12 .* msk;



%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.07 0.55 .86 .4];
axpos(2,:) = [0.07 0.05 .42 .4];
axpos(3,:) = [0.55 0.05 .42 .4];
cbpos = [0.96 0.55 0.01 .42];
axlabels = {'\textbf{A}','\textbf{B}','\textbf{C}'};
rho0 = 1027;
Cp = 4000;
arrowcolor = [0.301960784313725 0.35098039215686 0.93333333333333];
colororder = get(gca,'ColorOrder');
linewidth = 1.5;
ylim = [0 2000];
psimax = 4;
psistep = 0.2;
xlim = [-8.6 4];
icecolor = [186 242 239]/255;
psisteps = [0:psistep:-psistep psistep:psistep:psimax];
[ZZ,EE] = meshgrid(squeeze(-RF),eta);






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
colormap(gca,cmocean('amp',length(psisteps)));
cbhandle = colorbar;
set(cbhandle,'Position',cbpos);
title(cbhandle,'TW','FontSize',fontsize);
xlabel('Cross-shelf coordinate, \eta');
ylabel('Depth (m)');
% title('Total heat function');
set(gca,'FontSize',fontsize);
caxis([0 psimax]);

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

%%% Add reference latitudes
ax3 = axes('Position',get(ax1,'Position'));
set(ax3,'Color','None');
set(ax3,'XAxisLocation','Top');
set(ax3,'YAxisLocation','Right');
set(ax3,'YLim',get(ax1,'YLim'));
set(ax3,'YTick',[]);
box off;
set(ax3,'XLim',xlim-77);
set(ax3,'FontSize',fontsize);
set(get(ax3,'XLabel'),'String','Reference latitude');


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
colormap(gca,cmocean('amp',length(psisteps)));
xlabel('Cross-shelf coordinate, \eta');
ylabel('Depth (m)');
% title('``Warm'''' heat function');
set(gca,'FontSize',fontsize);
caxis([0 psimax]);

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

%%% Add reference latitudes
ax3 = axes('Position',get(ax1,'Position'));
set(ax3,'Color','None');
set(ax3,'XAxisLocation','Top');
set(ax3,'YAxisLocation','Right');
set(ax3,'YLim',get(ax1,'YLim'));
set(ax3,'YTick',[]);
box off;
set(ax3,'XLim',xlim-77);
set(ax3,'FontSize',fontsize);
set(get(ax3,'XLabel'),'String','Reference latitude');

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
colormap(gca,cmocean('amp',length(psisteps)));
xlabel('Cross-shelf coordinate, \eta');
% ylabel('Depth (m)');
% title('``Cold'''' heat function');
set(gca,'FontSize',fontsize);
caxis([0 psimax]);

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

%%% Add reference latitudes
ax3 = axes('Position',get(ax1,'Position'));
set(ax3,'Color','None');
set(ax3,'XAxisLocation','Top');
set(ax3,'YAxisLocation','Right');
set(ax3,'YLim',get(ax1,'YLim'));
set(ax3,'YTick',[]);
box off;
set(ax3,'XLim',xlim-77);
set(ax3,'FontSize',fontsize);
set(get(ax3,'XLabel'),'String','Reference latitude');


%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


figure1 = gcf;



% Create arrow
annotation(figure1,'arrow',[0.525 0.509000000000002],...
  [0.861771058315335 0.88876889848813],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.422 0.407000000000002],...
  [0.845572354211663 0.875809935205192],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.289 0.276000000000001],...
  [0.831533477321815 0.861771058315342],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.439 0.387],...
  [0.809935205183586 0.807775377969766],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.582 0.522000000000001],...
  [0.842332613390929 0.82397408207344],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.791 0.735],...
  [0.861771058315336 0.863930885529161],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.898000000000002 0.851],...
  [0.881209503239744 0.897408207343413],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.861 0.839000000000004],...
  [0.925485961123111 0.948164146868264],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.763 0.713],...
  [0.924406047516199 0.943844492440606],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.381 0.366000000000001],...
  [0.408207343412527 0.442764578833697],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.424 0.403000000000002],...
  [0.414686825053996 0.444924406047521],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.461 0.437000000000001],...
  [0.416846652267819 0.447084233261342],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.73 0.705000000000001],...
  [0.308855291576675 0.306695464362856],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.726 0.715],...
  [0.352051835853132 0.375809935205184],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.77 0.759],...
  [0.360691144708424 0.384449244060477],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.824000000000001 0.796000000000001],...
  [0.39416846652268 0.38012958963283],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.832000000000001 0.806000000000002],...
  [0.360691144708424 0.341252699784022],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');




% Create textbox
annotation(figure1,'textbox',...
  [0.0770000000000001 0.916904967602593 0.1 0.0242980561555075],...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'\textbf{Total heatfunction}'});

% Create textbox
annotation(figure1,'textbox',...
  [0.0740000000000001 0.422304535637151 0.1125 0.0242980561555076],...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'\textbf{``Warm'' heatfunction}'});

% Create textbox
annotation(figure1,'textbox',...
  [0.552 0.420144708423329 0.107 0.0242980561555076],...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'\textbf{``Cold'' heatfunction}'});

