% %%%
% %%% paper3_plotTSfluxes.m
% %%%
% %%% Plots heat/salt fluxes in quasi-latitude coordinates
% %%%
% 
% %%% Load experiment
% expdir = '../experiments';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% loadexp;
% 
% %%% Reference surface freezing temperature
% theta0 = -1.9;
% 
% %%% Reference salinity
% salt0 = 34.6;
% % salt0 = 34.72;
% 
% %%% Set true to deform coordinates in the cavity
% deform_cavity = false;
% 
% %%% Set true to use barotropic streamfunction as the coordinate system
% use_PsiBT = false;
% 
% %%% Set true to use grounding line coordinate
% gl_coord = true;
% 
% %%% Set true to use depth-averaged temperature as the coordinate system
% use_meanT = false;
% 
% %%% Set true to decompose eddy fluxes
% calc_eddy_decomp = false;
% 
% %%% Parameters
% rho0 = 1000;
% Cp = 4000;
% 
% %%% Load HeatFunction data file
% outfname = [expname,'_HeatFunction'];
% if (use_PsiBT)
%   outfname = [outfname,'_PsiBT'];
% else 
%   if (use_meanT)
%     outfname = [outfname,'_meanT'];
%   else 
%     if (deform_cavity)
%       outfname = [outfname,'_deform'];
%     elseif (gl_coord)
%       outfname = [outfname,'_GLcoord'];
%     end
%   end
% end
% load(fullfile('products',outfname));
% Neta = length(eta);
% 
% %%% Load positive and negative heatfunction components
% outfname = [expname,'_PosNegHeatFunction'];
% if (use_PsiBT)
%   outfname = [outfname,'_PsiBT'];
% else
%   if (use_meanT)
%     outfname = [outfname,'_meanT'];
%   else 
%     if (deform_cavity)
%       outfname = [outfname,'_deform'];
%     elseif (gl_coord)
%       outfname = [outfname,'_GLcoord'];
%     end
%   end
% end
% outfname = [outfname,'.mat'];
% load(fullfile('products',outfname));
% 
% %%% Mask for ice/land
% psiT_tot_mean = mean(psiT_tot-psi_tot*theta0,3)* rho0*Cp/1e12;
% msk = ones(size(psiT_tot_mean));
% msk_ice = NaN*msk;
% for j=1:Neta  
%   idx = find(psiT_tot_mean(j,:)==psiT_tot_mean(j,1));
%   idx(end) = [];
%   msk(j,idx) = NaN;
%   if (~isempty(idx))
%     msk_ice(j,1:idx(end)) = 1;
%   end
%   idx = find(abs(psiT_tot_mean(j,:))<1e-12,1,'first');
%   msk(j,idx+1:end) = NaN;
% end
% 
% 
% psiT_tot_plot = -mean(psiT_tot-psi_tot*theta0,3) * rho0*Cp/1e12 .* msk;
% psiT_mean_plot = -mean(psiT_mean-psi_tot*theta0,3) * rho0*Cp/1e12 .* msk;
% 
% % msk = ones(size(psiT_neg_mean));
% % msk_ice = NaN*msk;
% % for j=1:Neta  
% %   idx = find(psiT_neg_mean(j,:)==psiT_neg_mean(j,1));
% %   idx(end) = [];
% %   msk(j,idx) = NaN;
% %   if (~isempty(idx))
% %     msk_ice(j,1:idx(end)) = 1;
% %   end
% %   idx = find(abs(psiT_neg_mean(j,:))<1e-12,1,'first');
% %   msk(j,idx+1:end) = NaN;
% % end
% 
% 
% psiT_pos_plot = -mean(psiT_pos_mean,3) * rho0*Cp/1e12 .* msk;
% psiT_neg_plot = -mean(psiT_neg_mean,3) * rho0*Cp/1e12 .* msk;
% 


%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 18;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.07 0.56 .86 .38];
axpos(2,:) = [0.07 0.06 .42 .38];
axpos(3,:) = [0.55 0.06 .42 .38];
cbpos = [0.955 0.55 0.01 .38];
axlabels = {'\textbf{A}','\textbf{B}','\textbf{C}'};
rho0 = 1027;
Cp = 4000;
arrowcolor = [0.301960784313725 0.35098039215686 0.93333333333333];
colororder = get(gca,'ColorOrder');
linewidth = 1.5;
ylim = [0 2000];
psimax = 4;
psistep = 0.4;
xlim = [-8.6 3.7];
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
hh = clabel(C,h,'FontSize',fontsize-4,'LabelSpacing',200);
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
hh = clabel(C,h,'FontSize',fontsize-4,'LabelSpacing',80);
% set(hh,'fontsize',fontsize-4,'Color',[.05 .05 .05],'BackgroundColor','none','Edgecolor','none')     
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
hh = clabel(C,h,'FontSize',fontsize-4,'LabelSpacing',200);
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
annotation('textbox',[axpos(1,1)-0.045 axpos(1,2)-0.055 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.045 axpos(2,2)-0.055 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.045 axpos(3,2)-0.055 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


figure1 = gcf;




% Create arrow
annotation(figure1,'arrow',[0.647 0.636],...
  [0.281857451403889 0.305615550755941],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.749 0.738],...
  [0.343412526997841 0.367170626349894],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.841000000000001 0.815000000000002],...
  [0.352051835853132 0.33261339092873],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.822 0.793000000000002],...
  [0.371490280777538 0.355291576673867],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.463 0.439000000000001],...
  [0.407127429805616 0.437365010799139],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.429 0.408000000000002],...
  [0.408207343412527 0.438444924406052],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.4 0.378000000000001],...
  [0.401727861771058 0.434125269978406],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.869 0.847000000000004],...
  [0.910367170626351 0.933045356371504],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.795 0.745],...
  [0.913606911447085 0.933045356371492],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.798 0.742],...
  [0.868250539956805 0.87041036717063],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.489 0.473000000000002],...
  [0.845572354211664 0.872570194384459],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.253 0.240000000000001],...
  [0.772138228941686 0.802375809935213],...
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
annotation(figure1,'arrow',[0.898000000000002 0.851],...
  [0.881209503239744 0.897408207343413],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.73 0.705000000000001],...
  [0.308855291576675 0.306695464362856],...
  'Color',[0.301960784313725 0.35098039215686 0.93333333333333],...
  'LineWidth',2,...
  'HeadStyle','plain');




% Create textbox
annotation(figure1,'textbox',...
  [0.0770000000000001 0.916904967602593 0.1 0.0242980561555075],...
  'Interpreter','latex',...
  'FontSize',fontsize,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'\textbf{Total heatfunction}'});

% Create textbox
annotation(figure1,'textbox',...
  [0.0740000000000001 0.412304535637151 0.1125 0.0242980561555076],...
  'Interpreter','latex',...
  'FontSize',fontsize,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'\textbf{``Warm'' heatfunction}'});

% Create textbox
annotation(figure1,'textbox',...
  [0.552 0.410144708423329 0.107 0.0242980561555076],...
  'Interpreter','latex',...
  'FontSize',fontsize,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'\textbf{``Cold'' heatfunction}'});

