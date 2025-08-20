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

%%% Cross-shelf locations of ice front and shelf break for heat budget calculation
eta_icefront = -1.1;
eta_shelfbreak = 3.5;  

%%% Depth of ice front for heat budget calculation
zidx_icefront = 25;    

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

msk = ones(size(psiT_neg_mean));
msk_ice = NaN*msk;
for j=1:Neta  
  idx = find(psiT_neg_mean(j,:)==psiT_neg_mean(j,1));
  idx(end) = [];
  msk(j,idx) = NaN;
  if (~isempty(idx))
    msk_ice(j,1:idx(end)) = 1;
  end
  idx = find(abs(psiT_neg_mean(j,:))<1e-12,1,'first');
  msk(j,idx+1:end) = NaN;
end


psiT_pos_plot = -mean(psiT_pos_mean,3) * rho0*Cp/1e12 .* msk;
psiT_neg_plot = -mean(psiT_neg_mean,3) * rho0*Cp/1e12 .* msk;



%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.07 0.55 .86 .42];
axpos(2,:) = [0.07 0.05 .42 .42];
axpos(3,:) = [0.55 0.05 .42 .42];
cbpos = [0.96 0.55 0.01 .42];
axlabels = {'(a)','(b)','(c)'};
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
title('Total heat function');
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
title('``Warm'''' heat function');
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
title('``Cold'''' heat function');
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


%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


figure1 = gcf;

% Create arrow
annotation(figure1,'arrow',[0.577000000000001 0.522000000000001],...
  [0.86501079913607 0.8390928725702],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.524 0.508000000000002],...
  [0.874730021598272 0.901727861771067],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.636 0.622000000000002],...
  [0.943844492440605 0.966522678185756],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.86 0.838000000000004],...
  [0.938444924406048 0.961123110151201],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.78 0.73],...
  [0.939524838012959 0.958963282937366],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.44 0.388],...
  [0.818574514038877 0.816414686825057],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.792 0.736],...
  [0.874730021598273 0.876889848812098],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.902000000000002 0.855],...
  [0.893088552915769 0.909287257019438],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.285 0.272000000000001],...
  [0.844492440604752 0.874730021598279],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.417 0.402000000000002],...
  [0.862850971922246 0.893088552915775],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.817 0.789],...
  [0.404967602591794 0.390928725701944],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.822000000000001 0.796000000000002],...
  [0.375809935205184 0.356371490280782],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.729 0.704000000000001],...
  [0.319654427645789 0.31749460043197],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.661 0.65],...
  [0.339092872570195 0.362850971922247],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.726 0.715],...
  [0.366090712742981 0.389848812095033],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.764 0.753],...
  [0.375809935205184 0.399568034557237],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.461 0.437000000000001],...
  [0.436285097192225 0.466522678185748],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.373 0.358000000000001],...
  [0.43304535637149 0.46760259179266],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.422 0.401000000000002],...
  [0.434125269978402 0.464362850971927],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

