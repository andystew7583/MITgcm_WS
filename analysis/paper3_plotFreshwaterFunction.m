% %%%
% %%% paper3_plotTSfluxes.m
% %%%
% %%% Plots heat/salt fluxes in quasi-latitude coordinates
% %%%
% 
% %%% Load experiment
% expdir = '../experiments';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% % loadexp;
% 
% %%% Reference surface freezing temperature
% theta0 = -1.9;
% 
% %%% Reference salinity
% salt0 = 34.6;
% % salt0 = 34.72;
% 
% %%% Cross-shelf locations of ice front and shelf break for heat budget calculation
% eta_icefront = -1.1;
% eta_shelfbreak = 3.5;  
% 
% %%% Depth of ice front for heat budget calculation
% zidx_icefront = 25;    
% 
% %%% Set true to deform coordinates in the cavity
% deform_cavity = false;
% 
% %%% Set true to use barotropic streamfunction as the coordinate system
% use_PsiBT = false;
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
%     end
%   end
% end
% outfname = [outfname,'.mat'];
% load(fullfile('products',outfname));
% 
% %%% Mask for ice/land
% psiF_tot_mean = mean(psiS_tot-psi_tot*salt0,3)/ -salt0 * rhoFresh*t1year / 1e12;
% msk = ones(size(psiF_tot_mean));
% msk_ice = NaN*msk;
% for j=1:Neta  
%   idx = find(psiF_tot_mean(j,:)==psiF_tot_mean(j,1));
%   idx(end) = [];
%   msk(j,idx) = NaN;
%   if (~isempty(idx))
%     msk_ice(j,1:idx(end)) = 1;
%   end
%   idx = find(abs(psiF_tot_mean(j,:))<1e-12,1,'first');
%   msk(j,idx+1:end) = NaN;
% end
% 

psiF_tot_plot = -mean(psiS_tot-psi_tot*salt0,3) / -salt0 * rhoFresh*t1year / 1e12 .* msk;
psiF_mean_plot = -mean(psiS_mean-psi_tot*salt0,3) / -salt0 * rhoFresh*t1year / 1e12 .* msk;
psiF_eddy_plot = psiF_tot_plot - psiF_mean_plot;



%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(2,4);
axpos(1,:) = [0.07 0.55 .86 .42];
axpos(2,:) = [0.07 0.05 .86 .42];
cbpos1 = [0.95 0.55 0.01 .42];
cbpos2 = [0.95 0.05 0.01 .42];
axlabels = {'(a)','(b)'};
rho0 = 1027;
Cp = 4000;
arrowcolor = [0.301960784313725 0.35098039215686 0.93333333333333];
colororder = get(gca,'ColorOrder');
linewidth = 1.5;
ylim = [0 2000];
psiFmax = 800;
psiFstep = 50;
psiFsteps = [-psiFmax:psiFstep:-psiFstep psiFstep:psiFstep:psiFmax];
psimax = 1.2;
psistep = 0.1;
psisteps = [-psimax:psistep:-psistep psistep:psistep:psimax];
xlim = [-8.6 4];
icecolor = [186 242 239]/255;
[ZZ,EE] = meshgrid(squeeze(-RF),eta);






%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(206)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34  1000  926]);





%%% Total heat function
axes('Position',axpos(1,:));
pcolor(EE,ZZ,psiF_tot_plot);
shading flat;
hold on;
[C,h] = contour(EE,ZZ,psiF_tot_plot,[-psiFmax:psiFstep:-psiFstep psiFstep:psiFstep:psiFmax],'EdgeColor','k');
clabel(C,h);
hold off;
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(gca,cmocean('balance',2*length(psiFsteps)));
cbhandle = colorbar;
set(cbhandle,'Position',cbpos1);
title(cbhandle,'Gt/yr','FontSize',fontsize);
xlabel('Cross-shelf coordinate, \eta');
ylabel('Depth (m)');
title('Total freshwater function');
set(gca,'FontSize',fontsize);
caxis([-psiFmax psiFmax]);

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


psi_tot_plot = -mean(psi_eddy,3)/1e6.*msk;

%%% Total heat function
axes('Position',axpos(2,:));
pcolor(EE,ZZ,psi_tot_plot);
shading flat;
hold on;
[C,h] = contour(EE,ZZ,psi_tot_plot,[-psimax:psistep:-psistep psistep:psistep:psimax],'EdgeColor','k');
clabel(C,h);
hold off;
set(gca,'YDir','reverse');
set(gca,'XLim',xlim);
set(gca,'YLim',ylim);
set(gca,'Color',[.7 .7 .7])
% colormap redblue(32);
% colormap(cmocean('balance',round(2*psimax/psistep)));
colormap(gca,cmocean('balance',2*length(psisteps)));
cbhandle = colorbar;
set(cbhandle,'Position',cbpos2);
title(cbhandle,'Sv','FontSize',fontsize);
xlabel('Cross-shelf coordinate, \eta');
ylabel('Depth (m)');
title('Eulerian-mean streamfunction');
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



figure1 = gcf;

% arrowcolor = [0.650980392156863 0.650980392156863 0.650980392156863];
arrowcolor = [0 0 0];
arrowcolor = [14*16+4 7*16+3 14*16+4]/255;

% Create arrow
annotation(figure1,'arrow',[0.896 0.839],...
  [0.929805615550756 0.942764578833699],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.889 0.832],...
  [0.882289416846654 0.895248380129596],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.347 0.408],...
  [0.817494600431968 0.838012958963286],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.529 0.59],...
  [0.879049676025919 0.899568034557236],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.694 0.663000000000001],...
  [0.938444924406048 0.964362850971926],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.783 0.740000000000001],...
  [0.941684665226782 0.951403887688988],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.634 0.621000000000001],...
  [0.926565874730022 0.962203023758101],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.757 0.707],...
  [0.896328293736504 0.915766738660911],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.903 0.84],...
  [0.442764578833693 0.440604751619876],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.695 0.748000000000001],...
  [0.348812095032397 0.352051835853136],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.835000000000001 0.900000000000002],...
  [0.360691144708424 0.359611231101516],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.759 0.698000000000001],...
  [0.43952483801296 0.439524838012962],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.644 0.644000000000001],...
  [0.422246220302376 0.368250539956807],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.299 0.339],...
  [0.31317494600432 0.362850971922246],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.364 0.332],...
  [0.318574514038878 0.276457883369331],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.458 0.503000000000001],...
  [0.242980561555076 0.265658747300221],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.481 0.431000000000001],...
  [0.335853131749461 0.311015118790502],...
  'Color',arrowcolor,...
  'LineWidth',2,...
  'HeadStyle','plain');

