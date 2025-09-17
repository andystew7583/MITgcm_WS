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
salt0 = 34.73;
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

%%% Cross-shelf locations of ice front and shelf break for heat budget calculation
if (gl_coord)
  eta_icefront = 0;
else
  eta_icefront = -1.1;
end
eta_shelfbreak = 3.5;  

%%% Parameters
rho0 = 1000;
Cp = 4000;
rhoFresh = 1000;

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
psiF_tot_mean = mean(psiS_tot-psi_tot*salt0,3)/ -salt0 * rhoFresh*t1year / 1e12;
msk = ones(size(psiF_tot_mean));
msk_ice = NaN*msk;
for j=1:Neta  
  idx = find(psiF_tot_mean(j,:)==psiF_tot_mean(j,1));
  idx(end) = [];
  icefrontidx = find(eta==eta_icefront);
  if (~gl_coord || (j~=icefrontidx))  
    msk(j,idx) = NaN;
    if (~isempty(idx))
      msk_ice(j,1:idx(end)) = 1;
    end
  end  
  idx = find(abs(psiF_tot_mean(j,:))<1e-12,1,'first');
  msk(j,idx+1:end) = NaN;
end


%%% Compute mean latitude of eta isopleths
equiv_lats = NaN*eta;
eta_mid = 0.5*(eta(1:end-1)+eta(2:end));
msk_ocean = sum(hFacC,3)>0;
for j=2:Neta-1
  msk_ETA = (ETA>eta_mid(j-1)) & (ETA<=eta_mid(j));
  equiv_lats(j) = sum(YC.*msk_ETA.*msk_ocean) / sum(msk_ETA.*msk_ocean);
end
equiv_lats(equiv_lats==0) = NaN;


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
axpos(1,:) = [0.07 0.55 .86 .4];
axpos(2,:) = [0.07 0.05 .86 .4];
cbpos1 = [0.95 0.55 0.01 .42];
cbpos2 = [0.95 0.05 0.01 .42];
axlabels = {'\textbf{A}','\textbf{B}'};
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
xticks = [-8:2:4];
% eticks = [-83:1:-74];
% eticklabels = cell(1,length(eticks));
% for n=1:length(eticklabels)
%   eticklabels{n} = num2str(eticks(n));
% end
% eticklocs = interp1(equiv_lats(~isnan(equiv_lats)),eta(~isnan(equiv_lats)),eticks,'nearest');
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
% title('Total freshwater function');
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

ax3 = axes('Position',get(ax1,'Position'));
set(ax3,'Color','None');
set(ax3,'XAxisLocation','Top');
set(ax3,'YAxisLocation','Right');
set(ax3,'YLim',get(ax1,'YLim'));
set(ax3,'YTick',[]);
% set(ax3,'XTick',eticklocs);
% set(ax3,'XTickLabel',eticklabels);
box off;
set(ax3,'XLim',xlim-77);
set(ax3,'FontSize',fontsize);
set(get(ax3,'XLabel'),'String','Reference latitude');


psi_tot_plot = -mean(psi_tot,3)/1e6.*msk;

%%% Total overturning streamfunction function
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
% title('Eulerian-mean streamfunction');
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


ax3 = axes('Position',get(ax1,'Position'));
set(ax3,'Color','None');
set(ax3,'XAxisLocation','Top');
set(ax3,'YAxisLocation','Right');
set(ax3,'YLim',get(ax1,'YLim'));
set(ax3,'YTick',[]);
% set(ax3,'XTick',eticklocs);
% set(ax3,'XTickLabel',eticklabels);
box off;
set(ax3,'XLim',xlim-77);
set(ax3,'FontSize',fontsize);
set(get(ax3,'XLabel'),'String','Reference latitude');


%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');



figure1 = gcf;

% arrowcolor = [0.650980392156863 0.650980392156863 0.650980392156863];
arrowcolor = [0 0 0];
arrowcolor = [14*16+4 7*16+3 14*16+4]/255;

% Create arrow
annotation(figure1,'arrow',[0.661 0.661000000000001],...
  [0.414686825053996 0.360691144708427],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.741 0.688000000000001],...
  [0.909287257019438 0.938444924406052],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.664 0.662000000000001],...
  [0.892008639308855 0.93844492440605],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.536 0.597],...
  [0.857451403887691 0.877969762419008],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.37 0.408],...
  [0.874730021598272 0.838012958963286],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.798 0.755000000000001],...
  [0.931965442764579 0.941684665226785],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.907 0.85],...
  [0.916846652267819 0.929805615550762],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.455 0.500000000000001],...
  [0.22354211663067 0.246220302375815],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.481 0.431000000000001],...
  [0.322894168466524 0.298056155507565],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.294 0.334],...
  [0.293736501079914 0.34341252699784],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.362 0.33],...
  [0.302375809935206 0.260259179265659],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');


% Create arrow
annotation(figure1,'arrow',[0.889 0.832],...
  [0.882289416846654 0.895248380129596],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.903 0.84],...
  [0.442764578833693 0.440604751619876],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.695 0.748000000000001],...
  [0.348812095032397 0.352051835853136],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.835000000000001 0.900000000000002],...
  [0.360691144708424 0.359611231101516],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');

% Create arrow
annotation(figure1,'arrow',[0.759 0.698000000000001],...
  [0.43952483801296 0.439524838012962],...
  'Color',[0.894117647058824 0.450980392156863 0.894117647058824],...
  'LineWidth',2,...
  'HeadStyle','plain');


% Create textbox
annotation(figure1,'textbox',...
  [0.076 0.417984881209504 0.1535 0.0242980561555076],...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'\textbf{Eulerian-mean streamfunction}'});

% Create textbox
annotation(figure1,'textbox',...
  [0.0760000000000001 0.916904967602593 0.13 0.0242980561555075],...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off',...
  'EdgeColor','none',...
  'String',{'\textbf{Total freshwater function}'});
