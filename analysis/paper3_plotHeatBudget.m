%%%
%%% paper3_plotHeatBudget.m
%%%
%%% Plots heat/salt budget in quasi-latitude coordinates
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


%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.17 0.65 .66 .32];
axpos(2,:) = [0.07 0.05 .42 .42];
axpos(3,:) = [0.55 0.05 .42 .42];
cbpos = [0.96 0.55 0.01 .42];
axlabels = {'(a)','(b)','(c)'};
rho0 = 1027;
Cp = 4000;
arrowcolor = [0.301960784313725 0.55098039215686 0.93333333333333];
colororder = get(gca,'ColorOrder');
linewidth = 1.5;
ylim = [0 2000];
xlim = [-8.6 4];
icecolor = [186 242 239]/255;
[ZZ,EE] = meshgrid(squeeze(-RF),eta);






%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(204)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34  1000  926]);


eidx_icefront = find(abs(eta-eta_icefront)<1e-8);
zminidx_icefront = find(isnan(msk(eidx_icefront,:)),1);
eidx_shelfbreak = find(abs(eta-eta_shelfbreak)<1e-8);
zminidx_shelfbreak = find(isnan(msk(eidx_shelfbreak,:)),1);



%%% TODO add time series of vertical eddy heat flux and eddy energy
%%% production, pos/neg heat flux into cavity

%%% TODO shade region used to compute averages of vertical buoyancy fluxes
%%% Schematic 
axes('Position',axpos(1,:));
pcolor(EE,ZZ,psiT_tot_plot);
shading flat;
hold on;
clabel(C,h);
plot([eta_icefront,eta_shelfbreak],[-RC(zidx_icefront),-RC(zidx_icefront)],'k--','LineWidth',2);
plot([eta_icefront,eta_icefront],[-RC(zidx_icefront),-RC(zminidx_icefront)],'k--','LineWidth',2);
plot([eta_shelfbreak,eta_shelfbreak],[-RC(zidx_icefront),-RC(zminidx_shelfbreak)],'k--','LineWidth',2);
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


%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


figure1 = gcf;

% Create arrow
annotation(figure1,'arrow',[0.675 0.675],...
  [0.900647948164147 0.970842332613394],...
  'Color',arrowcolor,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create arrow
annotation(figure1,'arrow',[0.826000000000001 0.756000000000002],...
  [0.886609071274299 0.887688984881212],...
  'Color',arrowcolor,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

% Create arrow
annotation(figure1,'arrow',[0.594 0.52],...
  [0.83585313174946 0.835853131749461],...
  'Color',arrowcolor,...
  'LineWidth',30,...
  'HeadWidth',60,...
  'HeadStyle','plain',...
  'HeadLength',35);

