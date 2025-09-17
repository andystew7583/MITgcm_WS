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
% salt0 = 34.73;
% 
% %%% Set true to deform coordinates in the cavity
% deform_cavity = false;
% 
% %%% Set true to use barotropic streamfunction as the coordinate system
% use_PsiBT = false;
% 
% %%% Set true to use depth-averaged temperature as the coordinate system
% use_meanT = true;
% 
% %%% Set true to decompose eddy fluxes
% calc_eddy_decomp = false;
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
% %%% Some quantities computed separately for higher-resolution experiments
% if (strcmp(expname,'hires_seq_onetwentyfourth_notides_RTOPO2'))
% 
%   %%% Load eddy-induced components of heatfunction/transport
%   outfname = [expname,'_HeatFunctionEddyDecomp'];
%   if (use_PsiBT)
%     outfname = [outfname,'_PsiBT'];
%   else
%     if (use_meanT)
%       outfname = [outfname,'_meanT'];
%     else 
%       if (deform_cavity)
%         outfname = [outfname,'_deform'];
%       end
%     end
%   end
%   outfname = [outfname,'.mat'];
%   load(fullfile('products',outfname));
% 
%   psiT_eddy_adv = psiT_eddy_adv - psi_eddy*theta0;
%   psiT_eddy_stir = psiT_eddy - psiT_eddy_adv;
% 
% end
% %%% Latent heat of freezing
% Lf = 3.34e5;
% 
% 
% 
% 
% psiT_mod = psiT_tot-psi_tot*theta0;
% psiT_mean_mod = psiT_mean-psi_tot*theta0;
% 
% thflux_tot_plot = -mean(psiT_mod(:,1,:),3) * rho0*Cp/1e12;
% thflux_mean_plot = -mean(psiT_mean_mod(:,1,:),3)* rho0*Cp/1e12;
% thflux_eddy_plot = -mean(psiT_eddy(:,1,:),3)* rho0*Cp/1e12;
% thflux_eddy_adv_plot = -mean(psiT_eddy_adv(:,1,:),3)* rho0*Cp/1e12;
% thflux_eddy_stir_plot = -mean(psiT_eddy_stir(:,1,:),3)* rho0*Cp/1e12;
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%
% %%% MAKE PLOTS %%%
% %%%%%%%%%%%%%%%%%%
% 
% %%% Plotting options
% fontsize = 14;
% bathycntrs = [0 250 500 1000 2000 3000 4000];
% axpos = zeros(4,4);
% axpos(1,:) = [0.02 0.45 0.9 .5];
% axpos(2,:) = [0.1 0.06 .81 .39];
% cbpos = [0.84 0.56 0.01 .1];
% axlabels = {'\textbf{A}','\textbf{B}','\textbf{C}','\textbf{D}'};
% rho0 = rhoConst;
% Cp = 4000;
% colororder = get(gca,'ColorOrder');
% linewidth = 1.5;
% xlim = [-2.6 0];
% ylim = [0 2000];
% 
% %%% Plotting range for salinity figure
% latMin_b = min(min(YC));
% latMax_b = YC(1,end-spongethickness);
% lonMin_b = min(min(XC));
% lonMax_b = XC(end-spongethickness,1);







%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

%%% Set up the figure
figure(207)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417    34  600  926]);












%%%%%%%%%%%%%%%%%%%%%%
%%% ETA COORDINATE %%%
%%%%%%%%%%%%%%%%%%%%%%



%%% Plotting options
clim = [-2.5 -1];
% cmap = cmocean('balance',50);
% cmap = flip(haxby(20),1);
% cmap = cmap(8:23,:);
cmap = cmocean('thermal',15);
% cmap = cmocean('amp',16);
% cmap = pmkmp(16,'swtth');
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];


%%% Set up map plot
subplot('Position',axpos(1,:)+[.1 0 0 0]);
axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[latMin_b latMax_b], ...
  'MapLonLimit',[lonMin_b lonMax_b], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', [-70:10:-30],...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');  
axis off;
setm(gca,'MLabelParallel',-20)

%%% Plot MOC coordinate
ETA_plot = ETA;
ETA_plot(sum(hFacC,3)==0) = NaN;
pcolorm(YC,XC,ETA_plot);
shading interp
caxis(clim);
colormap(gca,cmap);


%%% Add colorbar and title
h = colorbar;
set(gca,'FontSize',fontsize);
set(h,'Position',cbpos)
title(h,'$^\circ$C','FontSize',fontsize,'interpreter','latex');
tightmap;

%%% Add bathymetry contours
hold on;
[cs,C] = contourm(YC,XC,SHELFICEtopo-bathy,bathycntrs,'EdgeColor',[.25 .25 .25]); 
% [cs,C] = contourm(YC,XC,ETA,clim(1):.1:clim(2),'EdgeColor',[.25 .25 .25]); 
hh = clabelm(cs,C);
set(hh,'fontsize',8,'Color',[.05 .05 .05],'BackgroundColor','none','Edgecolor','none')       
hold off;

%%% Add axis labels
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
set(gca,'Position',axpos(1,:));
hold off

% Create title
annotation(gcf,'textbox',...
  [0.24 0.952023758099355 0.7005 0.0242980561555075],...
  'String',{'Mean depth-averaged potential temperature, $\widetilde{\theta}$'},'EdgeColor','None','FontSize',fontsize+2,'interpreter','latex');







%%% Heat flux
axes('Position',axpos(2,:));
plot(eta,thflux_mean_plot,'Color',colororder(6,:),'LineWidth',linewidth);
hold on;
plot(eta,thflux_eddy_plot,'Color',colororder(5,:),'LineWidth',linewidth);
plot(eta,thflux_eddy_adv_plot,'Color',colororder(5,:),'LineWidth',linewidth,'LineStyle',':');
plot(eta,thflux_eddy_stir_plot,'Color',colororder(5,:),'LineWidth',linewidth,'LineStyle','--');
plot(eta,thflux_tot_plot,'Color',colororder(4,:),'LineWidth',linewidth);
% plot(eta,-surfQint/1e12,'Color',colororder(4,:),'LineWidth',linewidth);
% plot(eta,0*eta,'k--');
plot([4 4],[-5 15],'--','Color',[.3 .3 .3],'LineWidth',2);
plot([0 0],[-5 15],'--','Color',[.3 .3 .3],'LineWidth',2);
hold off;
set(gca,'FontSize',fontsize);
ylabel('Shoreward heat flux (TW)');
xlabel('$\widetilde{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2);
axis([-2.5 -1.5 -5 8]);
leghandle = legend('Mean','Eddy','Eddy advection','Eddy stirring','Total','Location','NorthWest');
set(leghandle,'FontSize',fontsize);
text(-4.5,-3,'FRIS','FontSize',fontsize);
text(1.3,-3,'Shelf','FontSize',fontsize);
grid on;
box on;




%%% Add panel labels
annotation('textbox',[axpos(1,1)+0.01 axpos(1,2)+0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.06 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
