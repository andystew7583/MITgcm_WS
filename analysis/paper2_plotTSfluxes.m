%%%
%%% paper2_plotTSfluxes.m
%%%
%%% Plots heat/salt fluxes in quasi-latitude coordinates
%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_RTOPO2';
expname_notides = 'hires_seq_onetwelfth_notides_RTOPO2';
% expname = 'hires_seq_onethird_RTOPO2';
% expname_notides = 'hires_seq_onethird_notides_RTOPO2';

%%% Options (see calcTSfluxes)
deform_cavity = false;
outfname = ['','_TSfluxes'];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];
outfname_tides = [expname,outfname];
outfname_notides = [expname_notides,outfname];

%%% Load pre-computed fluxes
load(fullfile('products',outfname_tides));
thflux_tot_tides = thflux_stand+thflux_fluc+mean(thflux_eddy,2);
thflux_tot_plot = thflux_tot_tides;
sltflux_tot_tides = sltflux_stand+sltflux_fluc+mean(sltflux_eddy,2);
sltflux_tot_plot = sltflux_tot_tides;
load(fullfile('products',outfname_notides));
thflux_mean_plot = thflux_stand;
thflux_fluc_plot = thflux_fluc;
thflux_eddy_plot = mean(thflux_eddy,2);
thflux_tot_notides = thflux_stand+thflux_fluc+mean(thflux_eddy,2);
thflux_tide_plot = thflux_tot_tides - thflux_tot_notides;
sltflux_mean_plot = sltflux_stand;
sltflux_fluc_plot = sltflux_fluc;
sltflux_eddy_plot = mean(sltflux_eddy,2);
sltflux_tot_notides = sltflux_stand+sltflux_fluc+mean(sltflux_eddy,2);
sltflux_tide_plot = sltflux_tot_tides - sltflux_tot_notides;
thflux_eddy_adv_plot = mean(thflux_eddy_adv,2);
thflux_eddy_stir_plot = mean(thflux_eddy_stir,2);
sltflux_eddy_adv_plot = mean(sltflux_eddy_adv,2);
sltflux_eddy_stir_plot = mean(sltflux_eddy_stir,2);






%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.07 0.53 .41 .4];
axpos(2,:) = [0.56 0.53 .41 .4];
axpos(3,:) = [0.07 0.06 .41 .4];
axpos(4,:) = [0.56 0.06 .41 .4];
cbpos = [0.95 0.05 0.015 .93];
axlabels = {'(a)','(b)'};
rho0 = 1027;
Cp = 4000;
colororder = get(gca,'ColorOrder');
linewidth = 1.5;



%%% Set up the figure
figure(211)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[382    39   1000   746]);

%%% Heat flux
subplot('Position',axpos(1,:));
plot(eta,rho0*Cp*thflux_mean_plot/1e12,'LineWidth',linewidth);
hold on;
plot(eta,rho0*Cp*thflux_fluc_plot/1e12,'LineWidth',linewidth);
plot(eta,rho0*Cp*thflux_eddy_plot/1e12,'LineWidth',linewidth);
plot(eta,rho0*Cp*thflux_tide_plot/1e12,'LineWidth',linewidth);
plot(eta,rho0*Cp*thflux_tot_plot/1e12,'LineWidth',linewidth);
plot(eta,0*eta,'k--');
plot([4 4],[-15 5],'--','Color',[.3 .3 .3]);
plot([0 0],[-15 5],'--','Color',[.3 .3 .3]);
hold off;
ylabel('Heat flux (TW)');
axis([-8 10 -15 5]);
set(gca,'FontSize',fontsize);
leghandle = legend('Mean','Seasonal/interannual','Eddy','Tide - No Tide','Total','Location','West');
set(leghandle,'FontSize',fontsize);
text(-4.5,-14,'FRIS','FontSize',fontsize);
text(1.3,-14,'Shelf','FontSize',fontsize);
text(5.5,-14,'Weddell Sea','FontSize',fontsize);

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

%%% Salt flux
subplot('Position',axpos(2,:));
plot(eta,rho0*sltflux_mean_plot/1e9,'LineWidth',linewidth);
hold on;
plot(eta,rho0*sltflux_fluc_plot/1e9,'LineWidth',linewidth);
plot(eta,rho0*sltflux_eddy_plot/1e9,'LineWidth',linewidth);
plot(eta,rho0*sltflux_tide_plot/1e9,'LineWidth',linewidth);
plot(eta,rho0*sltflux_tot_plot/1e9,'LineWidth',linewidth);
plot(eta,0*eta,'k--');
plot([4 4],[-.4 1.4],'--','Color',[.3 .3 .3]);
plot([0 0],[-.4 1.4],'--','Color',[.3 .3 .3]);
hold off;
ylabel('Salt flux (Gg/s)');
% xlabel('MOC coordinate, \eta');
axis([-8 10 -.4 1.4]);
set(gca,'FontSize',fontsize);
text(-4.5,-0.3,'FRIS','FontSize',fontsize);
text(1.3,-0.3,'Shelf','FontSize',fontsize);
text(5.5,-0.3,'Weddell Sea','FontSize',fontsize);

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

%%% Eddy heat flux
subplot('Position',axpos(3,:));
plot(eta,rho0*Cp*thflux_eddy_plot/1e12,'Color',colororder(3,:),'LineWidth',linewidth);
hold on;
plot(eta,rho0*Cp*thflux_eddy_adv_plot/1e12,'Color',colororder(6,:),'LineWidth',linewidth);
plot(eta,rho0*Cp*thflux_eddy_stir_plot/1e12,'Color',colororder(7,:),'LineWidth',linewidth);
plot(eta,0*eta,'k--');
plot([4 4],[-8 2],'--','Color',[.3 .3 .3]);
plot([0 0],[-8 2],'--','Color',[.3 .3 .3]);
hold off;
ylabel('Heat flux (TW)');
xlabel('MOC coordinate, \eta');
axis([-8 10 -8 2]);
set(gca,'FontSize',fontsize);
leghandle = legend('Eddy','Eddy advection','Eddy stirring','Location','West');
set(leghandle,'FontSize',fontsize);
text(-4.5,-7.5,'FRIS','FontSize',fontsize);
text(1.3,-7.5,'Shelf','FontSize',fontsize);
text(5.5,-7.5,'Weddell Sea','FontSize',fontsize);

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

%%% Eddy salt flux
subplot('Position',axpos(4,:));
plot(eta,rho0*sltflux_eddy_plot/1e9,'Color',colororder(3,:),'LineWidth',linewidth);
hold on;
plot(eta,rho0*sltflux_eddy_adv_plot/1e9,'Color',colororder(6,:),'LineWidth',linewidth);
plot(eta,rho0*sltflux_eddy_stir_plot/1e9,'Color',colororder(7,:),'LineWidth',linewidth);
plot(eta,0*eta,'k--');
plot([4 4],[-.1 .3],'--','Color',[.3 .3 .3]);
plot([0 0],[-.1 .3],'--','Color',[.3 .3 .3]);
hold off;
ylabel('Salt flux (Gg/s)');
xlabel('MOC coordinate, \eta');
axis([-8 10 -.1 .3]);
set(gca,'FontSize',fontsize);
text(-4.5,-.08,'FRIS','FontSize',fontsize);
text(1.3,-.08,'Shelf','FontSize',fontsize);
text(5.5,-.08,'Weddell Sea','FontSize',fontsize);

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