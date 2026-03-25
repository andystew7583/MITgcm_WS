%%%
%%% plotVerticalFluxes.m
%%% 
%%% Plots vertical eddy heat/salt/buoyancy fluxes.
%%%

%%% Load experiment data
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
deform_cavity = false;
loadexp;

%%% MOC grid
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
eta = -9:.1:11;
Neta = length(eta);

%%% T/S flux storage file
infname = [expname,'_TSfluxes'];
if (deform_cavity)
  infname = [infname,'_deform'];
end
infname = [infname,'.mat'];

%%% Load pre-computed vertical flux data
load(fullfile('products',infname),'wvel_tavg','theta_tavg','salt_tavg','wvelth_tavg','wvelslt_tavg','tflux_tavg','sflux_tavg');
tflux_tavg(sum(hFacC,3)==0) = NaN;
sflux_tavg(sum(hFacC,3)==0) = NaN;

%%% Calculate vertical fluxes
theta_f = -1.85;
salt_ref = 34.45;
theta_w = 0.5*(theta_tavg(:,:,1:Nr)+theta_tavg(:,:,[2:Nr 1]));
wvelth_ref = theta_f.*wvel_tavg;
wvelth_mean = theta_w.*wvel_tavg;
wvelth_eddy = wvelth_tavg - wvelth_mean;
salt_w = 0.5*(salt_tavg(:,:,1:Nr)+salt_tavg(:,:,[2:Nr 1]));
wvelslt_ref = salt_ref.*wvel_tavg;
wvelslt_mean = salt_w.*wvel_tavg;
wvelslt_eddy = wvelslt_tavg - wvelslt_mean;
theta_tavg(hFacC==0) = NaN;

                
%%% Plotting options                
fontsize = 14;
framepos = [417   455   986   445];
labelspacing = 200;
axpos = zeros(2,4);
axpos(1,:) = [0.05 0.11 0.45 0.84];
axpos(2,:) = [0.64 0.11 0.34 0.84];
cb1_pos = [0.52 0.11 0.015 0.84];
axlabels = {'(a)','(b)'};

kidx = 23;
msk = real((SHELFICEtopo>=0) & (ETA<3.5) & (bathy<0) & (bathy<zz(kidx)));
plotmsk = real(msk);
plotmsk(plotmsk==0) = NaN;
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
bathy_plot = bathy;
bathy_plot(hFacC(:,:,kidx)==0) = NaN;

iceshelf_msk = real((SHELFICEtopo<0) & (SHELFICEtopo>bathy));
iceshelf_msk(iceshelf_msk==0) = NaN;

figure(205);
clf;
set(gcf,'Position',framepos);

wvelth_eddy_plot = wvelth_eddy(:,:,23);
wvelth_eddy_plot(sum(hFacC,3)==0) = NaN;
subplot('Position',axpos(1,:));
pcolor(XC,YC,wvelth_eddy_plot*rho0*Cp);
shading interp;
hold on;
[C,h] = contour(XC,YC,-bathy_plot,bathycntrs,'EdgeColor','k');
clabel(C,h);
hold off;
cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
title(cbhandle,'W/m^2')
colormap(gca,cmocean('balance',40));
caxis([-100 100]);
set(gca,'XLim',[-63 -25]);
set(gca,'YLim',[-78.5 -71.3]);
set(gca,'FontSize',fontsize);
xlabel('Longitude');
ylabel('Latitude');
set(gca,'Color',[.8 .8 .8]);
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
% hold on;
handle = pcolor(XC,YC,-real(msk));
shading interp;
colormap(ax2,hot);
caxis(ax2,[-1 0]);
set(handle,'FaceAlpha',0.1);
set(ax2,'Color','None');
set(ax2,'XLim',[-63 -25]);
set(ax2,'YLim',[-78.5 -71.3]);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
ax3 = axes('Position',get(ax1,'Position'));
handle = pcolor(XC,YC,real(iceshelf_msk));
shading interp;
colormap(ax3,redblue);
set(handle,'FaceAlpha',1);
set(ax3,'Color','None');
set(ax3,'XLim',[-63 -25]);
set(ax3,'YLim',[-78.5 -71.3]);
set(ax3,'XTick',[]);
set(ax3,'YTick',[]);

sum(sum(rho0.*Cp.*wvelth_eddy(:,:,23).*msk.*RAC)) / sum(sum(msk.*RAC))
nansum(nansum(-tflux_tavg(:,:).*msk.*RAC)) / sum(sum(msk.*RAC))

subplot('Position',axpos(2,:));
plot(squeeze(sum(sum((wvelth_tavg-wvelth_ref).*RAC.*msk)))*rho0*Cp/sum(sum(msk.*RAC)),-squeeze(zz));
hold on;
plot(squeeze(sum(sum((wvelth_mean-wvelth_ref).*RAC.*msk)))*rho0*Cp/sum(sum(msk.*RAC)),-squeeze(zz));
plot(squeeze(sum(sum(wvelth_eddy.*RAC.*msk)))*rho0*Cp/sum(sum(msk.*RAC)),-squeeze(zz));
plot([-2 7],-[zz(kidx) zz(kidx)],'k--','LineWidth',1);
hold off
set(gca,'XLim',[-2 7]);
xlabel('Vertical heat flux (W/m^2)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy','Location','SouthEast');
set(gca,'YLim',[0 1000]);
% title(['Ice shelf face to \eta=3.5, total surface heat flux = ',num2str(-nansum(nansum(tflux_tavg.*RAC.*msk/1e12))),' TW']);
set(gca,'FontSize',fontsize);

annotation('textbox',[0.02 0.02 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.57 0.02 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');