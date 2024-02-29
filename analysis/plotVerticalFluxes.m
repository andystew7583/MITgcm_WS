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

kidx = 23;
msk = real((SHELFICEtopo>=0) & (ETA<3.5) & (bathy<0) & (bathy<zz(kidx)));
plotmsk = real(msk);
plotmsk(plotmsk==0) = NaN;
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
bathy_plot = bathy;
bathy_plot(hFacC(:,:,kidx)==0) = NaN;

figure(1);
clf;
pcolor(XC,YC,wvelth_eddy(:,:,23)*rho0*Cp);
shading interp;
hold on;
[C,h] = contour(XC,YC,-bathy_plot,bathycntrs,'EdgeColor','k');
clabel(C,h);
hold off;
colorbar;
colormap(gca,cmocean('balance',40));
caxis([-100 100]);
set(gca,'XLim',[-63 -25]);
set(gca,'YLim',[-78.5 -71.3]);
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
% hold off;

sum(sum(rho0.*Cp.*wvelth_eddy(:,:,23).*msk.*RAC)) / sum(sum(msk.*RAC))
nansum(nansum(-tflux_tavg(:,:).*msk.*RAC)) / sum(sum(msk.*RAC))



figure(2);
clf;
pcolor(XC,YC,wvelslt_eddy(:,:,23));
shading interp;
colorbar;
colormap redblue;
caxis([-1 1]*1e-5);

figure(3);
clf;
pcolor(XC,YC,wvelth_mean(:,:,23)-wvelth_ref(:,:,23));
shading interp;
colorbar;
colormap redblue;
caxis([-1 1]*1e-5);

figure(4);
clf;
pcolor(XC,YC,wvelslt_mean(:,:,23));
shading interp;
colorbar;
colormap redblue;
caxis([-1 1]*5e-4);

figure(5);
clf;
pcolor(XC,YC,wvelth_tavg(:,:,23)-wvelth_ref(:,:,23));
shading interp;
colorbar;
colormap redblue;
caxis([-1 1]*1e-5);

figure(6);
clf;
pcolor(XC,YC,wvelslt_tavg(:,:,23));
shading interp;
colorbar;
colormap redblue;
caxis([-1 1]*5e-4);

figure(7);
clf;
pcolor(XC,YC,wvel_tavg(:,:,23));
shading interp;
colorbar;
colormap redblue;
caxis([-1 1]*5e-4);

figure(8);
clf;
pcolor(XC,YC,theta_tavg(:,:,15));
shading interp;
colorbar;
colormap jet;
caxis([-2 -1]);

msk = (SHELFICEtopo>=0) & (ETA<3.5);

figure(9);
clf;
plot(squeeze(sum(sum((wvelth_tavg-wvelth_ref).*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
hold on;
plot(squeeze(sum(sum((wvelth_mean-wvelth_ref).*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
plot(squeeze(sum(sum(wvelth_eddy.*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
hold off
xlabel('Vertical heat flux (TW)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy');
set(gca,'YLim',[0 1000]);
title(['Ice shelf face to \eta=3.5, total surface heat flux = ',num2str(-nansum(nansum(tflux_tavg.*RAC.*msk/1e12))),' TW']);

figure(10);
clf;
plot(squeeze(sum(sum((wvelslt_tavg-wvelslt_ref).*RAC.*msk)))*rho0/1e9,-squeeze(zz));
hold on;
plot(squeeze(sum(sum((wvelslt_mean-wvelslt_ref).*RAC.*msk)))*rho0/1e9,-squeeze(zz));
plot(squeeze(sum(sum(wvelslt_eddy.*RAC.*msk)))*rho0/1e9,-squeeze(zz));
hold off
xlabel('Vertical salt flux (Gg/s)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy');
set(gca,'YLim',[0 1000]);
title(['Ice shelf face to \eta=3.5, total surface salt flux = ',num2str(-nansum(nansum(sflux_tavg.*RAC.*msk))/1e9),' Gg/s']);

msk = (ETA>0.5) & (ETA<3.5);

figure(11);
clf;
plot(squeeze(sum(sum((wvelth_tavg-wvelth_ref).*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
hold on;
plot(squeeze(sum(sum((wvelth_mean-wvelth_ref).*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
plot(squeeze(sum(sum(wvelth_eddy.*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
hold off
xlabel('Vertical heat flux (TW)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy');
set(gca,'YLim',[0 1000]);
title(['0.5<\eta<3.5, total surface heat flux = ',num2str(-nansum(nansum(tflux_tavg.*RAC.*msk))/1e12),' TW']);

figure(12);
clf;
plot(squeeze(sum(sum((wvelslt_tavg-wvelslt_ref).*RAC.*msk)))*rho0/1e9,-squeeze(zz));
hold on;
plot(squeeze(sum(sum((wvelslt_mean-wvelslt_ref).*RAC.*msk)))*rho0/1e9,-squeeze(zz));
plot(squeeze(sum(sum(wvelslt_eddy.*RAC.*msk)))*rho0/1e9,-squeeze(zz));
hold off
xlabel('Vertical salt flux (Gg/s)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy');
set(gca,'YLim',[0 1000]);
title(['0.5<\eta=3.5, total surface salt flux = ',num2str(-nansum(nansum(sflux_tavg.*RAC.*msk))/1e9),' Gg/s']);



msk = (SHELFICEtopo>=0) & (ETA<0.5);


figure(13);
clf;
plot(squeeze(sum(sum((wvelth_tavg-wvelth_ref).*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
hold on;
plot(squeeze(sum(sum((wvelth_mean-wvelth_ref).*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
plot(squeeze(sum(sum(wvelth_eddy.*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
hold off
xlabel('Vertical heat flux (TW)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy');
set(gca,'YLim',[0 1000]);
title(['Ice shelf face to \eta=0.5, total surface heat flux = ',num2str(-nansum(nansum(tflux_tavg.*RAC.*msk/1e12))),' TW']);

figure(14);
clf;
plot(squeeze(sum(sum((wvelslt_tavg-wvelslt_ref).*RAC.*msk)))*rho0/1e9,-squeeze(zz));
hold on;
plot(squeeze(sum(sum((wvelslt_mean-wvelslt_ref).*RAC.*msk)))*rho0/1e9,-squeeze(zz));
plot(squeeze(sum(sum(wvelslt_eddy.*RAC.*msk)))*rho0/1e9,-squeeze(zz));
hold off
xlabel('Vertical salt flux (Gg/s)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy');
set(gca,'YLim',[0 1000]);
title(['Ice shelf face to \eta=0.5, total surface salt flux = ',num2str(-nansum(nansum(sflux_tavg.*RAC.*msk))/1e9),' Gg/s']);

msk = (ETA>0.5) & (SHELFICEtopo>=0) & (bathy >= -1000) & (YC<YC(1,end-spongethickness+1)) & (XC<XC(end-spongethickness+1,1));

figure(15);
clf;
plot(squeeze(sum(sum((wvelth_tavg-wvelth_ref).*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
hold on;
plot(squeeze(sum(sum((wvelth_mean-wvelth_ref).*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
plot(squeeze(sum(sum(wvelth_eddy.*RAC.*msk)))*rho0*Cp/1e12,-squeeze(zz));
hold off
xlabel('Vertical heat flux (TW)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy');
set(gca,'YLim',[0 1000]);
title(['Cont. shelf outside cavity, total surface heat flux = ',num2str(-nansum(nansum(tflux_tavg.*RAC.*msk/1e12))),' TW']);

figure(16);
clf;
plot(squeeze(sum(sum((wvelslt_tavg-wvelslt_ref).*RAC.*msk)))*rho0/1e9,-squeeze(zz));
hold on;
plot(squeeze(sum(sum((wvelslt_mean-wvelslt_ref).*RAC.*msk)))*rho0/1e9,-squeeze(zz));
plot(squeeze(sum(sum(wvelslt_eddy.*RAC.*msk)))*rho0/1e9,-squeeze(zz));
hold off
xlabel('Vertical salt flux (Gg/s)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy');
set(gca,'YLim',[0 1000]);
title(['Cont. shelf outside cavity, total surface salt flux = ',num2str(-nansum(nansum(sflux_tavg.*RAC.*msk))/1e9),' Gg/s']);



%%% Make room in memory
% clear('wvel','wvelth_tavg','wvelslt_tavg','wvelth_mean','wvelth_eddy','wvelslt_mean','wvelslt_eddy','theta_f','salt_w','theta_w');




