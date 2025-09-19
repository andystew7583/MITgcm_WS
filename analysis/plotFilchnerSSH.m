% %%%
% %%% plotFilchnerSSH.m
% %%%
% %%% Plots variations of Filchner overflow transport and SSH variance.
% %%%
% 
% 
% %%% Options
% expdir = '../experiments';
% expname = 'hires_nest_onethirtieth_notides_RTOPO2';
% tmin = 1.01;
% tmax = 7.01;
% 
% %%% Load experiment
% loadexp;
% 
% %%% Select density variable in which to compute isopycnal fluxes
% densvar = 'PD0';
% % densvar = 'ND1';
% % densvar = 'ND2';
% % densvar = 'PT';
% 
% %%% Month index to start reading SSH data (first year is skipped in DSW
% %%% flux)
% startidx_ssh = 13;
% 
% outfname = [expname,'_ShelfHeatBudget.mat'];
% load(fullfile('products',outfname),'usq_eddy_int','vsq_eddy_int');
% EKE = 0.5*(usq_eddy_int+vsq_eddy_int);
% Hocean = sum(DRF.*hFacC,3);
% EKE = EKE./Hocean;
% 
% %%% Load pre-computed Filchner Trough transport
% outfname = [expname,'_FTtrans_',densvar];
% outfname = [outfname,'.mat'];
% load(fullfile('products',outfname));
% times_tmp = [366*86400 times];
% times = 0.5*(times_tmp(1:end-1)+times_tmp(2:end));
% 
% %%% Index of density level corresponding to max transport
% didx = find(dens_levs>=27.81,1);
% 
% %%% Load pre-computed SSH variance
% ncfname = fullfile('products',[expname,'_SSHvar.nc']);
% ssh_var_runmean = ncread(ncfname,'ssh_var_runmean');
% ssh_runmean = ncread(ncfname,'ssh_runmean');
% 
% %%% Load monthly-mean SSH variance
% ncfname = fullfile('products',[expname,'_SSHmonVar.nc']);
% ssh_mon_mean = ncread(ncfname,'ssh_mon_mean');
% ssh_mon_var = ncread(ncfname,'ssh_mon_var');
% ssh_mon_var(ssh_mon_var<0) = 0; %%% Sometimes <0 to machine precision
% ssh_mon_std = sqrt(ssh_mon_var);
% ssh_mon_test = ssh_mon_var.^(2/3);
% 
% %%% Compute correlation maps
% corr_mean = zeros(Nx,Ny);
% corr_var = zeros(Nx,Ny);
% corr_std = zeros(Nx,Ny);
% corr_EKE = zeros(Nx,Ny);
% corr_EKEprod = zeros(Nx,Ny);
% corr_EKE_surf = zeros(Nx,Ny);
% corr_EKE_50m = zeros(Nx,Ny);
% corr_test = zeros(Nx,Ny);
% for i=1:Nx
%   for j=1:Ny
%     corr_mean(i,j) = corr(squeeze(ssh_mon_mean(i,j,startidx_ssh:end)),T_ft(didx,:)');
%     corr_var(i,j) = corr(squeeze(ssh_mon_var(i,j,startidx_ssh:end)),T_ft(didx,:)');
%     corr_std(i,j) = corr(squeeze(ssh_mon_std(i,j,startidx_ssh:end)),T_ft(didx,:)');
%     corr_test(i,j) = corr(squeeze(ssh_mon_test(i,j,startidx_ssh:end)),T_ft(didx,:)');
%     corr_EKE(i,j) = corr(squeeze(EKE(i,j,startidx_ssh:end)),T_ft(didx,:)');
%     corr_EKEprod(i,j) = corr(squeeze(EKE(i,j,startidx_ssh:end)),T_ft(didx,:)');
%     corr_EKE_surf(i,j) = corr(squeeze(EKE_surf(i,j,startidx_ssh:end)),T_ft(didx,:)');    
%     corr_EKE_50m(i,j) = corr(squeeze(EKE_50m(i,j,startidx_ssh:end)),T_ft(didx,:)');
%   end
% end

%%% Average SSH and its variance and EKE over high-SSH-variance region
% xmin_avg = -35.75;
% xmax_avg = -35.25;
% ymin_avg = -74.45;
% ymax_avg = -74.35;
xmin_avg = -37;
xmax_avg = -34;
ymin_avg = -74.75;
ymax_avg = -74.05;
xidx = find(XC(:,1)>xmin_avg & XC(:,1)<xmax_avg);
yidx = find(YC(1,:)>ymin_avg & YC(1,:)<ymax_avg);
ssh_mon_var_avg = squeeze(sum(sum(ssh_mon_var(xidx,yidx,startidx_ssh:end).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);
ssh_mon_std_avg = squeeze(sum(sum(ssh_mon_std(xidx,yidx,startidx_ssh:end).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);
eke_mon_avg = squeeze(sum(sum(EKE(xidx,yidx,startidx_ssh:end).*DXG(xidx,yidx).*DYG(xidx,yidx).*(Hocean(xidx,yidx)),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx).*(Hocean(xidx,yidx)),1),2);
eke_prod_mon_avg = squeeze(sum(sum(PEtoEKE_zavg(xidx,yidx,startidx_ssh:end).*DXG(xidx,yidx).*DYG(xidx,yidx).*(Hocean(xidx,yidx)),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx).*(Hocean(xidx,yidx)),1),2);
eke_surf_mon_avg = squeeze(sum(sum(EKE_surf(xidx,yidx,startidx_ssh:end).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);
eke_50m_mon_avg = squeeze(sum(sum(EKE_50m(xidx,yidx,startidx_ssh:end).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);

%%% Apply smoothing
smoothlen = 4;
[DSW_flux_smooth,times_smooth] = runningmean(times,T_ft(didx,:)',smoothlen);
[ssh_mon_var_avg_smooth,times_smooth] = runningmean(times,ssh_mon_var_avg,smoothlen);
[ssh_mon_std_avg_smooth,times_smooth] = runningmean(times,ssh_mon_std_avg,smoothlen);
[eke_mon_avg_smooth,times_smooth] = runningmean(times,eke_mon_avg,smoothlen);
[eke_prod_mon_avg_smooth,times_smooth] = runningmean(times,eke_prod_mon_avg,smoothlen);
[eke_surf_mon_avg_smooth,times_smooth] = runningmean(times,eke_surf_mon_avg,smoothlen);
[eke_50m_mon_avg_smooth,times_smooth] = runningmean(times,eke_50m_mon_avg,smoothlen);


%%% Locations at which correlations are maximized
[iidx_test,jidx_test] = find(corr_test==max(corr_test(:)));
[iidx_std,jidx_std] = find(corr_std==max(corr_std(:)));
[iidx_EKE,jidx_EKE] = find(corr_EKE==max(corr_EKE(:)));

%%% Sample time series from where correlations are maximized
ssh_mon_test_samp = squeeze((ssh_mon_test(iidx_test,jidx_test,startidx_ssh:end)));
ssh_mon_std_samp = squeeze((ssh_mon_std(iidx_test,jidx_test,startidx_ssh:end)));
rfac_test = ssh_mon_test_samp \ (T_ft(didx,:)'/1e6);
rfac_std = ssh_mon_std_samp \ (T_ft(didx,:)'/1e6);
ssh_mon_test_fit = rfac_test*ssh_mon_test_samp;
ssh_mon_std_fit = rfac_std*ssh_mon_std_samp;
RMSE_std = sqrt(mean((ssh_mon_std_fit-T_ft(didx,:)'/1e6).^2));
RMSE_test = sqrt(mean((ssh_mon_test_fit-T_ft(didx,:)'/1e6).^2));

%%% Linear regression of transport on SSH and EKE, averaged in the high-SSH-variance region
[Toff_ssh_std_avg,rfac_ssh_std_avg,rfac2_ssh_std_avg] = linearfits(ssh_mon_std_avg,T_ft(didx,:)');
[Toff_ssh_var_avg,rfac_ssh_var_avg,rfac2_ssh_var_avg] = linearfits(ssh_mon_var_avg,T_ft(didx,:)');
[Toff_eke_avg,rfac_eke_avg,rfac2_eke_avg] = linearfits(eke_mon_avg,T_ft(didx,:)');
[Toff_eke_surf_avg,rfac_eke_surf_avg,rfac2_eke_surf_avg] = linearfits(eke_surf_mon_avg,T_ft(didx,:)');
[Toff_eke_50m_avg,rfac_eke_50m_avg,rfac2_eke_50m_avg] = linearfits(eke_50m_mon_avg,T_ft(didx,:)');
[Toff_eke_surf_avg,rfac_eke_surf_avg,rfac2_eke_surf_avg] = linearfits(eke_surf_mon_avg,T_ft(didx,:)');


[Toff_ssh_std_avg_smooth,rfac_ssh_std_avg_smooth,rfac2_ssh_std_avg_smooth] = linearfits(ssh_mon_std_avg_smooth,DSW_flux_smooth);
[Toff_ssh_var_avg_smooth,rfac_ssh_var_avg_smooth,rfac2_ssh_var_avg_smooth] = linearfits(ssh_mon_var_avg_smooth,DSW_flux_smooth);
[Toff_eke_avg_smooth,rfac_eke_avg_smooth,rfac2_eke_avg_smooth] = linearfits(eke_mon_avg_smooth,DSW_flux_smooth);
[Toff_eke_surf_avg_smooth,rfac_eke_surf_avg_smooth,rfac2_eke_surf_avg_smooth] = linearfits(eke_surf_mon_avg_smooth,DSW_flux_smooth);
[Toff_eke_50m_avg_smooth,rfac_eke_50m_avg_smooth,rfac2_eke_50m_avg_smooth] = linearfits(eke_50m_mon_avg_smooth,DSW_flux_smooth);
[Toff_eke_prod_avg_smooth,rfac_eke_prod_avg_smooth,rfac2_eke_prod_avg_smooth] = linearfits(eke_prod_mon_avg_smooth,DSW_flux_smooth);

%%% Plotting options
fontsize = 14;
fignum = 0;
framepos = [1000         536         809         678];

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(eke_prod_mon_avg,T_ft(didx,:)/1e6);
hold on;
plot([0 3e-7],(Toff_eke_prod_avg+rfac_eke_prod_avg*[0 3e-7])/1e6,'k:');
plot([0 3e-7],rfac2_eke_prod_avg*[0 3e-7]/1e6,'k--');
hold off;
xlabel('Volume-averaged monthly EKE production (m^2/s^3)');
ylabel('DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
axis([0 3e-7 0 2]);
print('-dpng','-r150','Figures/Filchner/EKEprodvsDSWflux_scatter.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(eke_50m_mon_avg,T_ft(didx,:)/1e6);
hold on;
plot([0 0.08],(Toff_eke_50m_avg+rfac_eke_50m_avg*[0 0.08])/1e6,'k:');
plot([0 0.08],rfac2_eke_50m_avg*[0 0.08]/1e6,'k--');
hold off;
xlabel('Area-averaged monthly top-50m EKE (m^2/s^2)');
ylabel('DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
axis([0 0.08 0 2]);
print('-dpng','-r150','Figures/Filchner/EKE50mvsDSWflux_scatter.png');

% fignum = fignum + 1;
% figure(fignum);
% clf;
% set(gcf,'Position',framepos);
% scatter(eke_surf_mon_avg,T_ft(didx,:)/1e6);
% hold on;
% plot([0 0.035],(Toff_eke_surf_avg+rfac_eke_surf_avg*[0 0.035])/1e6,'k:');
% plot([0 0.035],rfac2_eke_surf_avg*[0 0.035]/1e6,'k--');
% hold off;
% xlabel('Area-averaged monthly surface EKE (m^2/s^2)');
% ylabel('DSW overflow transport (Sv)');
% set(gca,'FontSize',fontsize);
% axis([0 0.035 0 2]);
% print('-dpng','-r150','Figures/Filchner/EKEsurfvsDSWflux_scatter.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(ssh_mon_std_avg,T_ft(didx,:)/1e6);
hold on;
plot([0 0.08],(Toff_ssh_std_avg+rfac_ssh_std_avg*[0 0.08])/1e6,'k:');
plot([0 0.08],rfac2_ssh_std_avg*[0 0.08]/1e6,'k--');
hold off;
xlabel('Area-averaged monthly SSH std. (m)');
ylabel('DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
axis([0 0.08 0 2]);
print('-dpng','-r150','Figures/Filchner/SSHSTDvsDSWflux_scatter.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(ssh_mon_var_avg,T_ft(didx,:)/1e6);
hold on;
plot([0 0.003],(Toff_ssh_var_avg+rfac_ssh_var_avg*[0 0.003])/1e6,'k:');
plot([0 0.003],rfac2_ssh_var_avg*[0 0.003]/1e6,'k--');
hold off;
xlabel('Area-averaged monthly SSH var. (m)');
ylabel('DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
axis([0 0.003 0 2]);
print('-dpng','-r150','Figures/Filchner/SSHVARvsDSWflux_scatter.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(eke_mon_avg,T_ft(didx,:)/1e6);
hold on;
plot([0 0.02],(Toff_eke_avg+rfac_eke_avg*[0 0.02])/1e6,'k:');
plot([0 0.02],rfac2_eke_avg*[0 0.02]/1e6,'k--');
hold off;
xlabel('Volume-averaged monthly EKE (m^2/s^2)');
ylabel('DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
axis([0 0.02 0 2]);

print('-dpng','-r150','Figures/Filchner/EKEvsDSWflux_scatter.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(ssh_mon_std_avg,eke_mon_avg);
% hold on;
% plot(eke_mon_avg_sort,(Toff_eke_avg+rfac_eke_avg*eke_mon_avg_sort)/1e6,'k--');
% hold off;
ylabel('Volume-averaged monthly EKE (m^2/s^2)');
xlabel('Area-averaged monthly SSH std. (m)');
set(gca,'FontSize',fontsize);
% axis([0 6e-3 0 0.2]);
axis([0 0.08 0 0.035]);
print('-dpng','-r150','Figures/Filchner/EKEvsDSWflux_scatter.png');


%%% Smoothed correlations

xlim = [0 2e-3];
fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(ssh_mon_var_avg_smooth,DSW_flux_smooth/1e6);
hold on;
plot(xlim,(Toff_ssh_var_avg_smooth+rfac_ssh_var_avg_smooth*xlim)/1e6,'k:');
plot(xlim,rfac2_ssh_var_avg_smooth*xlim/1e6,'k--');
hold off;
xlabel('Smoothed area-averaged SSH var. (m)');
ylabel('Smoothed DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
set(gca,'XLim',xlim);
set(gca,'YLim',[0 2]);
print('-dpng','-r150','Figures/Filchner/SSHVARvsDSWflux_smoothed_scatter.png');

xlim = [0 0.08];
fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(ssh_mon_std_avg_smooth,DSW_flux_smooth/1e6);
hold on;
plot(xlim,(Toff_ssh_std_avg_smooth+rfac_ssh_std_avg_smooth*xlim)/1e6,'k:');
plot(xlim,rfac2_ssh_std_avg_smooth*xlim/1e6,'k--');
hold off;
xlabel('Smoothed area-averaged SSH std. (m)');
ylabel('Smoothed DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
set(gca,'XLim',xlim);
set(gca,'YLim',[0 2]);
print('-dpng','-r150','Figures/Filchner/SSHSTDvsDSWflux_smoothed_scatter.png');

xlim = [0 0.02];
fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(eke_mon_avg_smooth,DSW_flux_smooth/1e6);
hold on;
plot(xlim,(Toff_eke_avg_smooth+rfac_eke_avg_smooth*xlim)/1e6,'k:');
plot(xlim,rfac2_eke_avg_smooth*xlim/1e6,'k--');
hold off;
xlabel('Smoothed volume-averaged monthly EKE (m^2/s^2)');
ylabel('DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
set(gca,'XLim',xlim);
set(gca,'YLim',[0 2]);
print('-dpng','-r150','Figures/Filchner/EKEvsDSWflux_smoothed_scatter.png');

xlim = [0 0.05];
fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(eke_50m_mon_avg_smooth,DSW_flux_smooth/1e6);
hold on;
plot(xlim,(Toff_eke_50m_avg_smooth+rfac_eke_50m_avg_smooth*xlim)/1e6,'k:');
plot(xlim,rfac2_eke_50m_avg_smooth*xlim/1e6,'k--');
hold off;
xlabel('Smoothed volume-averaged monthly upper-50m EKE (m^2/s^2)');
ylabel('DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
set(gca,'XLim',xlim);
set(gca,'YLim',[0 2]);
print('-dpng','-r150','Figures/Filchner/EKE50mvsDSWflux_smoothed_scatter.png');

xlim = [0 2e-7];
fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
scatter(eke_prod_mon_avg_smooth,DSW_flux_smooth/1e6);
hold on;
plot(xlim,(Toff_eke_prod_avg_smooth+rfac_eke_prod_avg_smooth*xlim)/1e6,'k:');
plot(xlim,rfac2_eke_prod_avg_smooth*xlim/1e6,'k--');
hold off;
xlabel('Smoothed volume-averaged monthly EKE prod. (m^2/s^3)');
ylabel('DSW overflow transport (Sv)');
set(gca,'FontSize',fontsize);
set(gca,'XLim',xlim);
set(gca,'YLim',[0 2]);
print('-dpng','-r150','Figures/Filchner/EKEprodvsDSWflux_smoothed_scatter.png');

% figure(9);
% scatter(ssh_mon_test_fit,T_ft(didx,:)/1e6);
% hold on;
% plot(sort(ssh_mon_test_fit),sort(ssh_mon_test_fit),'k--')
% hold off;
% axis([0 2 0 2])

% figure(10);
% scatter(squeeze((ssh_mon_std(iidx_std,jidx_std,startidx_ssh:end))),T_ft(didx,:)/1e6);
% axis([0 0.1 0 2])
% 
% figure(11);
% clf;
% plot(times/t1year,T_ft(didx,:)/1e6,'r');
% colororder = get(gca,'ColorOrder');
% box off;
% set(gca,'YLim',[0 2]);
% ax = gca;
% ax2 = axes('Position',get(ax,'Position'));
% plot(ax2,times/t1year,rfac_std*squeeze(ssh_mon_std(iidx_std,jidx_std,startidx_ssh:end)),'g');
% set(ax2,'Color','None');
% set(ax2,'YAxisLocation','Right');
% set(ax2,'YLim',[0 2]);
% box off;


DSW_flux = T_ft(didx,:)';
R2_ssh_var = 1 - sum((DSW_flux-rfac2_ssh_var_avg*ssh_mon_var_avg).^2) / sum((DSW_flux-mean(DSW_flux)).^2);
R2_eke = 1 - sum((DSW_flux-rfac2_eke_avg*eke_mon_avg).^2) / sum((DSW_flux-mean(DSW_flux)).^2);
fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
plot(times/t1year,T_ft(didx,:)/1e6);
colororder = get(gca,'ColorOrder');
box off;
set(gca,'YLim',[0 2]);
ax = gca;
hold on;
plot(times/t1year,(rfac2_ssh_var_avg*ssh_mon_var_avg)/1e6);
plot(times/t1year,(rfac2_eke_avg*eke_mon_avg)/1e6);
hold off;
set(ax,'FontSize',fontsize);
xlabel(ax,'Time (years)');
ylabel(ax,'DSW transport, monthly mean (Sv)');
legend('Diagnosed','SSH var. linear regression','Volume-averaged EKE linear regression')
text(4,0.4,['R^2 = ',num2str(R2_ssh_var)],'Color',colororder(2,:));
text(4,0.3,['R^2 = ',num2str(R2_eke)],'Color',colororder(3,:));
print('-dpng','-r150','Figures/Filchner/TimeSeries_monthly.png');


R2_ssh_var_smooth = 1 - sum((DSW_flux_smooth-rfac2_ssh_var_avg_smooth*ssh_mon_var_avg_smooth).^2) / sum((DSW_flux_smooth-mean(DSW_flux_smooth)).^2);
R2_eke_smooth = 1 - sum((DSW_flux_smooth-rfac2_eke_avg_smooth*eke_mon_avg_smooth).^2) / sum((DSW_flux_smooth-mean(DSW_flux_smooth)).^2);
fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
plot(times_smooth/t1year,DSW_flux_smooth/1e6);
colororder = get(gca,'ColorOrder');
box off;
set(gca,'YLim',[0 2]);
ax = gca;
hold on;
plot(times_smooth/t1year,(rfac2_ssh_var_avg_smooth*ssh_mon_var_avg_smooth)/1e6);
plot(times_smooth/t1year,(rfac2_eke_avg_smooth*eke_mon_avg_smooth)/1e6);
hold off;
set(ax,'FontSize',fontsize);
xlabel(ax,'Time (years)');
ylabel(ax,'DSW transport, smoothed (Sv)');
legend('Diagnosed','SSH var. linear regression','Volume-averaged EKE linear regression')
text(4,0.4,['R^2 = ',num2str(R2_ssh_var_smooth)],'Color',colororder(2,:));
text(4,0.3,['R^2 = ',num2str(R2_eke_smooth)],'Color',colororder(3,:));
print('-dpng','-r150','Figures/Filchner/TimeSeries_smoothed.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,corr_var);
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
plot([xmin_avg xmax_avg xmax_avg xmin_avg xmin_avg],[ymin_avg ymin_avg ymax_avg ymax_avg ymin_avg],'k--','LineWidth',1.5);
hold off;
caxis([-1 1]);
colormap(cmocean('balance',20));
title('Correlation between Filchner overflow and SSH var.')
print('-dpng','-r150','Figures/Filchner/CorrMap_SSHVAR.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,corr_mean);
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
plot([xmin_avg xmax_avg xmax_avg xmin_avg xmin_avg],[ymin_avg ymin_avg ymax_avg ymax_avg ymin_avg],'k--','LineWidth',1.5);
hold off;
caxis([-1 1]);
colormap(cmocean('balance',20));
title('Correlation between Filchner overflow and mean SSH')
print('-dpng','-r150','Figures/Filchner/CorrMap_SSHMEAN.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,corr_std);
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
plot([xmin_avg xmax_avg xmax_avg xmin_avg xmin_avg],[ymin_avg ymin_avg ymax_avg ymax_avg ymin_avg],'k--','LineWidth',1.5);
hold off;
caxis([-1 1]);
colormap(cmocean('balance',20));
title('Correlation between Filchner overflow and SSH std.')
print('-dpng','-r150','Figures/Filchner/CorrMap_SSHSTD.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,corr_EKE);
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
plot([xmin_avg xmax_avg xmax_avg xmin_avg xmin_avg],[ymin_avg ymin_avg ymax_avg ymax_avg ymin_avg],'k--','LineWidth',1.5);
hold off;
caxis([-1 1]);
colormap(cmocean('balance',20));
title('Correlation between Filchner overflow and depth-averaged EKE.')
print('-dpng','-r150','Figures/Filchner/CorrMap_EKE.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,corr_EKE_50m);
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
plot([xmin_avg xmax_avg xmax_avg xmin_avg xmin_avg],[ymin_avg ymin_avg ymax_avg ymax_avg ymin_avg],'k--','LineWidth',1.5);
hold off;
caxis([-1 1]);
colormap(cmocean('balance',20));
title('Correlation between Filchner overflow and upper-50m EKE')
print('-dpng','-r150','Figures/Filchner/CorrMap_EKE50m.png');



fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,log10(mean(EKE,3)));
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
plot([xmin_avg xmax_avg xmax_avg xmin_avg xmin_avg],[ymin_avg ymin_avg ymax_avg ymax_avg ymin_avg],'k--','LineWidth',1.5);
hold off;
% caxis([0 1]*3.5e-2);
caxis([-4 -1]);
colormap(cmocean('amp',100));
title('Depth-averaged EKE (log_1_0 m^2/s^2)')
print('-dpng','-r150','Figures/Filchner/EKE.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,mean(ssh_mon_var,3));
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
plot([xmin_avg xmax_avg xmax_avg xmin_avg xmin_avg],[ymin_avg ymin_avg ymax_avg ymax_avg ymin_avg],'k--','LineWidth',1.5);
hold off;
caxis([0 1]*3e-3);
colormap(cmocean('amp',100));
title('Mean monthly SSH variance (m^2)')
print('-dpng','-r150','Figures/Filchner/SSHVAR.png');

fignum = fignum + 1;
figure(fignum);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,mean(ssh_mon_std,3));
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
plot([xmin_avg xmax_avg xmax_avg xmin_avg xmin_avg],[ymin_avg ymin_avg ymax_avg ymax_avg ymin_avg],'k--','LineWidth',1.5);
hold off;
caxis([0 1]*0.08);
colormap(cmocean('amp',8));
title('Mean monthly standard deviation of SSH (m)')
print('-dpng','-r150','Figures/Filchner/SSHSTD.png');






% figure(13);
% scatter(T_ft_smooth(didx,1:end-lagidx),ssh_var_runmean_avg(1+lagidx:end));
% scatter(T_ft_smooth(didx,1:end-2),sqrt(ssh_var_runmean_avg(1+2:end)));


% % [r,lags] = xcorr(Tdsw_10day,ssh_var_runmean_avg)
% [r,lags] = xcorr(T_ft_smooth(didx,:)',ssh_var_runmean_avg)
% figure(16);
% plot(lags,r)
% 
% rr = zeros(length(dens_levs),1);
% for m=1:length(rr)
%   rr(m) = corr(T_ft_smooth(m,1:end-lagidx)',ssh_var_runmean_avg(1+lagidx:end));
% end
% figure(17);
% plot(dens_levs,rr);


%%% Convenience function to compute running means
function [vec_smooth,times_smooth] = runningmean (times,vec,smoothlen)

  %%% To store output
  times_smooth = times;
  vec_smooth = vec;

  %%% Compute running means
  for n=smoothlen+1:1:length(vec)-smoothlen
    vec_smooth(n) = mean(vec(n-smoothlen:n+smoothlen));
  end

  %%% Remove end points where running mean can't be computed
  times_smooth(length(vec)-smoothlen:length(vec)) = [];  
  vec_smooth(length(vec)-smoothlen:length(vec)) = [];
  times_smooth(1:smoothlen) = [];
  vec_smooth(1:smoothlen) = [];

end

%%% Convenience function to compute linear fits
function [Toff,rfac,rfac2] = linearfits (vec2,vec1)

  [vec2_sort,I] = sort(vec2);
  vec1_sort = vec1(I);
  p = polyfit(vec2_sort,vec1_sort,1);
  rfac = p(1);
  Toff = p(2);
  rfac2 = squeeze((vec2)) \ (vec1);

end