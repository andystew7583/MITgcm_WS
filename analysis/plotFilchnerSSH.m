%%%
%%% plotFilchnerSSH.m
%%%
%%% Plots variations of Filchner overflow transport and SSH variance.
%%%


%%% Options
expdir = '../experiments';
expname = 'hires_nest_onethirtieth_notides_RTOPO2';
tmin = 1.01;
tmax = 7.01;

%%% Load experiment
loadexp;

%%% Select density variable in which to compute isopycnal fluxes
densvar = 'PD0';
% densvar = 'ND1';
% densvar = 'ND2';
% densvar = 'PT';

%%% Month index to start reading SSH data (first year is skipped in DSW
%%% flux)
startidx_ssh = 13;

outfname = [expname,'_ShelfHeatBudget.mat'];
load(fullfile('products',outfname),'usq_eddy_int','vsq_eddy_int');
EKE = 0.5*(usq_eddy_int+vsq_eddy_int);

%%% Load pre-computed Filchner Trough transport
outfname = [expname,'_FTtrans_',densvar];
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));
times_tmp = [366*86400 times];
times = 0.5*(times_tmp(1:end-1)+times_tmp(2:end));

%%% Index of density level corresponding to max transport
didx = find(dens_levs>=27.81,1);

%%% Load pre-computed SSH variance
ncfname = fullfile('products',[expname,'_SSHvar.nc']);
ssh_var_runmean = ncread(ncfname,'ssh_var_runmean');
ssh_runmean = ncread(ncfname,'ssh_runmean');

%%% Load monthly-mean SSH variance
ncfname = fullfile('products',[expname,'_SSHmonVar.nc']);
ssh_mon_mean = ncread(ncfname,'ssh_mon_mean');
ssh_mon_var = ncread(ncfname,'ssh_mon_var');
ssh_mon_var(ssh_mon_var<0) = 0; %%% Sometimes <0 to machine precision
ssh_mon_std = sqrt(ssh_mon_var);
ssh_mon_test = ssh_mon_var.^(2/3);

%%% Compute correlations
corr_mean = zeros(Nx,Ny);
corr_var = zeros(Nx,Ny);
corr_std = zeros(Nx,Ny);
corr_EKE = zeros(Nx,Ny);

corr_test = zeros(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    corr_mean(i,j) = corr(squeeze(ssh_mon_mean(i,j,startidx_ssh:end)),T_ft(didx,:)');
    corr_var(i,j) = corr(squeeze(ssh_mon_var(i,j,startidx_ssh:end)),T_ft(didx,:)');
    corr_std(i,j) = corr(squeeze(ssh_mon_std(i,j,startidx_ssh:end)),T_ft(didx,:)');

    corr_test(i,j) = corr(squeeze(ssh_mon_test(i,j,startidx_ssh:end)),T_ft(didx,:)');
    corr_EKE(i,j) = corr(squeeze(EKE(i,j,startidx_ssh:end)),T_ft(didx,:)');
  end
end



%%% Average SSH and its variance over high-SSH region
xidx = find(XC(:,1)>-35.75 & XC(:,1)<-35.25);
yidx = find(YC(1,:)>-74.45 & YC(1,:)<-74.35);
ssh_mon_var_avg = squeeze(sum(sum(ssh_mon_var(xidx,yidx,:).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);
ssh_mon_std_avg = squeeze(sum(sum(ssh_mon_std(xidx,yidx,:).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);




[iidx_test,jidx_test] = find(corr_test==max(corr_test(:)));
[iidx_std,jidx_std] = find(corr_std==max(corr_std(:)));
[iidx_EKE,jidx_EKE] = find(corr_EKE==max(corr_EKE(:)));

ssh_mon_test_samp = squeeze((ssh_mon_test(iidx_test,jidx_test,startidx_ssh:end)));

ssh_mon_std_samp = squeeze((ssh_mon_std(iidx_test,jidx_test,startidx_ssh:end)));

rfac_test = ssh_mon_test_samp \ (T_ft(didx,:)'/1e6);
rfac_std = ssh_mon_std_samp \ (T_ft(didx,:)'/1e6);
rfac_avg = squeeze((ssh_mon_std_avg(startidx_ssh:end))) \ (T_ft(didx,:)'/1e6);

ssh_mon_test_fit = rfac_test*ssh_mon_test_samp;

ssh_mon_std_fit = rfac_std*ssh_mon_std_samp;

RMSE_std = sqrt(mean((ssh_mon_std_fit-T_ft(didx,:)'/1e6).^2))
RMSE_test = sqrt(mean((ssh_mon_test_fit-T_ft(didx,:)'/1e6).^2))

fontsize = 14;

figure(9);
scatter(ssh_mon_test_fit,T_ft(didx,:)/1e6);
hold on;
plot(sort(ssh_mon_test_fit),sort(ssh_mon_test_fit),'k--')
hold off;
axis([0 2 0 2])

figure(10);
scatter(squeeze((ssh_mon_std(iidx_std,jidx_std,startidx_ssh:end))),T_ft(didx,:)/1e6);
axis([0 0.1 0 2])

figure(11);
clf;
plot(times/t1year,T_ft(didx,:)/1e6,'r');
colororder = get(gca,'ColorOrder');
box off;
set(gca,'YLim',[0 2]);
ax = gca;
ax2 = axes('Position',get(ax,'Position'));
plot(ax2,times/t1year,rfac_std*squeeze(ssh_mon_std(iidx_std,jidx_std,startidx_ssh:end)),'g');
set(ax2,'Color','None');
set(ax2,'YAxisLocation','Right');
set(ax2,'YLim',[0 2]);
box off;

figure(12);
clf;
plot(times/t1year,T_ft(didx,:)/1e6,'r');
colororder = get(gca,'ColorOrder');
box off;
set(gca,'YLim',[0 2]);
ax = gca;
ax2 = axes('Position',get(ax,'Position'));
plot(ax2,times/t1year,rfac_avg*ssh_mon_std_avg(startidx_ssh:end)','g');
set(ax2,'Color','None');
set(ax2,'YAxisLocation','Right');
box off;
% set(ax(2),'YLim',[0 0.07]);
% plotyy(tt,Tdsw,tt,sqrt(ssh_var_runmean_avg));
% plotyy(tt,Tdsw_10day,tt,smooth(sqrt(ssh_var_runmean_avg),20));
set(ax,'FontSize',fontsize);
set(ax2,'FontSize',fontsize);
set(ax2,'YLim',[0 2]);
xlabel(ax,'Time (years)');
ylabel(ax,'DSW transport, 15-day mean (Sv)');
ylabel(ax2,'10-day SSH variance (m^2)');
set(ax2,'Position',get(ax,'Position'));


figure(13);
pcolor(XC,YC,corr_var);
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
hold off;
caxis([-1 1]);

figure(14);
pcolor(XC,YC,corr_mean);
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
hold off;
caxis([-1 1]);

figure(15);
pcolor(XC,YC,std(ssh_mon_var,1,3));
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
hold off;
caxis([-1 1]*1e-3);


figure(16);
pcolor(XC,YC,corr_std);
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
hold off;
caxis([-1 1]);

figure(17);
pcolor(XC,YC,corr_EKE);
shading flat;
colorbar;
hold on;
contour(XC,YC,bathy,[-300 -2000 -1000 -500],'EdgeColor','k');
hold off;
caxis([-1 1]);




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


