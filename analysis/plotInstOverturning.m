%%%
%%% calcInstOverturning.m
%%%
%%% Calculates the overturning circulation in density surfaces, similar to 
%%% that calculated using the MITgcm 'layers' package, from instantaneous
%%% model output
%%%


%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
tmin = 1.01;
tmax = 7.01;

%%% Load experiment
loadexp;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Select density variable in which to compute isopycnal fluxes
densvar = 'PD0';
% densvar = 'ND1';
% densvar = 'ND2';
% densvar = 'PT';

%%% Load pre-computed MOC
outfname = [expname,'_MOC_inst_',densvar];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

%%% Load pre-computed Filchner Trough transport
outfname = [expname,'_FTtrans_inst_',densvar];
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

%%% Load pre-computed SSH variance
ncfname = fullfile('products',[expname,'_SSHvar.nc']);
ssh_var_runmean = ncread(ncfname,'ssh_var_runmean');
ssh_runmean = ncread(ncfname,'ssh_runmean');


%%% Mean streamfunction
psi_plot = mean(psi,3)/1e6;
ylim = [27.3 28.2];
psistep = 0.25;
fontsize = 14;
psimax = 4;
figure(10);
clf;
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
hold on;
[C,h] = contour(EE,DD,psi_plot,[-6 -5 -4 -3 -2 -1 -0.5 -0.25 -.1 .1 0.25 0.5 1],'EdgeColor','k');
clabel(C,h);
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
plot(eta,27.85*ones(size(eta)),'k:')
hold off;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 10.5]);
% grid on;
colormap(gca,cmocean('balance',round(2*psimax/(psistep))));%,'pivot',0));
xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)')
set(gca,'FontSize',fontsize);
cbhandle = colorbar;
% set(cbhandle,'Position',cb3_pos);
title(cbhandle,'Sv');
text(-6.5,27.4,'Overturning streamfunction, \psi^\sigma','FontSize',fontsize+2);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);

%%% Variability of streamfunction
psi_plot = std(psi,[],3)/1e6;
ylim = [27.3 28.2];
psistep = 0.25;
fontsize = 14;
psimax = 4;
figure(11);
clf;
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
hold on;
[C,h] = contour(EE,DD,psi_plot,[-6 -5 -4 -3 -2 -1 -0.5 -0.25 -.1 .1 0.25 0.5 1],'EdgeColor','k');
clabel(C,h);
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
plot(eta,27.85*ones(size(eta)),'k:')
hold off;
caxis([0 psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 10.5]);
% grid on;
colormap(gca,cmocean('amp',round(psimax/(psistep))));%,'pivot',0));
xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)')
set(gca,'FontSize',fontsize);
cbhandle = colorbar;
% set(cbhandle,'Position',cb3_pos);
title(cbhandle,'Sv');
text(-6.5,27.4,'std(\psi^\sigma)','FontSize',fontsize+2);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);

% eidx1 = find(eta>=4.3,1);
% eidx2 = find(eta>=4.3,1);
% didx = find(dens_levs>=27.85,1);
eidx1 = find(eta>=3.6,1);
eidx2 = find(eta>=3.6,1);
didx = find(dens_levs>=27.81,1);
Tdsw = -squeeze(mean(psi(eidx1:eidx2,didx,:),1))/1e6;
Tdsw_smooth = 0*Tdsw;
T_ft_smooth = 0*T_ft;
smooth_len = 10;
lagidx = 6;
Nt = length(Tdsw_smooth);
for m=1:Nt
  Tdsw_smooth(m) = mean(Tdsw(max(m-smooth_len,1):min(m+smooth_len,Nt)));
  T_ft_smooth(:,m) = mean(T_ft(:,max(m-smooth_len,1):min(m+smooth_len,Nt)),2);
end


%%% Average SSH and its variance over high-SSH region
xidx = find(XC(:,1)>-36.4 & XC(:,1)<-33.1);
yidx = find(YC(1,:)>-74.5 & YC(1,:)<-74.2);
% xidx = find(XC(:,1)>-36.4 & XC(:,1)<-33.1);
% yidx = find(YC(1,:)>-74.6 & YC(1,:)<-74.1);
% xidx = find(XC(:,1)>-40 & XC(:,1)<-30);
% yidx = find(YC(1,:)>-75 & YC(1,:)<-74);
ssh_var_runmean_avg = squeeze(sum(sum(ssh_var_runmean(xidx,yidx,:).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);
ssh_runmean_avg = squeeze(sum(sum(ssh_runmean(xidx,yidx,:).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);


figure(12);
clf;
% plot(-smooth(squeeze(mean(psi(eidx1:eidx2,didx,:),1)),20)/1e6);
tt = ((1:Nt)-1)/2;
plot(tt,Tdsw_smooth);
colororder = get(gca,'ColorOrder');
hold on;
plot(tt,T_ft_smooth(didx,:)/1e6,'r');
hold off;
box off;
ax = gca;
ax2 = axes('Position',get(ax,'Position'));
plot(ax2,tt,ssh_var_runmean_avg,'g');
set(ax2,'Color','None');
set(ax2,'YAxisLocation','Right');
box off;
% set(ax(2),'YLim',[0 0.07]);
% plotyy(tt,Tdsw,tt,sqrt(ssh_var_runmean_avg));
% plotyy(tt,Tdsw_10day,tt,smooth(sqrt(ssh_var_runmean_avg),20));
text(ax,10,3,['r=',num2str(corr(Tdsw_smooth,ssh_var_runmean_avg))],'FontSize',fontsize);
text(ax,10,2.5,['r=',num2str(corr(T_ft_smooth(didx,:)',ssh_var_runmean_avg))],'FontSize',fontsize);
text(ax,10,2,['r=',num2str(corr(T_ft_smooth(didx,1:end-lagidx)',ssh_var_runmean_avg(1+lagidx:end)))],'FontSize',fontsize);

set(ax,'FontSize',fontsize);
set(ax2,'FontSize',fontsize);
xlabel(ax,'2011 yearday');
ylabel(ax,'DSW transport, 15-day mean (Sv)');
ylabel(ax2,'10-day SSH variance (m^2)');


figure(13);
scatter(T_ft_smooth(didx,1:end-lagidx),ssh_var_runmean_avg(1+lagidx:end));
scatter(T_ft_smooth(didx,1:end-2),sqrt(ssh_var_runmean_avg(1+2:end)));


% [r,lags] = xcorr(Tdsw_10day,ssh_var_runmean_avg)
[r,lags] = xcorr(T_ft_smooth(didx,:)',ssh_var_runmean_avg)
figure(16);
plot(lags,r)

rr = zeros(length(dens_levs),1);
for m=1:length(rr)
  rr(m) = corr(T_ft_smooth(m,1:end-lagidx)',ssh_var_runmean_avg(1+lagidx:end));
end
figure(17);
plot(dens_levs,rr);

figure(14);

clf;
% plot(-smooth(squeeze(mean(psi(eidx1:eidx2,didx,:),1)),20)/1e6);
tt = ((1:Nt)-1)/2;
ax = plotyy(tt,Tdsw_smooth,tt,sqrt(ssh_runmean_avg));
% plotyy(tt,Tdsw,tt,sqrt(ssh_var_runmean_avg));
% plotyy(tt,Tdsw_10day,tt,smooth(sqrt(ssh_var_runmean_avg),20));

figure(15);
scatter(T_ft_smooth(didx,:),ssh_runmean_avg);



figure(18);
clf;
tt = ((1:Nt)-1)/2;
plot(tt(1:end-lagidx),T_ft_smooth(didx,1:end-lagidx)/1e6);
colororder = get(gca,'ColorOrder');
box off;
ax = gca;
ax2 = axes('Position',get(ax,'Position'));
plot(ax2,tt(1:end-lagidx),ssh_var_runmean_avg(1+lagidx:end)*1e4,'Color',colororder(2,:));
set(ax2,'Color','None');
set(ax2,'YAxisLocation','Right');
box off;
text(ax,10,2,['r=',num2str(corr(T_ft_smooth(didx,1:end-lagidx)',ssh_var_runmean_avg(1+lagidx:end)))],'FontSize',fontsize);
set(ax,'FontSize',fontsize);
set(ax2,'FontSize',fontsize);
set(ax,'YLim',[0.3 2.5]);
set(ax2,'YLim',[0 70]);
xlabel(ax,'2011 yearday');
ylabel(ax,'DSW transport, 10-day mean (Sv)');
ylabel(ax2,'10-day SSH variance (cm^2)');
title(['3-day lag between Filchner Trough transport (latitude ',num2str(ft_lat),') and SSH variance']);



