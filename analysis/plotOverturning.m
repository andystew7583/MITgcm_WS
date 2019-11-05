 %%%
%%% plotOverturning.m
%%%
%%% Plots the overturning circulation in potential density space.
%%%

%%% Load experiment and pre-computed MOC data
loadexp;
load([expname,'_MOC_pt.mat']);

%%% Plotting options
mac_plots = 1;
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.15 0.17 0.68 0.78];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end
psimin = -1.5;
psimax = 1.5;

%%% Calculate zonal-mean density
pt_xtavg = squeeze(nanmean(pt_tavg(:,:,:)));
pt_f_xtavg = squeeze(nanmean(pt_f(:,:,:)));
vflux_xint = squeeze(nanmean(vflux(:,:,:)))*Lx;
vflux_m_xint = squeeze(nanmean(vflux_m(:,:,:)))*Lx;

%%% Overturning in y/pt space
figure(6);
clf;
axes('FontSize',16);
[PT YY] = meshgrid(ptlevs,yy);
contourf(YY,PT,psi_pt,30,'EdgeColor','None');
colorbar;
colormap redblue;
caxis([psimin -psimin]);

%%% Mean overturning in y/pt space
figure(7);
clf;
axes('FontSize',16);
[PT YY] = meshgrid(ptlevs,yy);
contourf(YY,PT,psim_pt,30);
colorbar;
colormap redblue;
caxis([psimin -psimin]);

%%% Eddy overturning in y/pt space
figure(8);
clf;
axes('FontSize',16);
[PT YY] = meshgrid(ptlevs,yy);
contourf(YY,PT,psie_pt,30);
colorbar;
colormap redblue;
caxis([psimin -psimin]);

%%% Isopycnal fluxes in y/pt space
figure(9);
clf;
axes('FontSize',16);
[PT YY] = meshgrid((ptlevs(2:end)+ptlevs(1:end-1))/2,yy);
contourf(YY,PT,vflux_xint,30);
colorbar;
colormap redblue;

%%% y/z grid for streamfunction plots
% makePsiGrid;
zz_psi = zz_f;
yy_psi = yy;
[ZZ_psi,YY_psi] = meshgrid(zz_psi,yy_psi);
for j=1:Ny      
  
  %%% Adjust height of top point
  ZZ_psi(j,1) = 0;
  
  %%% Calculate depth of bottom cell  
  hFacS_col = squeeze(hFacS_f(1,j,:));  
  kmax = length(hFacS_col(hFacS_col>0)); 
  if (kmax > 0)
    ZZ_psi(j,kmax) = - sum(delRf.*hFacS_col');
  end
  
end

%%% Plot the residual overturning in y/z space
handle = figure(10);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_z,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_z,[psimin:0.025:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Mean overturning in y/z space
figure(11);
clf;
axes('FontSize',16);
contourf(YY_psi,ZZ_psi,psim_z,30);
colorbar;
colormap redblue;
caxis([psimin -psimin]);

%%% Plot the eddy streamfunction in y/z space
handle = figure(12);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi,ZZ_psi,psie_z,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
% plot(yy/1000,-hb/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi,ZZ_psi,psie_z,[psimin:0.05:-0.05 0.05:0.05:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
if (mac_plots)
  annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{eddy}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
else
  annotation('textbox',[0.7 0.05 0.3 0.05],'String','$\psi_{\mathrm{eddy}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end

%%% Fine-grid potential temperature
figure(13);
clf;
axes('FontSize',16);
[ZZ YY] = meshgrid(zz_f,yy);
contourf(YY,ZZ,pt_f_xtavg,ptlevs);
colorbar;
colormap jet;

%%% Coarse-grid potential temperature
figure(14);
clf;
axes('FontSize',16);
[ZZ YY] = meshgrid(zz,yy);
contourf(YY,ZZ,pt_xtavg,ptlevs);
colorbar;
colormap jet;