%%%
%%% paper2_plotOverturning.m
%%%
%%% Plots the overturning circulation variability in density space.
%%%

%%% Needed for EOF
addpath CDT/cdt;

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_RTOPO2';
loadexp;

%%% Options (see calcOverturning)
calc_psi_eddy = true;
deform_cavity = false;
use_layers = true;
densvar = 'PD0';
psimax = 6;
psistep = 0.25;
% psimax = 2;
% psistep = 0.1;
ylim = [27.2 28.2];
% ylim = [27.3 29];

%%% Construct output file name
outfname = [expname,'_MOC_',densvar];
if (calc_psi_eddy)
  if (use_layers)
    estr = '_layers';
  else
    estr = '_TRM';
  end
else
  estr = '_noeddy';
end
outfname = [outfname,estr];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));
[DD,EE] = meshgrid(dens_levs,eta);

%%% Compute EOFs
Nt = length(times);
Neta = length(eta);
Nd = length(dens_levs);
psi_tot = psi_mean+psi_eddy;
psi_tot_noseason = psi_tot - repmat(squeeze(mean(reshape(psi_tot,[Neta Nd 12,size(psi_tot,3)/12]),4)),[1 1 size(psi_tot,3)/12]);
[eof_maps,pc,expvar] = eof((psi_tot)/1e6);
[eof_maps_noseason,pc_noseason,expvar_noseason] = eof((psi_tot_noseason)/1e6);
[eof_maps_offshelf,pc_offshelf,expvar_offshelf] = eof((psi_tot)/1e6,'mask',EE>=4);
[eof_maps_shelf,pc_shelf,expvar_shelf] = eof((psi_tot)/1e6,'mask',EE<4);
[eof_maps_cavity,pc_cavity,expvar_cavity] = eof((psi_tot)/1e6,'mask',EE<0);

%%% Load pre-computed time-mean T and S
outfname = [expname,'_TSfluxes'];
outfname = [outfname,'.mat'];
load(fullfile('products',outfname),'theta_tavg');
load(fullfile('products',outfname),'salt_tavg');

%%% Load surface flux time series
outfname = [expname,'_surfFluxes'];
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

%%% Mean "sea surface" salinity and temperature for flux calculation
sst = NaN*ones(Nx,Ny);
sss = NaN*ones(Nx,Ny);
ssp = NaN*ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    kmin = find(hFacC(i,j,:)>0,1,'first');
    if (~isempty(kmin))
      sss(i,j) = salt_tavg(i,j,kmin);
      sst(i,j) = theta_tavg(i,j,kmin);
      ssp = -gravity*rhoConst*RC(kmin)/1e4;
    end
  end
end
[alpha,beta] = calcAlphaBeta(sss,sst,ssp);
Cp = 4e3;

%%% Loop through outputs and compute buoyancy fluxes
bflux_int = zeros(1,Nt);
bflux_int_shelf = zeros(1,Nt);
bflux_int_offshelf = zeros(1,Nt);
taux_int = zeros(1,Nt);
taux_ews_int = zeros(1,Nt);
for n=1:Nt
  
  %%% Composite surface heat flux
  htFlxDn = tflux(:,:,n);
  tmp = SHIhtFlx(:,:,n);
  htFlxDn(hFacC(:,:,1)==0) = -tmp(hFacC(:,:,1)==0);
  htFlxDn(sum(hFacC,3)==0) = NaN;

  %%% Composite virtual surface salt flux
  sltFlxDn = sflux(:,:,n);
  tmp = SHIfwFlx(:,:,n);
  sltFlxDn(hFacC(:,:,1)==0) = tmp(hFacC(:,:,1)==0).*sss(hFacC(:,:,1)==0);
  sltFlxDn(sum(hFacC,3)==0) = NaN;

  %%% Calculate buoyancy flux
  bfluxT = gravity*(alpha.*htFlxDn/Cp/rhoConst);
  bfluxS = - gravity*(beta.*sltFlxDn/rhoConst);
  bflux = bfluxT + bfluxS;
  
  %%% Area-integrated buoyancy flux
  bflux_int(n) = nansum(nansum(bflux.*RAC));   
  bflux_int_shelf(n) = nansum(nansum(bflux.*RAC.*(ETA<4)));
  bflux_int_offshelf(n) = nansum(nansum(bflux.*RAC.*(ETA>=4)));
  
  %%% Area-averaged wind stress
  taux_int(n) = sum(sum(oceTAUX(:,:,n).*RAW.*hFacW(:,:,1)))./sum(sum(RAW.*hFacW(:,:,1)));
  msk_ews = (XG>-20) & (YC<-66);
  taux_ews_int(n) = sum(sum(oceTAUX(:,:,n).*RAW.*hFacW(:,:,1).*msk_ews))./sum(sum(RAW.*hFacW(:,:,1).*msk_ews));

end
bflux_int_noseason = bflux_int - repmat(squeeze(mean(reshape(bflux_int,[1 12 length(bflux_int)/12]),3)),[1 length(bflux_int)/12]);% + mean(bflux_int);
bflux_int_shelf_noseason = bflux_int_shelf - repmat(squeeze(mean(reshape(bflux_int_shelf,[1 12 length(bflux_int_shelf)/12]),3)),[1 length(bflux_int_shelf)/12]);% + mean(bflux_int_shelf);
bflux_int_offshelf_noseason = bflux_int_offshelf - repmat(squeeze(mean(reshape(bflux_int_offshelf,[1 12 length(bflux_int_offshelf)/12]),3)),[1 length(bflux_int_offshelf)/12]);% + mean(bflux_int_offshelf);
taux_int_noseason = taux_int - repmat(squeeze(mean(reshape(taux_int,[1 12 length(taux_int)/12]),3)),[1 length(taux_int)/12]);
taux_ews_int_noseason = taux_ews_int - repmat(squeeze(mean(reshape(taux_ews_int,[1 12 length(taux_ews_int)/12]),3)),[1 length(taux_ews_int)/12]);

%%% Streamfunction strength estimates
psimax_cavity = zeros(1,Nt);
psimax_shelf = zeros(1,Nt);
psimax_full = zeros(1,Nt);
for n=1:Nt
  tmp = psi_tot(:,:,n);
  psimax_full(n) = - min(min(tmp));
  psimax_shelf(n) = -min(min(tmp(EE<4)));
  psimax_cavity(n) = max(max(tmp(EE<0)));
end







%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE PLOT 1 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
                                
fontsize = 14;
framepos = [417    26   1030   959];
axpos = zeros(6,4);
axpos(1,:) = [0.08 0.68 0.38 0.25];
axpos(2,:) = [0.55 0.68 0.38 0.25];
axpos(3,:) = [0.08 0.36 0.38 0.25];
axpos(4,:) = [0.55 0.36 0.38 0.25];
axpos(5,:) = [0.08 0.05 0.38 0.25];
axpos(6,:) = [0.55 0.05 0.38 0.25];
cb1_pos = [0.94 0.68 0.015 0.25];
cb2_pos = [0.94 0.36 0.015 0.25];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

figure(207);                
clf;
set(gcf,'Position',framepos); 
colororder = get(gca,'ColorOrder');
                

%%% Streamfunction standard deviation
subplot('Position',axpos(1,:));
psi_plot = std((psi_tot)/1e6,[],3);
pcolor(EE,DD,psi_plot);
shading interp;
hold on
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
hold off
caxis([0 2.5]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 11]);
colormap(gca,cmocean('amp',25));
% colormap(gca,hot(25));
% xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)')
text(-6.5,27.3,'Overturning streamfunction s.d.','FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);
cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
title(cbhandle,'Sv');
handle = title(gca,'Seasonal cycle included','Fontsize',fontsize+2);
set(handle,'Position',get(handle,'Position')+[0 -.17 0]);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);

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











%%% Streamfunction standard deviation, seasonal cycle removed
subplot('Position',axpos(2,:));
psi_plot = std((psi_tot_noseason)/1e6,[],3);
pcolor(EE,DD,psi_plot);
shading interp;
hold on
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
hold off
caxis([0 2.5]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 11]);
colormap(gca,cmocean('amp',25));
% colormap(gca,hot(25));
% xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)')
text(-6.5,27.3,'Overturning streamfunction s.d.','FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);
cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
title(cbhandle,'Sv');
handle = title(gca,'Seasonal cycle removed','Fontsize',fontsize+2);
set(handle,'Position',get(handle,'Position')+[0 -.17 0]);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);


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











%%% Leading EOF
psi_plot = eof_maps(:,:,1)*std(pc(1,:));
subplot('Position',axpos(3,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
hold on
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
hold off
caxis([-2.5 2.5]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 11]);
colormap(gca,cmocean('balance',50));
xlabel('MOC coordinate, \eta');
% ylabel('Potential density (kg/m^3)')
text(-6.5,27.3,['First EOF (',num2str(expvar(1),'%.1f'),'% of variance)'],'FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);


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








%%% Leading EOF, no seasonal cycle
psi_plot = eof_maps_noseason(:,:,1)*std(pc_noseason(1,:));
subplot('Position',axpos(4,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
hold on
plot([4 4],ylim,'--','Color',[.3 .3 .3]);
plot([0 0],ylim,'--','Color',[.3 .3 .3]);
hold off
caxis([-2.5 2.5]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 11]);
colormap(gca,cmocean('balance',50));
xlabel('MOC coordinate, \eta');
% ylabel('Potential density (kg/m^3)')
text(-6.5,27.3,['First EOF (',num2str(expvar_noseason(1),'%.1f'),'% of variance)'],'FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);
text(-4.5,28.15,'FRIS','FontSize',fontsize);
text(1.3,28.15,'Shelf','FontSize',fontsize);
text(5.5,28.15,'Weddell Sea','FontSize',fontsize);


cbhandle = colorbar;
set(cbhandle,'Position',cb2_pos);
title(cbhandle,'Sv');

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







%%% Time series of streamfunction strength
subplot('Position',axpos(5,:));
plot(times/t1year+2007,psimax_full/1e6);
hold on;
plot(times/t1year+2007,0*times,'k--');
hold off;
xlabel('Time (years)');
ylabel('Overturning (Sv)')
set(gca,'FontSize',fontsize);
set(gca,'YColor',colororder(1,:));
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
plot([0 0],[0 0],'Color',colororder(1,:));
hold on;
plot(times/t1year+2007,pc(1,:)/std(pc(1,:)),'Color',colororder(2,:));
plot(times/t1year+2007,-(bflux_int)/std(bflux_int),'Color',colororder(3,:));
% plot(times/t1year,taux_ews_int/std(taux_ews_int),'Color',colororder(3,:));
hold off;
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Bottom');
set(ax2,'YAxisLocation','Right');
box off;
set(ax2,'XLim',get(ax1,'XLim'));
set(ax2,'FontSize',fontsize);
% ylabel(ax2,'PC1/normalized buoyancy loss');
% set(get(ax2,'YLabel'),'Rotation',270)
% set(get(ax2,'YLabel'),'Position',get(get(ax2,'YLabel'),'Position')+[0.4 0 0])
leghandle = legend('Overturning','PC1','Buoyancy loss','Location','NorthWest');
set(leghandle,'Orientation','Horizontal');






%%% Time series of streamfunction strength, interannual
subplot('Position',axpos(6,:));
plot(times/t1year+2007,smooth(psimax_full,12)/1e6)
xlabel('Time (years)');
% ylabel('Overturning (Sv)')
set(gca,'FontSize',fontsize);
set(gca,'YColor',colororder(1,:));
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
plot([0 0],[0 0],'Color',colororder(1,:));
hold on;
plot(times/t1year+2007,smooth(pc_noseason(1,:)/std(pc_noseason(1,:)),12),'Color',colororder(2,:));
% plot(times/t1year,-smooth((bflux_int_shelf-mean(bflux_int_shelf))/std(bflux_int_shelf),12),'Color',colororder(3,:));
% plot(times/t1year,-smooth((bflux_int_shelf_noseason-mean(bflux_int_shelf_noseason))/std(bflux_int_shelf_noseason),12),'Color',colororder(3,:));
plot(times/t1year+2007,-smooth((bflux_int_shelf_noseason)/std(bflux_int_shelf_noseason),12,'moving'),'Color',colororder(3,:));
% plot(times/t1year,-smooth((bflux_int_noseason-mean(bflux_int_noseason))/std(bflux_int_noseason),12),'Color',colororder(3,:));
% plot(times/t1year,-(taux_ews_int_noseason-mean(taux_ews_int_noseason))/std(taux_ews_int_noseason),'Color',colororder(3,:));
% plot(times/t1year,taux_ews_int/std(taux_ews_int),'Color',colororder(3,:));
plot(times/t1year+2007,0*times,'k--');
hold off;
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Bottom');
set(ax2,'YAxisLocation','Right');
box off;
set(ax2,'XLim',get(ax1,'XLim'));
set(ax2,'FontSize',fontsize);
ylabel(ax2,'PC1; normalized buoyancy loss');
set(get(ax2,'YLabel'),'Rotation',270)
set(get(ax2,'YLabel'),'Position',get(get(ax2,'YLabel'),'Position')+[0.4 0 0])


pvar = -smooth((bflux_int_shelf_noseason)/std(bflux_int_shelf_noseason),12,'moving');
% pvar = -((bflux_int_shelf_noseason)/std(bflux_int_shelf_noseason))';
evar = smooth(psimax_full,12)/1e6;
% evar = smooth(pc_noseason(1,:)/std(pc_noseason(1,:)),12);
[r,lags] = xcorr(pvar,evar,'unbiased');
maxlag = lags(find(r==max(r)))
for maxlag = 6:18
  maxlag
  corr(pvar(1:end-maxlag),evar(maxlag+1:end))
end








%%% Panel labels
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.07 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.07 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.07 axpos(4,2)-0.04 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(5,1)-0.05 axpos(5,2)-0.04 0.03 0.03],'String',axlabels{5},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(6,1)-0.05 axpos(6,2)-0.04 0.03 0.03],'String',axlabels{6},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');









%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE PLOT 2 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
                                
fontsize = 14;
framepos = [417    26   1030   640];
axpos = zeros(6,4);
axpos(1,:) = [0.08 0.52 0.38 0.35];
axpos(2,:) = [0.54 0.52 0.38 0.35];
axpos(3,:) = [0.08 0.08 0.38 0.35];
axpos(4,:) = [0.54 0.08 0.38 0.35];
cb1_pos = [0.94 0.52 0.015 0.35];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

figure(208);                
clf;
set(gcf,'Position',framepos);            
                









%%% Leading EOF, cavity only
psi_plot = eof_maps_cavity(:,:,1)*std(pc_cavity(1,:));
subplot('Position',axpos(1,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
caxis([-.5 .5]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 0]);
colormap(gca,cmocean('balance',25));
xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)')
text(-6.8,27.3,['First EOF (',num2str(expvar_cavity(1),'%.1f'),'% of variance)'],'FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);
handle = title(gca,'FRIS cavity (\eta < 0)','Fontsize',fontsize+2);
set(handle,'Position',get(handle,'Position')+[0 -.2 0]);

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








%%% Leading EOF, shelf only
psi_plot = -eof_maps_shelf(:,:,1)*std(pc_shelf(1,:));
subplot('Position',axpos(2,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
caxis([-.5 .5]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-3 4]);
colormap(gca,cmocean('balance',50));
xlabel('MOC coordinate, \eta');
% ylabel('Potential density (kg/m^3)')
text(-2.8,27.3,['First EOF (',num2str(expvar_shelf(1),'%.1f'),'% of variance)'],'FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);
handle = title(gca,'Continental shelf (\eta < 4)','Fontsize',fontsize+2);
set(handle,'Position',get(handle,'Position')+[0 -.2 0]);

cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
title(cbhandle,'Sv');

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







%%% Time series of streamfunction strength
subplot('Position',axpos(3,:));
plot(times/t1year+2007,psimax_cavity/1e6);
hold on;
plot(times/t1year+2007,0*times,'k--');
hold off;
xlabel('Year');
ylabel('Overturning (Sv)')
set(gca,'FontSize',fontsize);
set(gca,'YColor',colororder(1,:));
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
plot([0 0],[0 0],'Color',colororder(1,:));
hold on;
plot(times/t1year+2007,pc_cavity(1,:)/std(pc_cavity(1,:)),'Color',colororder(2,:));
plot(times/t1year+2007,-(bflux_int_shelf)/std(bflux_int_shelf),'Color',colororder(3,:));
hold off;
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Bottom');
set(ax2,'YAxisLocation','Right');
box off;
set(ax2,'XLim',get(ax1,'XLim'));
set(ax2,'FontSize',fontsize);
% ylabel(ax2,'PC1/normalized buoyancy loss');
% set(get(ax2,'YLabel'),'Rotation',270)
% set(get(ax2,'YLabel'),'Position',get(get(ax2,'YLabel'),'Position')+[0.4 0 0])
leghandle = legend('Overturning','PC1','Buoyancy loss','Location','SouthWest');
set(leghandle,'Orientation','Horizontal');










%%% Time series of streamfunction strength, interannual
subplot('Position',axpos(4,:));
plot(times/t1year+2007,smooth(psimax_shelf,12)/1e6)
xlabel('Year');
% ylabel('Overturning (Sv)')
set(gca,'FontSize',fontsize);
set(gca,'YColor',colororder(1,:));
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
plot([0 0],[0 0],'Color',colororder(1,:));
hold on;
plot(times/t1year+2007,-smooth(pc_shelf(1,:)/std(pc_shelf(1,:)),12),'Color',colororder(2,:));
% plot(times/t1year,-smooth((bflux_int_shelf-mean(bflux_int_shelf))/std(bflux_int_shelf),12),'Color',colororder(3,:));
% plot(times/t1year,-smooth((bflux_int_shelf_noseason-mean(bflux_int_shelf_noseason))/std(bflux_int_shelf_noseason),12),'Color',colororder(3,:));
% plot(times/t1year,-(bflux_int_shelf_noseason)/std(bflux_int_shelf_noseason),'Color',colororder(3,:));
plot(times/t1year+2007,-smooth((bflux_int_shelf_noseason)/std(bflux_int_shelf_noseason),12,'moving'),'Color',colororder(3,:));
% plot(times/t1year,-smooth((bflux_int_noseason-mean(bflux_int_noseason))/std(bflux_int_noseason),12),'Color',colororder(3,:));
% plot(times/t1year,-(taux_ews_int_noseason-mean(taux_ews_int_noseason))/std(taux_ews_int_noseason),'Color',colororder(3,:));
% plot(times/t1year,taux_ews_int/std(taux_ews_int),'Color',colororder(3,:));
plot(times/t1year+2007,0*times,'k--');
hold off;
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Bottom');
set(ax2,'YAxisLocation','Right');
box off;
set(ax2,'XLim',get(ax1,'XLim'));
set(ax2,'FontSize',fontsize);
ylabel(ax2,'PC1; normalized buoyancy loss');
set(get(ax2,'YLabel'),'Rotation',270)
set(get(ax2,'YLabel'),'Position',get(get(ax2,'YLabel'),'Position')+[0.4 0 0])




pvar = -smooth((bflux_int_shelf_noseason)/std(bflux_int_shelf_noseason),12,'moving');
evar = smooth(psimax_shelf,12)/1e6;
% evar = -smooth(pc_shelf(1,:)/std(pc_shelf(1,:)),12);
[r,lags] = xcorr(pvar,evar,'none');
maxlag = lags(find(r==max(r)))
for maxlag = 6:18
  maxlag
  corr(pvar(2:end-maxlag),evar(maxlag+2:end))
end








%%% Panel labels
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.07 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.05 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.05 axpos(4,2)-0.04 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');








