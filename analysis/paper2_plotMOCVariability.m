% %%%
% %%% paper2_plotOverturning.m
% %%%
% %%% Plots the overturning circulation variability in density space.
% %%%
% 
% %%% Needed for EOF
% addpath CDT/cdt;
% 
% %%% Load experiment
% expdir = '../experiments';
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
% loadexp;
% 
% %%% Options (see calcOverturning)
% calc_psi_eddy = true;
% deform_cavity = false;
% use_layers = true;
% densvar = 'PD0';
% psimax = 6;
% psistep = 0.25;
% % psimax = 2;
% % psistep = 0.1;
% ylim = [27.3 28.2];
% % ylim = [27.3 29];
% 
% %%% Construct output file name
% outfname = [expname,'_MOC_',densvar];
% if (calc_psi_eddy)
%   if (use_layers)
%     estr = '_layers';
%   else
%     estr = '_TRM';
%   end
% else
%   estr = '_noeddy';
% end
% outfname = [outfname,estr];
% if (deform_cavity)
%   outfname = [outfname,'_deform'];
% end
% outfname = [outfname,'.mat'];
% load(fullfile('products',outfname));
% [DD,EE] = meshgrid(dens_levs,eta);
% 
% %%% Compute EOFs
% Nt = length(times);
% Neta = length(eta);
% Nd = length(dens_levs);
% psi_tot = psi_mean+psi_eddy;
% psi_tot_noseason = psi_tot - repmat(squeeze(mean(reshape(psi_tot,[Neta Nd 12,size(psi_tot,3)/12]),4)),[1 1 size(psi_tot,3)/12]);
% [eof_maps,pc,expvar] = eof((psi_tot)/1e6);
% [eof_maps_noseason,pc_noseason,expvar_noseason] = eof((psi_tot_noseason)/1e6);
[eof_maps_offshelf,pc_offshelf,expvar_offshelf] = eof((psi_tot)/1e6,'mask',EE>=4);
[eof_maps_shelf,pc_shelf,expvar_shelf] = eof((psi_tot)/1e6,'mask',EE<4);
[eof_maps_cavity,pc_cavity,expvar_cavity] = eof((psi_tot)/1e6,'mask',EE<0);
% 
% %%% Load pre-computed time-mean T and S
% outfname = [expname,'_TSfluxes'];
% outfname = [outfname,'.mat'];
% load(fullfile('products',outfname),'theta_tavg');
% load(fullfile('products',outfname),'salt_tavg');
% 
% %%% Load surface flux time series
% outfname = [expname,'_surfFluxes'];
% outfname = [outfname,'.mat'];
% load(fullfile('products',outfname));
% 
% %%% Mean "sea surface" salinity and temperature for flux calculation
% sst = NaN*ones(Nx,Ny);
% sss = NaN*ones(Nx,Ny);
% ssp = NaN*ones(Nx,Ny);
% for i=1:Nx
%   for j=1:Ny
%     kmin = find(hFacC(i,j,:)>0,1,'first');
%     if (~isempty(kmin))
%       sss(i,j) = salt_tavg(i,j,kmin);
%       sst(i,j) = theta_tavg(i,j,kmin);
%       ssp = -gravity*rhoConst*RC(kmin)/1e4;
%     end
%   end
% end
% [alpha,beta] = calcAlphaBeta(sss,sst,ssp);
% Cp = 4e3;
% 
% %%% Loop through outputs and compute buoyancy fluxes
% bflux_int = zeros(1,Nt);
% bflux_int_shelf = zeros(1,Nt);
% bflux_int_offshelf = zeros(1,Nt);
% for n=1:Nt
%   
%   %%% Composite surface heat flux
%   htFlxDn = tflux(:,:,n);
%   tmp = SHIhtFlx(:,:,n);
%   htFlxDn(hFacC(:,:,1)==0) = -tmp(hFacC(:,:,1)==0);
%   htFlxDn(sum(hFacC,3)==0) = NaN;
% 
%   %%% Composite virtual surface salt flux
%   sltFlxDn = sflux(:,:,n);
%   tmp = SHIfwFlx(:,:,n);
%   sltFlxDn(hFacC(:,:,1)==0) = tmp(hFacC(:,:,1)==0).*sss(hFacC(:,:,1)==0);
%   sltFlxDn(sum(hFacC,3)==0) = NaN;
% 
%   %%% Calculate buoyancy flux
%   bfluxT = gravity*(alpha.*htFlxDn/Cp/rhoConst);
%   bfluxS = - gravity*(beta.*sltFlxDn/rhoConst);
%   bflux = bfluxT + bfluxS;
%   
%   %%% Area-integrated buoyancy flux
%   bflux_int(n) = nansum(nansum(bflux.*RAC));   
%   bflux_int_shelf(n) = nansum(nansum(bflux.*RAC.*(ETA<4)));
%   bflux_int_offshelf(n) = nansum(nansum(bflux.*RAC.*(ETA>=4)));
% 
% end
% bflux_int_noseason = bflux_int - repmat(squeeze(mean(reshape(bflux_int,[1 12 length(bflux_int)/12]),3)),[1 length(bflux_int)/12]);
% bflux_int_shelf_noseason = bflux_int_shelf - repmat(squeeze(mean(reshape(bflux_int_shelf,[1 12 length(bflux_int)/12]),3)),[1 length(bflux_int)/12]);
% bflux_int_offshelf_noseason = bflux_int_offshelf - repmat(squeeze(mean(reshape(bflux_int_offshelf,[1 12 length(bflux_int)/12]),3)),[1 length(bflux_int)/12]);
% 
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
axpos(2,:) = [0.54 0.68 0.38 0.25];
axpos(3,:) = [0.08 0.36 0.38 0.25];
axpos(4,:) = [0.54 0.36 0.38 0.25];
axpos(5,:) = [0.08 0.05 0.38 0.25];
axpos(6,:) = [0.54 0.05 0.38 0.25];
cb1_pos = [0.94 0.68 0.015 0.25];
cb2_pos = [0.94 0.36 0.015 0.25];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

figure(207);                
clf;
set(gcf,'Position',framepos);            
                

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
% xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)')
text(-6.5,27.4,'Overturning streamfunction s.d., \psi','FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);
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
% xlabel('MOC coordinate, \eta');
ylabel('Potential density (kg/m^3)')
text(-6.5,27.4,'Overturning streamfunction s.d., \psi','FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);
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











%%% Leading EOF
psi_plot = -eof_maps(:,:,1)*std(pc(1,:));
subplot('Position',axpos(3,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
caxis([-2.5 2.5]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 11]);
colormap(gca,cmocean('balance',50));
xlabel('MOC coordinate, \eta');
% ylabel('Potential density (kg/m^3)')
text(-6.5,27.4,['First EOF (',num2str(expvar(1),'%.1f'),'% of variance)'],'FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);

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
psi_plot = -eof_maps_noseason(:,:,1)*std(pc_noseason(1,:));
subplot('Position',axpos(4,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
caxis([-2.5 2.5]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 11]);
colormap(gca,cmocean('balance',50));
xlabel('MOC coordinate, \eta');
% ylabel('Potential density (kg/m^3)')
text(-6.5,27.4,['First EOF (',num2str(expvar_noseason(1),'%.1f'),'% of variance)'],'FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);

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
plot(times/t1year,-pc(1,:)/std(pc(1,:)));
xlabel('Time (years)');
ylabel('(Sv)')
set(gca,'FontSize',fontsize);
% text(-6.5,27.4,'Overturning streamfunction, \psi','FontSize',fontsize+2);

% ax1 = gca;
% ax2 = axes('Position',get(ax1,'Position'));
% % plot(tt_onethird_notides/t1year+2007-18,0*tt_onethird_notides,'k-','LineWidth',0.5);
% % hold on;
% % plot(tt_onethird_notides/t1year+2007-18,0*tt_onethird_notides,'k--','LineWidth',0.5);
% % hold off;
% set(ax2,'Color','None');
% set(ax2,'XAxisLocation','Top');
% set(ax2,'YLim',get(ax1,'YLim'));
% set(ax2,'YTick',[]);
% box off;
% set(ax2,'XLim',get(ax1,'XLim')-78.16);
% set(ax2,'FontSize',fontsize);
% set(get(ax2,'XLabel'),'String','Pseudo-Latitude');







%%% Time series of streamfunction strength, shelf only
subplot('Position',axpos(6,:));
plot(times/t1year,-pc_noseason(1,:)/std(pc_noseason(1,:)));
xlabel('Time (years)');
ylabel('(Sv)')
set(gca,'FontSize',fontsize);









%%% Panel labels
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.07 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.07 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.07 axpos(4,2)-0.04 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(5,1)-0.07 axpos(5,2)-0.04 0.03 0.03],'String',axlabels{5},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(6,1)-0.07 axpos(6,2)-0.04 0.03 0.03],'String',axlabels{6},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');









%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE PLOT 2 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
                                
fontsize = 14;
framepos = [417    26   1030   600];
axpos = zeros(6,4);
axpos(1,:) = [0.08 0.55 0.38 0.38];
axpos(2,:) = [0.54 0.55 0.38 0.38];
axpos(3,:) = [0.08 0.05 0.38 0.38];
axpos(4,:) = [0.54 0.05 0.38 0.38];
cb1_pos = [0.94 0.55 0.015 0.38];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

figure(208);                
clf;
set(gcf,'Position',framepos);            
                









%%% Leading EOF, cavity only
psi_plot = -eof_maps_cavity(:,:,1)*std(pc_cavity(1,:));
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
text(-6.5,27.4,['First EOF (',num2str(expvar_cavity(1),'%.1f'),'% of variance)'],'FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);

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
text(-6.5,27.4,['First EOF (',num2str(expvar_shelf(1),'%.1f'),'% of variance)'],'FontSize',fontsize+2);
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);

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
plotyyy(times/t1year,psimax_cavity,times/t1year,-pc_cavity(1,:)/std(pc_cavity(1,:)),times/t1year,(mean(bflux_int_shelf)-bflux_int_shelf)/std(flux_int_shelf));
xlabel('Time (years)');
ylabel('(Sv)')
set(gca,'FontSize',fontsize);
% text(-6.5,27.4,'Overturning streamfunction, \psi','FontSize',fontsize+2);

% ax1 = gca;
% ax2 = axes('Position',get(ax1,'Position'));
% % plot(tt_onethird_notides/t1year+2007-18,0*tt_onethird_notides,'k-','LineWidth',0.5);
% % hold on;
% % plot(tt_onethird_notides/t1year+2007-18,0*tt_onethird_notides,'k--','LineWidth',0.5);
% % hold off;
% set(ax2,'Color','None');
% set(ax2,'XAxisLocation','Top');
% set(ax2,'YLim',get(ax1,'YLim'));
% set(ax2,'YTick',[]);
% box off;
% set(ax2,'XLim',get(ax1,'XLim')-78.16);
% set(ax2,'FontSize',fontsize);
% set(get(ax2,'XLabel'),'String','Pseudo-Latitude');







%%% Time series of streamfunction strength, shelf only
subplot('Position',axpos(4,:));
plotyy(times/t1year,psimax_shelf,times/t1year,-pc_shelf(1,:)/std(pc_shelf(1,:)));
xlabel('Time (years)');
ylabel('(Sv)')
set(gca,'FontSize',fontsize);









%%% Panel labels
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.07 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.07 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.07 axpos(4,2)-0.04 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');








