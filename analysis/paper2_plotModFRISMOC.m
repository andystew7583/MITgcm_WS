%%%
%%% paper2_plotOverturning.m
%%%
%%% Plots the overturning circulation in density space.
%%%


%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_RTOPO2';
loadexp;

%%% Options (see calcOverturning)
calc_psi_eddy = true;
deform_cavity = true;
use_layers = true;
densvar = 'PD0';
psimax = 6;
psistep = 0.25;
% psimax = 2;
% psistep = 0.1;
ylim = [27.2 28.1];
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

%%% Load MOC data file
load(fullfile('products',outfname));
psi_plot = mean(psi_mean+psi_eddy,3)/1e6;
% psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);

%%% For plotting eta grid
bathy_plot = bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
ETA_plot = ETA;
ETA_plot(SHELFICEtopo-bathy<=0) = NaN;

%%% For plotting pot dens
idx = find(XC(:,1)<-38,1,'last');
dens_min = 27.5;
dens_max = 28.1;
dens_plot = dens_tavg;
dens_plot(hFacC==0)= NaN;
dens_plot = squeeze(dens_plot(idx,:,:));
dens_plot(dens_plot<dens_min) = dens_min;
dens_plot(dens_plot>dens_max) = dens_max;




%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE THE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
fontsize = 14;
framepos = [417    26   791   959];
labelspacing = 200;
axpos = zeros(6,4);
axpos(1,:) = [0.08 0.55 0.84 0.4];
axpos(2,:) = [0.08 0.05 0.84 0.4];
cb1_pos = [0.93 0.55 0.02 0.4];
cb2_pos = [0.93 0.05 0.02 0.4];
axlabels = {'(a)','(b)','(c)'};

figure(214);                
clf;
set(gcf,'Position',framepos);            
                

%%% Eta grid
subplot('Position',axpos(1,:));
pcolor(XC,YC,SHELFICEtopo-bathy_plot);
shading interp;
colormap(gca,flip(haxby,1))
hold on
% [C,h]=contour(XC,YC,ETA_plot,[-9:1:20],'EdgeColor','k');
[C,h]=contour(XC,YC,ETA_plot,[0.5 1 1.5 2:20],'EdgeColor','k');
plot(-38*ones(1,Ny),YC(1,:),'k--','Color',[.3 .3 .3],'LineWidth',2);
hold off;
clabel(C,h);
caxis([0 5000]);
ylabel('Latitude','interpreter','latex');
xlabel('Longitude','interpreter','latex');
text(-13,-82,'Modified MOC coordinate, $\tilde{\eta}$','FontSize',fontsize+2,'interpreter','latex');
set(gca,'FontSize',fontsize);
set(gca,'Color',[.85 .85 .85]);
cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
title(cbhandle,'m');


%%% Streamfunction
subplot('Position',axpos(2,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
hold on;
[C,h] = contour(EE,DD,psi_plot,[-6 -5 -4 -3 -2 -1 -0.5 -0.25 0.25 0.5 1],'EdgeColor','k');
clabel(C,h);
hold off;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[0 10.5]);
colormap(gca,cmocean('balance',round(2*psimax/(psistep))));%,'pivot',0));
xlabel('Modified MOC coordinate, $\tilde{\eta}$','interpreter','latex');
ylabel('Potential density (kg/m$^3$)','interpreter','latex')
set(gca,'FontSize',fontsize);
cbhandle = colorbar;
set(cbhandle,'Position',cb2_pos);
title(cbhandle,'Sv');
text(7,28.05,'Overturning streamfunction, $\psi^\sigma$','FontSize',fontsize+2,'interpreter','latex');


%%% Panel labels
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.07 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

