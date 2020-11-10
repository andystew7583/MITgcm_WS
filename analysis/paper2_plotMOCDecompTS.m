%%%
%%% paper2_plotMOCDecompTS.m
%%%
%%% Plots the m overturning circulation in density space.
%%%


%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%

% %%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_RTOPO2';
loadexp;
expname_notides = 'hires_seq_onetwelfth_notides_RTOPO2';

%%% Options (see calcOverturning)
calc_psi_eddy = true;
deform_cavity = false;
use_layers = true;
densvar = 'PD0';

%%% Load streamfunction
outfname = ['','_MOC_TS'];
if (calc_psi_eddy)  
  estr = '_TRM';  
else
  estr = '_noeddy';
end
outfname = [outfname,estr];
outfname = [outfname,'.mat'];
outfname_tides = [expname,outfname];
outfname_notides = [expname_notides,outfname];
NS = length(SS);
NT = length(TT);

%%% Grids
[TTT,SSS]=meshgrid(TT,SS);
DDD = densjmd95(SSS,TTT,-RC(1)*gravity*rhoConst/1e4*ones(size(SSS))) - 1000;



%%% Load MOC data file
load(fullfile('products',outfname_tides));
psi_mean_tides = -mean(psi_TS_mean_intT,3)/1e6;
psi_eddy_tides = -mean(psi_TS_eddy_intT,3)/1e6;
psi_mean_plot = -mean(psi_TS_mean_stand_intT,3)/1e6;
% psi_fluc_plot = -mean(psi_TS_mean_fluc_intT,3)/1e6;
psi_fluc_plot = -(mean(psi_TS_mean_intT,3)/1e6-mean(psi_TS_mean_stand_intT,3)/1e6);
psi_eddy_plot = -mean(psi_TS_eddy_intT,3)/1e6;
psi_tot_tides = -mean(psi_TS_mean_intT+psi_TS_eddy_intT,3)/1e6;
load(fullfile('products',outfname_notides));
psi_mean_notides = -mean(psi_TS_mean_intT,3)/1e6;
psi_eddy_notides = -mean(psi_TS_eddy_intT,3)/1e6;
psi_tot_notides = -mean(psi_TS_mean_intT+psi_TS_eddy_intT,3)/1e6;
psi_tide_plot = psi_tot_tides-psi_tot_notides;





%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.06 0.54 .4 .43];
axpos(2,:) = [0.53 0.54 .4 .43];
axpos(3,:) = [0.06 0.05 .4 .43];
axpos(4,:) = [0.53 0.05 .4 .43];
cbpos = [0.95 0.05 0.015 .93];
axlabels = {'(a)','(b)','(c)','(d)'};
psimax = 2;
psistep = 0.1;
colorcntrs = [-psimax:psistep:psimax];
linecntrs = [-1 -.5 -0.25 -.1 .1 0.25 0.5 1];
satlinecntrs = [-8 -6 -4 -3 -2]; 







%%% Set up the figure
figure(210)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[382          55        1120         930]);




%%% Mean streamfunction
subplot('Position',axpos(1,:));
pcolor(SSS,TTT,psi_mean_plot);
shading interp;
hold on;
[C,h] = contour(SSS,TTT,psi_mean_plot,linecntrs,'EdgeColor','k');
clabel(C,h,'Color','k');
[C,h] = contour(SSS,TTT,psi_mean_plot,satlinecntrs,'EdgeColor','w');
clabel(C,h,'Color','w');
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor',[0.5 0.5 0.5]);
clabel(C,h,'Color',[0.5 0.5 0.5]);
hold off;
colormap(gca,redblue(round(2*psimax/psistep)));
caxis([-psimax psimax]);
% xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',fontsize);
text(34.03,1,'Mean component','FontSize',fontsize+2,'fontweight','bold');
text(34.8,-1.8,'HSSW','FontSize',fontsize);
text(34.71,0.6,'WDW','FontSize',fontsize);
text(34.75,-2.3,'ISW','FontSize',fontsize);
text(34.4,-2,'WW','FontSize',fontsize);
text(34.1,-0.5,'AASW','FontSize',fontsize);

%%% Fluctuating streamfunction
subplot('Position',axpos(2,:));
pcolor(SSS,TTT,psi_fluc_plot);
shading interp;
hold on;
[C,h] = contour(SSS,TTT,psi_fluc_plot,linecntrs,'EdgeColor','k');
clabel(C,h,'Color','k');
[C,h] = contour(SSS,TTT,psi_fluc_plot,satlinecntrs,'EdgeColor','w');
clabel(C,h,'Color','w');
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor',[0.5 0.5 0.5]);
clabel(C,h,'Color',[0.5 0.5 0.5]);
hold off;
colormap(gca,redblue(round(2*psimax/psistep)));
caxis([-psimax psimax]);
% xlabel('S (g/kg)');
% ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',fontsize)
text(34.03,1,'Seasonal/interannual component','FontSize',fontsize+2,'fontweight','bold');
text(34.8,-1.8,'HSSW','FontSize',fontsize);
text(34.71,0.6,'WDW','FontSize',fontsize);
text(34.75,-2.3,'ISW','FontSize',fontsize);
text(34.4,-2,'WW','FontSize',fontsize);
text(34.1,-0.5,'AASW','FontSize',fontsize);

%%% Eddy streamfunction
subplot('Position',axpos(3,:));
pcolor(SSS,TTT,psi_eddy_plot);
shading interp;
hold on;
[C,h] = contour(SSS,TTT,psi_eddy_plot,linecntrs,'EdgeColor','k');
clabel(C,h,'Color','k');
[C,h] = contour(SSS,TTT,psi_eddy_plot,satlinecntrs,'EdgeColor','w');
clabel(C,h,'Color','w');
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor',[0.5 0.5 0.5]);
clabel(C,h,'Color',[0.5 0.5 0.5]);
hold off;
colormap(gca,redblue(round(2*psimax/psistep)));
caxis([-psimax psimax]);
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',fontsize);
text(34.03,1,'Eddy component','FontSize',fontsize+2,'fontweight','bold');
text(34.8,-1.8,'HSSW','FontSize',fontsize);
text(34.71,0.6,'WDW','FontSize',fontsize);
text(34.75,-2.3,'ISW','FontSize',fontsize);
text(34.4,-2,'WW','FontSize',fontsize);
text(34.1,-0.5,'AASW','FontSize',fontsize);

%%% Tidal streamfunction
subplot('Position',axpos(4,:));
pcolor(SSS,TTT,psi_tide_plot);
shading interp;
hold on;
[C,h] = contour(SSS,TTT,psi_tide_plot,linecntrs,'EdgeColor','k');
clabel(C,h,'Color','k');
[C,h] = contour(SSS,TTT,psi_tide_plot,satlinecntrs,'EdgeColor','w');
clabel(C,h,'Color','w');
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor',[0.5 0.5 0.5]);
clabel(C,h,'Color',[0.5 0.5 0.5]);
hold off;
colormap(gca,redblue(round(2*psimax/psistep)));
caxis([-psimax psimax]);
xlabel('S (g/kg)');
% ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',fontsize);
text(34.03,1,'Tide - No Tide','FontSize',fontsize+2,'fontweight','bold');
text(34.8,-1.8,'HSSW','FontSize',fontsize);
text(34.71,0.6,'WDW','FontSize',fontsize);
text(34.75,-2.3,'ISW','FontSize',fontsize);
text(34.4,-2,'WW','FontSize',fontsize);
text(34.1,-0.5,'AASW','FontSize',fontsize);

%%% Add colorbar
cbhandle = colorbar;
set(cbhandle,'Position',cbpos);
title(cbhandle,'Sv');



%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.04 axpos(4,2)-0.05 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
