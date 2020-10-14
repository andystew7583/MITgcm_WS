%%%
%%% plotTSstreamfunction.m
%%%
%%% Plots overturning circulation in T/S space.
%%%

addpath CDT/cdt;

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onesixth_RTOPO2';
loadexp;

%%% Set true to compute eddy-induced transports
calc_psi_eddy = true; 

%%% Construct output file name
outfname = [expname,'_MOC_TS'];
if (calc_psi_eddy)  
  estr = '_TRM';  
else
  estr = '_noeddy';
end
outfname = [outfname,estr];
outfname = [outfname,'.mat'];

%%% Load streamfunction data
load(fullfile('products',outfname));
NS = length(SS);
NT = length(TT);

%%% Grids
[TTT,SSS]=meshgrid(TT,SS);
DDD = densjmd95(SSS,TTT,-RC(1)*gravity*rhoConst/1e4*ones(size(SSS))) - 1000;

%%% Plotting options
% psimax = 2;
% psistep = 0.1;
psimax = 8;
psistep = 0.25;
colorcntrs = [-psimax:psistep:psimax];
% linecntrs = [-psimax:psistep:-psistep psistep:psistep:psimax];
linecntrs = [];

%%% Make a plot
psi_plot = mean(-psi_TS_mean_intT,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(1);
clf;
set(gcf,'Position',[469   365   918   617]);
contourf(SSS,TTT,psi_plot,colorcntrs,'EdgeColor','None');
hold on;
contour(SSS,TTT,psi_plot,linecntrs,'EdgeColor','k');
hold off;
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
colormap(redblue(round(2*psimax/psistep)));
caxis([-psimax psimax]);
colorbar;
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',14);
set(gca,'Position',[0.1035    0.1100    0.7680    0.8430]);

%%% Make a plot
psi_plot = mean(-psi_TS_mean_intT-psi_TS_eddy_intT,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(2);
clf;
set(gcf,'Position',[469   365   918   617]);
contourf(SSS,TTT,psi_plot,colorcntrs,'EdgeColor','None');
hold on;
contour(SSS,TTT,psi_plot,linecntrs,'EdgeColor','k');
hold off;
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
colormap(redblue(round(2*psimax/psistep)));
caxis([-psimax psimax]);
colorbar;
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',14);
set(gca,'Position',[0.1035    0.1100    0.7680    0.8430]);

%%% Make a plot
psi_plot = mean(-psi_TS_eddy_intT,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(3);
clf;
set(gcf,'Position',[469   365   918   617]);
contourf(SSS,TTT,psi_plot,colorcntrs,'EdgeColor','None');
hold on;
contour(SSS,TTT,psi_plot,linecntrs,'EdgeColor','k');
hold off;
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
colormap(redblue(round(2*psimax/psistep)));
caxis([-psimax psimax]);
colorbar;
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',14);
set(gca,'Position',[0.1035    0.1100    0.7680    0.8430]);

%%% Make a plot
psi_plot = std(-psi_TS_mean_intT-psi_TS_eddy_intT,[],3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(4);
clf;
set(gcf,'Position',[469   365   918   617]);
pcolor(SSS,TTT,psi_plot);
shading interp;
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
colormap(hot);
caxis([0 2]);
colorbar;
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',14);
set(gca,'Position',[0.1035    0.1100    0.7680    0.8430]);

psi_TS_noseason = psi_TS_mean_intT + psi_TS_eddy_intT;
psi_TS_season = mean(reshape(psi_TS_noseason,[NS,NT,12,length(times)/12]),4);
psi_TS_noseason = psi_TS_noseason - repmat(psi_TS_season,[1 1 length(times)/12]);

%%% Make a plot
% psi_plot = mean(-psi_TS_season,3)/1e6;
psi_plot = -psi_TS_noseason(:,:,1)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(5);
clf;
set(gcf,'Position',[469   365   918   617]);
contourf(SSS,TTT,psi_plot,colorcntrs,'EdgeColor','None');
hold on;
contour(SSS,TTT,psi_plot,linecntrs,'EdgeColor','k');
hold off;
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
colormap(redblue(round(2*psimax/psistep)));
caxis([-psimax psimax]);
colorbar;
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',14);
set(gca,'Position',[0.1035    0.1100    0.7680    0.8430]);

%%% Make a plot
psi_plot = std(-psi_TS_noseason,[],3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(6);
clf;
set(gcf,'Position',[469   365   918   617]);
pcolor(SSS,TTT,psi_plot);
shading interp;
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
colormap(hot);
caxis([0 2]);
colorbar;
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',14);
set(gca,'Position',[0.1035    0.1100    0.7680    0.8430]);


%%% Make a plot
psi_plot = std(-psi_TS_season,[],3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
figure(7);
clf;
set(gcf,'Position',[469   365   918   617]);
pcolor(SSS,TTT,psi_plot);
shading interp;
hold on;
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
hold off;
clabel(C,h);
colormap(hot);
caxis([0 2]);
colorbar;
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',14);
set(gca,'Position',[0.1035    0.1100    0.7680    0.8430]);



Nt = length(times);
% [eof_maps,pc,expvar] = eof((-psi_TS_mean_intT-psi_TS_eddy_intT)/1e6);
[eof_maps,pc,expvar] = eof((-psi_TS_noseason)/1e6);


psi_plot = eof_maps(:,:,end);
figure(7);
clf;
set(gcf,'Position',[469   365   918   617]);
pcolor(SSS,TTT,psi_plot);
shading interp;
% hold on;
% [C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor','k');
% hold off;
% clabel(C,h);
colormap(redblue);
caxis([-.1 .1]);
colorbar;
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',14);
set(gca,'Position',[0.1035    0.1100    0.7680    0.8430]);

figure(9);
plot(times/t1year,pc(end,:));

