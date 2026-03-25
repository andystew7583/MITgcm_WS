%%%
%%% compareFRISMeltTimeSeries.m
%%%
%%% Compares FRIS melt rates between experiments.
%%%

%%% Load experiment
expdir = '../experiments';
expname_tides = 'hires_seq_onethird_RTOPO2';
expname_notides = 'hires_seq_onethird_notides_RTOPO2';
expname = expname_tides;
loadexp; %%% To get grids

%%% Physical parameters
rhoi = 920;

%%% Load melt time series
[tt_tides,SHImelt_tides,SHImelt_mean_tides] = calcFRISMeltTimeSeries (expdir,expname_tides); 
[tt_notides,SHImelt_notides,SHImelt_mean_notides] = calcFRISMeltTimeSeries (expdir,expname_notides); 
SHImelt_tides = -SHImelt_tides/1e12*86400*365;
SHImelt_notides = -SHImelt_notides/1e12*86400*365;
tt_tides = tt_tides/86400/365;
tt_notides = tt_notides/86400/365;
SHImelt_mean_tides = -SHImelt_mean_tides/rhoi*86400*365;
SHImelt_mean_notides = -SHImelt_mean_notides/rhoi*86400*365;

%%% Plotting setup
colororder = get(gca,'colororder');
idx_mean = 213:324;
fontsize = 14;
framepos = [311         422        1099         451];
plotloc = [0.0855    0.1419    0.8854    0.7783];

%%% Make the plot
figure(20);
set(gcf,'Position',framepos);
plot(tt_tides,SHImelt_tides,'Color',colororder(1,:));
hold on;
plot(tt_notides,SHImelt_notides,'Color',colororder(2,:));
plot(tt_tides(idx_mean),mean(SHImelt_tides(idx_mean))*ones(size(idx_mean)),'--','Color',colororder(1,:));
plot(tt_notides(idx_mean),mean(SHImelt_notides(idx_mean))*ones(size(idx_mean)),'--','Color',colororder(2,:));
hold off;
xlabel('Time (years)');
ylabel('Melt (Gt/yr)');
title('FRIS melt rate');
legend('Tides','No tides','Location','SouthEast');
set(gca,'FontSize',fontsize);
set(gca,'Position',plotloc);
set(gca,'XLim',[min(tt_tides) max(tt_tides)]);
set(gca,'YLim',[0 250]);


%%% Plotting setup
colororder = get(gca,'colororder');
fontsize = 14;
framepos = [311         422        600         451];
plotloc = [0.1255    0.1419    0.7254    0.7783];

%%% Melt rate with tides
figure(21);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,SHImelt_mean_tides);
shading interp;
colorbar;
caxis([-5 5]);
axis([-83 -35 -83 -75]);
colormap(cmocean('balance',50));set(gca,'FontSize',fontsize);
set(gca,'Position',plotloc);
xlabel('Longitude');
ylabel('Latitude');
title('Melt rate with tides (m/yr)');

%%% Melt rate diffeeeeeerence
figure(22);
clf;
set(gcf,'Position',framepos);
pcolor(XC,YC,SHImelt_mean_tides-SHImelt_mean_notides);
shading interp;
colorbar;
caxis([-5 5]);
axis([-83 -35 -83 -75]);
colormap(cmocean('balance',50));set(gca,'FontSize',fontsize);
set(gca,'Position',plotloc);
xlabel('Longitude');
ylabel('Latitude');
title('Melt rate, tides-noTides (m/yr)');
