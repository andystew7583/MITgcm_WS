%%%
%%% plotFRISMeltTimeSeries.m
%%%
%%% Plots a time series of total FRIS melt rates for one of our
%%% experiments.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;

%%% Melt time series
[tt,SHImelt,SHImelt_mean] = calcFRISMeltTimeSeries (expdir,expname); 

%%% Plot the time series
figure(10);
plot(tt/86400/365,SHImelt/1e12*86400*365);
xlabel('Time (years)');
ylabel('FRIS melt rate (Gt/yr)');

figure(12);
pcolor(XC,YC,-SHImelt_mean/920*86400*365);
shading interp;
colorbar;
caxis([-5 5]);
colormap redblue;
