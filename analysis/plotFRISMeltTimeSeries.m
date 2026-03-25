%%%
%%% plotFRISMeltTimeSeries.m
%%%
%%% Plots a time series of total FRIS melt rates for one of our
%%% experiments.
%%%

%%% Options
loadexp;

compute_eff_melt = false;

compute_upstream = false;

%%% Melt time series
[tt,SHImelt,SIprod,SHImelt_mean, ...
  SIprod_mean,XC,YC,bathy,SHELFICEtopo, ...
  SHImelt_eff,SHImelt_tend,SHImelt_diff, ...
  SHImelt_tflux,SHImelt_conv,FWupstream] ...
  = calcFRISMeltTimeSeries (expdir,expname,compute_eff_melt,compute_upstream); 

%%% Plot the time series
figure(30);
plot(tt/86400/365,SHImelt/1e12*86400*365);
if (compute_eff_melt)
  hold on;
  plot(tt/86400/365,SHImelt_eff/1e12*86400*365);;
  hold off;
  legend('True melt','Effective melt');
end
xlabel('Time (years)');cd .
ylabel('FRIS melt rate (Gt/yr)');

figure(31);
plot(tt/86400/365,SHImelt/1e12*86400*365);
hold on;
plot(tt/86400/365,SIprod/1e12*86400*365);
plot(tt/86400/365,(SHImelt+SIprod)/1e12*86400*365);
plot(tt/86400/365,smooth(SIprod,24)/1e12*86400*365);
hold off;
legend('FRIS melt','Sea ice production','Net freshwater flux');
xlabel('Time (years)');
ylabel('FRIS melt rate (Gt/yr)');
grid on;

figure(32);
pcolor(XC,YC,-SHImelt_mean/920*86400*365);
shading interp;
colorbar;
caxis([-20 20]);
colormap redblue;

figure(33);
pcolor(XC,YC,-SIprod_mean/920*86400*365);
shading interp;
colorbar;
caxis([-20 20]);
colormap redblue;
