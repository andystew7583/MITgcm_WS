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
  = calcFRISMeltTimeSeries (expdir,expname,compute_eff_melt,compute_upstream,9,3); 

%%% Plot the time series
figure(30);
plot(tt/t1year,SHImelt/1e12*t1year);
if (compute_eff_melt)
  hold on;
  plot(tt/t1year,SHImelt_eff/1e12*t1year);
  hold off;
  legend('True melt','Effective melt');
end
xlabel('Time (years)');cd .
ylabel('FRIS melt rate (Gt/yr)');

figure(31);
plot(tt/t1year,SHImelt/1e12*t1year);
hold on;
plot(tt/t1year,SIprod/1e12*t1year);
plot(tt/t1year,(SHImelt+SIprod)/1e12*t1year);
plot(tt/t1year,smooth(SIprod,24)/1e12*t1year);
plot(tt/t1year,smooth(SHImelt+SIprod,24)/1e12*t1year);
hold off;
legend('FRIS melt','Sea ice production','Net freshwater flux');
xlabel('Time (years)');
ylabel('FRIS melt rate (Gt/yr)');
grid on;

figure(32);
pcolor(XC,YC,-SHImelt_mean/920*t1year);
shading interp;
colorbar;
caxis([-20 20]);
colormap redblue;

figure(33);
pcolor(XC,YC,-SIprod_mean/920*t1year);
shading interp;
colorbar;
caxis([-20 20]);
colormap redblue;
