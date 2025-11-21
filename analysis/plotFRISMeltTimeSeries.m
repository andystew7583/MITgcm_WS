%%%
%%% plotFRISMeltTimeSeries.m
%%%
%%% Plots a time series of total FRIS melt rates for one of our
%%% experiments.
%%%

%%% Options
loadexp;

compute_eff_melt = false;

%%% Melt time series
[tt,SHImelt,SIprod,SHImelt_mean,XC,YC,bathy,SHELFICEtopo,SHImelt_eff,SHImelt_tend,SHImelt_diff,SHImelt_tflux] = calcFRISMeltTimeSeries (expdir,expname,compute_eff_melt); 

%%% Plot the time series
figure(10);
plot(tt/86400/365,SHImelt/1e12*86400*365);
if (compute_eff_melt)
  hold on;
  plot(tt/86400/365,SHImelt_eff/1e12*86400*365);;
  hold off;
  legend('True melt','Effective melt');
end
xlabel('Time (years)');cd .
ylabel('FRIS melt rate (Gt/yr)');

figure(11);
plot(tt/86400/365,SHImelt/1e12*86400*365);
hold on;
plot(tt/86400/365,SIprod/1e12*86400*365);
plot(tt/86400/365,(SHImelt+SIprod)/1e12*86400*365);
hold off;
legend('FRIS melt','Sea ice production','Net freshwater flux');
xlabel('Time (years)');cd .
ylabel('FRIS melt rate (Gt/yr)');

figure(12);
pcolor(XC,YC,-SHImelt_mean/920*86400*365);
shading interp;
colorbar;
caxis([-5 5]);
colormap redblue;
