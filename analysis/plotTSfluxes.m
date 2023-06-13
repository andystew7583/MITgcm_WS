%%%
%%% plotTSfluxes.m
%%%
%%% Plots heat/salt fluxes in quasi-latitude coordinates
%%%

%%% Load experiment
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';

rho0 = 1000;
Cp = 4000;


%%% Options (see calcTSfluxes)
deform_cavity = false;

%%% Load data file
outfname = [expname,'_TSfluxes'];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

thflux_tot = thflux_stand+thflux_fluc+mean(thflux_eddy,2);
sltflux_tot = sltflux_stand+sltflux_fluc+mean(sltflux_eddy,2);

%%% Make plots

figure(2);
plot(eta,thflux_stand*rho0*Cp);
hold on;
plot(eta,thflux_fluc*rho0*Cp);
plot(eta,mean(thflux_eddy*rho0*Cp,2));
plot(eta,mean(thflux_tot*rho0*Cp,2));
% plot(eta,tflux_eta/1027/4000);
hold off;
legend('Mean','Seasonal','Eddy','Total');


% figure(2);
% plot(eta,std(thflux_tot,[],2));
% hold on;
% plot(eta,std(thflux_mean,[],2));
% plot(eta,std(thflux_eddy,[],2));
% hold off;pl

figure(3);
plot(eta,sltflux_stand*rho0);
hold on;
plot(eta,sltflux_fluc*rho0);
plot(eta,mean(sltflux_eddy*rho0,2));
plot(eta,mean(sltflux_tot*rho0,2));
hold off;
legend('Mean','Seasonal','Eddy','Total');


% figure(4);
% plot(eta,std(sltflux_tot,[],2));
% hold on;
% plot(eta,std(sltflux_mean,[],2));
% plot(eta,std(sltflux_eddy,[],2));
% hold off;
% legend('Total','Mean','Eddy');
