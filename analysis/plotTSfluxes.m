%%%
%%% plotTSfluxes.m
%%%
%%% Plots heat/salt fluxes in quasi-latitude coordinates
%%%

%%% Load experiment
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_RTOPO2';
expname = 'hires_seq_onetwelfth_notides_RTOPO2';


%%% Options (see calcTSfluxes)
deform_cavity = false;

%%% Load data file
outfname = [expname,'_TSfluxes'];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

%%% Make plots

figure(1);
plot(eta,thflux_stand);
hold on;
plot(eta,thflux_fluc);
plot(eta,mean(thflux_eddy,2));
plot(eta,mean(thflux_tot,2));
% plot(eta,tflux_eta/1027/4000);
hold off;
legend('Mean','Seasonal','Eddy','Total');


figure(2);
plot(eta,std(thflux_tot,[],2));
hold on;
plot(eta,std(thflux_mean,[],2));
plot(eta,std(thflux_eddy,[],2));
hold off;
legend('Total','Mean','Eddy');

figure(3);
plot(eta,sltflux_stand);
hold on;
plot(eta,sltflux_fluc);
plot(eta,mean(sltflux_eddy,2));
plot(eta,mean(sltflux_tot,2));
hold off;
legend('Mean','Seasonal','Eddy','Total');


figure(4);
plot(eta,std(sltflux_tot,[],2));
hold on;
plot(eta,std(sltflux_mean,[],2));
plot(eta,std(sltflux_eddy,[],2));
hold off;
legend('Total','Mean','Eddy');
