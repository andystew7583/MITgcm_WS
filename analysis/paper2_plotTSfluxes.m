%%%
%%% paper2_plotTSfluxes.m
%%%
%%% Plots heat/salt fluxes in quasi-latitude coordinates
%%%

% %%% Load experiment
% expdir = '../experiments';
% expname = 'hires_seq_onetwelfth_RTOPO2';
% expname_notides = 'hires_seq_onetwelfth_notides_RTOPO2';
% 
% %%% Options (see calcTSfluxes)
% deform_cavity = false;
% outfname = ['','_TSfluxes'];
% if (deform_cavity)
%   outfname = [outfname,'_deform'];
% end
% outfname = [outfname,'.mat'];
% outfname_tides = [expname,outfname];
% outfname_notides = [expname_notides,outfname];
% 
% %%% Load pre-computed fluxes
% load(fullfile('products',outfname_tides));
% thflux_tot_tides = thflux_stand+thflux_fluc+mean(thflux_eddy,2);
% thflux_tot_plot = thflux_tot_tides;
% sltflux_tot_tides = sltflux_stand+sltflux_fluc+mean(sltflux_eddy,2);
% sltflux_tot_plot = sltflux_tot_tides;
% load(fullfile('products',outfname_notides));
% thflux_mean_plot = thflux_stand;
% thflux_fluc_plot = thflux_fluc;
% thflux_eddy_plot = mean(thflux_eddy,2);
% thflux_tot_notides = thflux_stand+thflux_fluc+mean(thflux_eddy,2);
% thflux_tide_plot = thflux_tot_tides - thflux_tot_notides;
% sltflux_mean_plot = sltflux_stand;
% sltflux_fluc_plot = sltflux_fluc;
% sltflux_eddy_plot = mean(sltflux_eddy,2);
% sltflux_tot_notides = sltflux_stand+sltflux_fluc+mean(sltflux_eddy,2);
% sltflux_tide_plot = sltflux_tot_tides - sltflux_tot_notides;






%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(4,4);
axpos(1,:) = [0.07 0.54 .9 .43];
axpos(2,:) = [0.07 0.05 .9 .43];
cbpos = [0.95 0.05 0.015 .93];
axlabels = {'(a)','(b)'};
rho0 = 1027;
Cp = 4000;



%%% Set up the figure
figure(211)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[382    39   770   946]);

%%% Heat flux
subplot('Position',axpos(1,:));
plot(eta,rho0*Cp*thflux_mean_plot/1e12);
hold on;
plot(eta,rho0*Cp*thflux_fluc_plot/1e12);
plot(eta,rho0*Cp*thflux_eddy_plot/1e12);
plot(eta,rho0*Cp*thflux_tide_plot/1e12);
plot(eta,rho0*Cp*thflux_tot_plot/1e12);
plot(eta,0*eta,'k--');
hold off;
ylabel('Heat flux (TW)');
axis([-8 10 -15 5]);
set(gca,'FontSize',fontsize);
leghandle = legend('Mean','Seasonal/interannual','Eddy','Tide - No Tide','Total','Location','SouthWest');
set(leghandle,'FontSize',fontsize);

%%% Salt flux
subplot('Position',axpos(2,:));
plot(eta,rho0*sltflux_mean_plot/1e9);
hold on;
plot(eta,rho0*sltflux_fluc_plot/1e9);
plot(eta,rho0*sltflux_eddy_plot/1e9);
plot(eta,rho0*sltflux_tide_plot/1e9);
plot(eta,rho0*sltflux_tot_plot/1e9);
plot(eta,0*eta,'k--');
hold off;
ylabel('Salt flux (Gg/s)');
xlabel('MOC coordinate, \eta');
axis([-8 10 -.4 1.4]);
set(gca,'FontSize',fontsize);
