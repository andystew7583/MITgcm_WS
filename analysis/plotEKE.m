%%%
%%% plotEKE.m
%%%
%%% Plots EKE and EKE production.
%%%

%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% expname = 'hires_seq_onesixth_RTOPO2';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;

%%% Load pre-computed products
outfname = [expname,'_EKE.mat'];
load(fullfile('products',outfname));



%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%

figure(51);
pcolor(XC,YC,log10(EKE_zavg));
% pcolor(XC,YC,EKE_zavg);
shading interp;
colorbar;
colormap jet;

figure(52);
pcolor(XC,YC,PEtoEKE_zavg.*sum(hFacC.*DRF,3));
shading interp;
colorbar
colormap redblue;
% caxis([-2 2]*1e-4);

figure(53);
pcolor(XC,YC,MKEtoEKE_zavg.*sum(hFacC.*DRF,3));
shading interp;
colorbar
colormap redblue;
% caxis([-2 2]*1e-4);

figure(54)
pcolor(XC,YC,log10(max(EKE(:,:,:),[],3)));
shading interp
colorbar;
colormap jet;

