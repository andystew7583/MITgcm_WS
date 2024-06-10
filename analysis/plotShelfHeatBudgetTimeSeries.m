%%%
%%% plotShelfHeatBudgetTimeSeries.m
%%%
%%% Plots heat fluxes and eddy-related quantities over the continental shelf.
%%%

addpath CDT/cdt;

%%% Load experiment
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_RTOPO2';
% expname = 'hires_seq_onetwelfth_RTOPO2';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to decompose eddy fluxes
calc_eddy_decomp = false;

%%% Define coordinate system for integrating to compute heatfunction
if (use_PsiBT)

  infname = [expname,'_TSfluxes'];
  load(fullfile('products',infname),'uvel_tavg');

  %%% Calculate depth-averaged zonal velocity
  UU = sum(uvel_tavg.*repmat(DRF,[Nx Ny 1]).*hFacW,3);
  clear('uvel_tavg');
  
  %%% Calculate barotropic streamfunction
  Psi = zeros(Nx+1,Ny+1);
  Psi(2:Nx+1,2:Ny+1) = -cumsum(UU.*DYG,2);
  Psi = Psi(1:Nx,1:Ny);
  
  %%% Interpolate to cell centers
  ETA = 0.25*(Psi(1:Nx,1:Ny)+Psi([2:Nx 1],1:Ny)+Psi(1:Nx,[2:Ny 1])+Psi([2:Nx 1],[2:Ny 1]))/1e6;
  
  %%% Streamunction grid for flux calculation
  eta = -2:.1:10;
  Neta = length(eta);

else

  ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
  eta = -9:.1:11;
  Neta = length(eta);

  %%% Bounds and horizontal mask for heat budget volume
  eta_icefront = -1.1;
  eta_shelfbreak = 3.5;  
  eidx_icefront = find(abs(eta-eta_icefront)==min(abs(eta-eta_icefront)));
  eidx_shelfbreak = find(abs(eta-eta_shelfbreak)==min(abs(eta-eta_shelfbreak)));
  zidx_icefront = 25;     
  msk = (ETA > eta_icefront) & (ETA < eta_shelfbreak);
  msk_slope = (ETA > eta_shelfbreak & ETA < 5); 
  msk_outershelf = (ETA > 2 & ETA < eta_shelfbreak); 

end

%%% Alternative geographical mask
% msk = (SHELFICEtopo>=0) & (bathy >= -1000) & (YC<YC(1,end-spongethickness+1)) & (XC<XC(end-spongethickness+1,1));

%%% Parameters
rho0 = 1000;
Cp = 4000;
ylim = [0 2000];
theta0 = -1.9;

%%% Load HeatFunction data file
outfname = [expname,'_HeatFunction'];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else
  if (deform_cavity)
    outfname = [outfname,'_deform'];
  end
end
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));

%%% Load shelf heat budget diagnostics
outfname = [expname,'_ShelfHeatBudget.mat'];
load(fullfile('products',outfname));

%%% TODO remove
% usq_eddy_int = usq_eddy_int(:,:,1:Ntime);
% vsq_eddy_int = vsq_eddy_int(:,:,1:Ntime);
% salt_int_lower = salt_int_lower(:,:,1:Ntime);
% salt_int_upper = salt_int_upper(:,:,1:Ntime);
% theta_int_lower = theta_int_lower(:,:,1:Ntime);
% theta_int_upper = theta_int_upper(:,:,1:Ntime);
% salt_bnd = salt_bnd(:,:,1:Ntime);
% theta_bnd = theta_bnd(:,:,1:Ntime);
% PEtoEKE_int = PEtoEKE_int(:,:,1:Ntime);
% MKEtoEKE_int = MKEtoEKE_int(:,:,1:Ntime);
% wt_tot_flux = wt_tot_flux(:,:,1:Ntime);
% ws_tot_flux = ws_tot_flux(:,:,1:Ntime);
% wt_mean_flux = wt_mean_flux(:,:,1:Ntime);
% ws_mean_flux = ws_mean_flux(:,:,1:Ntime);
% wt_eddy_flux = wt_eddy_flux(:,:,1:Ntime);
% ws_eddy_flux = ws_eddy_flux(:,:,1:Ntime);
% w_flux = w_flux(:,:,1:Ntime);
% tt = tt(1:Ntime);

%%% Compute volume averages
EKE_avg = zeros(1,length(tt));
EKE_slope = zeros(1,length(tt));
EKE_outershelf = zeros(1,length(tt));
PEtoEKE_avg = zeros(1,length(tt));
PEtoEKE_slope = zeros(1,length(tt));
PEtoEKE_outershelf = zeros(1,length(tt));
MKEtoEKE_avg = zeros(1,length(tt));
salt_upper_avg = zeros(1,length(tt));
salt_lower_avg = zeros(1,length(tt));
theta_upper_avg = zeros(1,length(tt));
theta_lower_avg = zeros(1,length(tt));
volW_lower = sum(sum(sum(RAW.*msk.*DRF(:,:,zidx_top:end).*hFacW(:,:,zidx_top:end))));
volS_lower = sum(sum(sum(RAS.*msk.*DRF(:,:,zidx_top:end).*hFacS(:,:,zidx_top:end))));
volC_lower = sum(sum(sum(RAC.*msk.*DRF(:,:,zidx_top:end).*hFacC(:,:,zidx_top:end))));
volW_upper = sum(sum(sum(RAW.*msk.*DRF(:,:,1:zidx_top-1).*hFacW(:,:,1:zidx_top-1))));
volS_upper = sum(sum(sum(RAS.*msk.*DRF(:,:,1:zidx_top-1).*hFacS(:,:,1:zidx_top-1))));
volC_upper = sum(sum(sum(RAC.*msk.*DRF(:,:,1:zidx_top-1).*hFacC(:,:,1:zidx_top-1))));
volC_slope = sum(sum(sum(RAC.*msk_slope.*DRF.*hFacC)));
volW_slope = sum(sum(sum(RAW.*msk_slope.*DRF.*hFacW)));
volS_slope = sum(sum(sum(RAS.*msk_slope.*DRF.*hFacS)));
volC_outershelf = sum(sum(sum(RAC.*msk_outershelf.*DRF.*hFacC)));
volW_outershelf = sum(sum(sum(RAW.*msk_outershelf.*DRF.*hFacW)));
volS_outershelf = sum(sum(sum(RAS.*msk_outershelf.*DRF.*hFacS)));
for n=1:length(tt)
  EKE_avg(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk)) / (volW_upper+volW_lower)  ...
             + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk)) / (volS_upper+volW_lower);
  EKE_slope(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk_slope)) / (volW_slope)  ...
             + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk_slope)) / (volS_slope);  
  EKE_outershelf(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk_outershelf)) / (volW_outershelf)  ...
             + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk_outershelf)) / (volS_outershelf);  
  PEtoEKE_avg(n) = sum(sum(PEtoEKE_int(:,:,n).*RAC.*msk)) / (volC_upper+volC_lower);
  PEtoEKE_slope(n) = sum(sum(PEtoEKE_int(:,:,n).*RAC.*msk_slope)) / (volC_slope);
  PEtoEKE_outershelf(n) = sum(sum(PEtoEKE_int(:,:,n).*RAC.*msk_outershelf)) / (volC_outershelf);
  MKEtoEKE_avg(n) = sum(sum(MKEtoEKE_int(:,:,n).*RAC.*msk)) / (volC_upper+volC_lower);             
  salt_lower_avg(n) = sum(sum(salt_int_lower(:,:,n).*RAC.*msk)) / volC_lower; 
  salt_upper_avg(n) = sum(sum(salt_int_upper(:,:,n).*RAC.*msk)) / volC_upper; 
  theta_lower_avg(n) = sum(sum(theta_int_lower(:,:,n).*RAC.*msk)) / volC_lower; 
  theta_pos_lower_avg(n) = sum(sum(theta_pos_int_lower(:,:,n).*RAC.*msk)) / volC_lower; 
  theta_neg_lower_avg(n) = sum(sum(theta_neg_int_lower(:,:,n).*RAC.*msk)) / volC_lower; 
  theta_upper_avg(n) = sum(sum(theta_int_upper(:,:,n).*RAC.*msk)) / volC_upper; 
end


wt_tot_mod = wt_tot_flux - w_flux.*theta0;
wt_mean_mod = wt_mean_flux - w_flux.*theta0;
hflux_vert = squeeze(sum(sum(wt_tot_mod.*RAC.*msk,1),2)*rho0*Cp/1e12);
hflux_mean_vert = squeeze(sum(sum(wt_mean_mod.*RAC.*msk,1),2)*rho0*Cp/1e12);
hflux_eddy_vert = squeeze(sum(sum(wt_eddy_flux.*RAC.*msk,1),2)*rho0*Cp/1e12);

psiT_mod = psiT_tot-psi_tot*theta0;

hflux_icefront = squeeze(psiT_mod(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_shelfbreak = squeeze(psiT_mod(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

hflux_eddy_icefront = squeeze(psiT_eddy(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_eddy_shelfbreak = squeeze(psiT_eddy(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

htend = diff(theta_lower_avg)'/(times(2)-times(1))*volC_lower/1e12*rho0*Cp;

figure(6);
plot(times(1:end-1)/t1year,hflux_icefront(1:end-1));
hold on;
plot(times(1:end-1)/t1year,-hflux_shelfbreak(1:end-1));
plot(times(1:end-1)/t1year,-hflux_vert(1:end-1));
plot(times(1:end-1)/t1year,-hflux_eddy_vert(1:end-1),'y--');
plot(times(1:end-1)/t1year,-hflux_mean_vert(1:end-1),'y:');
plot(times(1:end-1)/t1year,htend);
plot(times(1:end-1)/t1year,hflux_icefront(1:end-1)-hflux_shelfbreak(1:end-1)-hflux_vert(1:end-1) - htend,'k--');
plot(times(1:end-1)/t1year,0*times(1:end-1),'k:');
plot(times/t1year,hflux_eddy_icefront,'b--');
plot(times/t1year,-hflux_eddy_shelfbreak,'r--');
hold off

figure(7);
clf;
plot(times/t1year,theta_lower_avg);
hold on;
plot(times/t1year,theta_pos_lower_avg);
plot(times/t1year,theta_neg_lower_avg);
% plot(times/t1year,theta_upper_avg);
hold off;

figure(8);
clf;
plot(times(2:end)/t1year,diff(theta_lower_avg)/(times(2)-times(1))*volC_lower/1e12*rho0*Cp);
hold on;
plot(times(2:end)/t1year,hflux_icefront(2:end)-hflux_shelfbreak(2:end)-hflux_vert(2:end),'k--');
plot(times(2:end)/t1year,hflux_icefront(1:end-1)-hflux_shelfbreak(1:end-1)-hflux_vert(1:end-1),'k:');
hold off;


figure(11);
plot(times/t1year,salt_lower_avg);
hold on;
plot(times/t1year,salt_upper_avg);
hold off;

figure(12);
plot(times/t1year,PEtoEKE_avg);
hold on;
plot(times/t1year,MKEtoEKE_avg);
hold off

figure(13);
plot(times/t1year,EKE_avg);

figure(14);
plotyy(times/t1year,PEtoEKE_slope,times/t1year,hflux_eddy_vert);

figure(15);
plotyy(times/t1year,PEtoEKE_avg,times/t1year,hflux_eddy_vert);


figure(16);
plotyy(times/t1year,PEtoEKE_outershelf,times/t1year,hflux_eddy_vert);

figure(17);
plotyy(times/t1year,-EKE_outershelf.*theta_lower_avg,times/t1year,hflux_eddy_vert);

figure(18);
plotyy(times/t1year,EKE_outershelf,times/t1year,hflux_eddy_vert);

figure(100);scatter(PEtoEKE_outershelf,hflux_eddy_vert)

figure(101);scatter(EKE_outershelf.*(theta_lower_avg+1.9),hflux_eddy_vert)