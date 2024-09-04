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
% loadexp;

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
% salt0 = 34.72;
salt0 = 34.6;
fontsize = 14;

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

%%% Some quantities computed separately for higher-resolution experiments
if (strcmp(expname,'hires_seq_onetwentyfourth_notides_RTOPO2'))

  startyear = 2008;

  %%% Load eddy-induced components of heatfunction/transport
  outfname = [expname,'_HeatFunctionEddyDecomp'];
  if (use_PsiBT)
    outfname = [outfname,'_PsiBT'];
  else
    if (deform_cavity)
      outfname = [outfname,'_deform'];
    end
  end
  outfname = [outfname,'.mat'];
  load(fullfile('products',outfname));
  
  %%% Load positive and negative heatfunction components
  outfname = [expname,'_PosNegHeatFunction'];
  if (use_PsiBT)
    outfname = [outfname,'_PsiBT'];
  else
    if (deform_cavity)
      outfname = [outfname,'_deform'];
    end
  end
  outfname = [outfname,'.mat'];
  load(fullfile('products',outfname));

else

  startyear = 2007;

end

%%% Load shelf heat budget diagnostics
outfname = [expname,'_ShelfHeatBudget.mat'];
load(fullfile('products',outfname));

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
theta_pos_lower_avg = zeros(1,length(tt));
theta_neg_lower_avg = zeros(1,length(tt));
volW_lower = sum(sum(sum(RAW.*msk.*DRF(:,:,zidx_icefront:end).*hFacW(:,:,zidx_icefront:end))));
volS_lower = sum(sum(sum(RAS.*msk.*DRF(:,:,zidx_icefront:end).*hFacS(:,:,zidx_icefront:end))));
volC_lower = sum(sum(sum(RAC.*msk.*DRF(:,:,zidx_icefront:end).*hFacC(:,:,zidx_icefront:end))));
volW_upper = sum(sum(sum(RAW.*msk.*DRF(:,:,1:zidx_icefront-1).*hFacW(:,:,1:zidx_icefront-1))));
volS_upper = sum(sum(sum(RAS.*msk.*DRF(:,:,1:zidx_icefront-1).*hFacS(:,:,1:zidx_icefront-1))));
volC_upper = sum(sum(sum(RAC.*msk.*DRF(:,:,1:zidx_icefront-1).*hFacC(:,:,1:zidx_icefront-1))));
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
hflux_eddy_adv_vert = squeeze(sum(sum(w_eddy_flux.*(theta_bnd-theta0).*RAC.*msk,1),2)*rho0*Cp/1e12);
hflux_eddy_stir_vert = hflux_eddy_vert - hflux_eddy_adv_vert;

psiT_mod = psiT_tot-psi_tot*theta0;
psiT_mean_mod = psiT_mean-psi_tot*theta0;

hflux_icefront = squeeze(psiT_mod(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_shelfbreak = squeeze(psiT_mod(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

hflux_mean_icefront = squeeze(psiT_mean_mod(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_mean_shelfbreak = squeeze(psiT_mean_mod(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

hflux_eddy_icefront = squeeze(psiT_eddy(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_eddy_shelfbreak = squeeze(psiT_eddy(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

hflux_pos_icefront = squeeze(psiT_pos_mean(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_pos_shelfbreak = squeeze(psiT_pos_mean(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

hflux_neg_icefront = squeeze(psiT_neg_mean(eidx_icefront,zidx_icefront,:))*rho0*Cp/1e12;
hflux_neg_shelfbreak = squeeze(psiT_neg_mean(eidx_shelfbreak,zidx_icefront,:))*rho0*Cp/1e12;

Nt = length(tt);
M = zeros(Nt,Nt);
M(1,1) = 3/4;
M(1,2) = 1/8;
M(1,Nt) = 1/8;
for n=2:Nt-1
  M(n,n-1) = 1/8;
  M(n,n) = 3/4;
  M(n,n+1) = 1/8;
end
M(Nt,Nt-1) = 1/8;
M(Nt,Nt) = 3/4;
M(Nt,1) = 1/8;
% for n=1:Nt-1
%   M(n,n) = 0.5;
%   M(n,n+1) = 0.5;
% end
% M(Nt,Nt) = 1;
% % M(Nt,1) = 0.5;
theta_lower_inst = M \ theta_lower_avg';

htend = (theta_lower_inst([2:end 1]) - theta_lower_inst([end 1:end-1]))/(2*(times(2)-times(1)))*volC_lower/1e12*rho0*Cp;
htend = htend(1:end-1);


% htend = (theta_lower_avg([2:end 1]) - theta_lower_avg([end 1:end-1]))/(2*(times(2)-times(1)))*volC_lower/1e12*rho0*Cp;
% htend = htend(1:end-1)';

% htend = diff(theta_lower_avg)'/(times(2)-times(1))*volC_lower/1e12*rho0*Cp;

hflux_int = -mean(psiT_mod(:,1,:),3) * rho0*Cp/1e12;
hflux_mean_int = -mean(psiT_mean_mod(:,1,:),3)* rho0*Cp/1e12;
hflux_eddy_int = -mean(psiT_eddy(:,1,:),3)* rho0*Cp/1e12;






ws_tot_mod = ws_tot_flux - w_flux.*salt0;
ws_mean_mod = ws_mean_flux - w_flux.*salt0;
sflux_vert = squeeze(sum(sum(ws_tot_mod.*RAC.*msk,1),2)*rho0/1e9);
sflux_mean_vert = squeeze(sum(sum(ws_mean_mod.*RAC.*msk,1),2)*rho0/1e9);
sflux_eddy_vert = squeeze(sum(sum(ws_eddy_flux.*RAC.*msk,1),2)*rho0/1e9);
sflux_eddy_adv_vert = squeeze(sum(sum(w_eddy_flux.*(salt_bnd-salt0).*RAC.*msk,1),2)*rho0/1e9);
sflux_eddy_stir_vert = sflux_eddy_vert - sflux_eddy_adv_vert;

psiS_mod = psiS_tot-psi_tot*salt0;
psiS_mean_mod = psiS_mean-psi_tot*salt0;

sflux_icefront = squeeze(psiS_mod(eidx_icefront,zidx_icefront,:))*rho0/1e9;
sflux_shelfbreak = squeeze(psiS_mod(eidx_shelfbreak,zidx_icefront,:))*rho0/1e9;

sflux_mean_icefront = squeeze(psiS_mean_mod(eidx_icefront,zidx_icefront,:))*rho0/1e9;
sflux_mean_shelfbreak = squeeze(psiS_mean_mod(eidx_shelfbreak,zidx_icefront,:))*rho0/1e9;

sflux_eddy_icefront = squeeze(psiS_eddy(eidx_icefront,zidx_icefront,:))*rho0/1e9;
sflux_eddy_shelfbreak = squeeze(psiS_eddy(eidx_shelfbreak,zidx_icefront,:))*rho0/1e9;

% sflux_pos_icefront = squeeze(psiS_pos_mean(eidx_icefront,zidx_icefront,:))*rho0/1e9;
% sflux_pos_shelfbreak = squeeze(psiS_pos_mean(eidx_shelfbreak,zidx_icefront,:))*rho0/1e9;
% 
% sflux_neg_icefront = squeeze(psiS_neg_mean(eidx_icefront,zidx_icefront,:))*rho0/1e9;
% sflux_neg_shelfbreak = squeeze(psiS_neg_mean(eidx_shelfbreak,zidx_icefront,:))*rho0/1e9;

% stend = diff(salt_lower_avg)'/(times(2)-times(1))*volC_lower/1e9*rho0;
stend = (salt_lower_avg([2:end 1]) - salt_lower_avg([end 1:end-1]))/(2*(times(2)-times(1)))*volC_lower/1e9*rho0;
stend = stend(1:end-1)';






figure(5);
plot(eta,hflux_int,'LineWidth',1.5);
hold on;
plot(eta,hflux_mean_int,'LineWidth',1.5);
plot(eta,hflux_eddy_int,'LineWidth',1.5);
hold off;
set(gca,'XLim',[-9 4]);
xlabel('MOC coordinate \eta');
ylabel('Heat flux (TW)');
set(gca,'FontSize',14);
legend('Total','Mean','Eddy');

thecolororder = get(gca,'ColorOrder');

figure(6);
plot(startyear+times(1:end-1)/t1year,hflux_icefront(1:end-1),'LineWidth',1.5);
hold on;
plot(startyear+times(1:end-1)/t1year,-hflux_shelfbreak(1:end-1),'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,-hflux_vert(1:end-1),'LineWidth',1.5);
% plot(times(1:end-1)/t1year,-hflux_mean_vert(1:end-1),'y:','LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,htend,'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,hflux_icefront(1:end-1)-hflux_shelfbreak(1:end-1)-hflux_vert(1:end-1) - htend,'k--','LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,0*times(1:end-1),'k:','LineWidth',1.5);
hold off
legend('Q_c_a_v_i_t_y','Q_s_h_e_l_f','Q_v_e_r_t','Tendency','Sum');
xlabel('Year');
ylabel('Heat flux (TW)');
set(gca,'Position',[0.05 0.1 0.9 0.85]);
set(gca,'FontSize',14);


figure(7);
plot(startyear+times(1:end-1)/t1year,hflux_icefront(1:end-1),'LineWidth',1.5);
hold on;
plot(startyear+times(1:end-1)/t1year,-hflux_shelfbreak(1:end-1),'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,-hflux_vert(1:end-1),'LineWidth',1.5);
% plot(times(1:end-1)/t1year,-hflux_mean_vert(1:end-1),'y:','LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,htend,'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,hflux_icefront(1:end-1)-hflux_shelfbreak(1:end-1)-hflux_vert(1:end-1) - htend,'k--','LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,0*times(1:end-1),'k:','LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,hflux_eddy_icefront(1:end-1),'--','LineWidth',1.5,'Color',thecolororder(1,:));
plot(startyear+times(1:end-1)/t1year,-hflux_eddy_shelfbreak(1:end-1),'--','LineWidth',1.5,'Color',thecolororder(2,:));
plot(startyear+times(1:end-1)/t1year,-hflux_eddy_vert(1:end-1),'--','LineWidth',1.5,'Color',thecolororder(3,:));
hold off
legend('Q_c_a_v_i_t_y','Q_s_h_e_l_f','Q_v_e_r_t','Tendency','Sum','Q_c_a_v_i_t_y (eddy)','Q_s_h_e_l_f (eddy)','Q_v_e_r_t (eddy)');
xlabel('Year');
ylabel('Heat flux (TW)');
set(gca,'Position',[0.05 0.1 0.9 0.85]);
set(gca,'FontSize',14);

figure(9);
plot(startyear+times(1:end-1)/t1year,hflux_eddy_vert(1:end-1));
hold on;
plot(startyear+times(1:end-1)/t1year,hflux_eddy_adv_vert(1:end-1));
plot(startyear+times(1:end-1)/t1year,hflux_eddy_stir_vert(1:end-1));
hold off;
title('Vertical eddy heat flux (TW)');
set(gca,'Position',[0.1 0.15 0.8 0.75]);
legend('Total','Advection','Stirring');
set(gca,'FontSize',14);
xlabel('Year');

% figure(7);
% clf;
% plot(times/t1year,theta_lower_avg-theta0);
% hold on;
% plot(times/t1year,theta_pos_lower_avg);
% plot(times/t1year,theta_neg_lower_avg);
% % plot(times/t1year,theta_upper_avg);
% hold off;
% 
% figure(8);
% clf;
% plot(times(2:end)/t1year,diff(theta_lower_avg)/(times(2)-times(1))*volC_lower/1e12*rho0*Cp);
% hold on;
% plot(times(2:end)/t1year,hflux_icefront(2:end)-hflux_shelfbreak(2:end)-hflux_vert(2:end),'k--');
% plot(times(2:end)/t1year,hflux_icefront(1:end-1)-hflux_shelfbreak(1:end-1)-hflux_vert(1:end-1),'k:');
% hold off;

figure(10);
plot(startyear+times(1:end-1)/t1year,hflux_icefront(1:end-1));
hold on;
plot(startyear+times(1:end-1)/t1year,hflux_mean_icefront(1:end-1));
plot(startyear+times(1:end-1)/t1year,hflux_eddy_icefront(1:end-1));
plot(startyear+times(1:end-1)/t1year,hflux_pos_icefront(1:end-1));
plot(startyear+times(1:end-1)/t1year,hflux_neg_icefront(1:end-1));
hold off;
title('Ice front heat flux (TW)');
set(gca,'Position',[0.1 0.15 0.8 0.75]);
legend('Total','Mean','Eddy','Pos','Neg');
set(gca,'FontSize',14);
xlabel('Year');

figure(11);
plot(times(1:end-1)/t1year,hflux_shelfbreak(1:end-1));
hold on;
plot(times(1:end-1)/t1year,hflux_mean_shelfbreak(1:end-1));
plot(times(1:end-1)/t1year,hflux_eddy_shelfbreak(1:end-1));
plot(times(1:end-1)/t1year,hflux_pos_shelfbreak(1:end-1));
plot(times(1:end-1)/t1year,hflux_neg_shelfbreak(1:end-1));
plot(times(1:end-1)/t1year,hflux_neg_shelfbreak(1:end-1)+hflux_pos_shelfbreak(1:end-1),'k--');
hold off;
title('Shelf break heat flux');
legend('Total','Mean','Eddy','Pos','Neg','Pos+Neg');

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
ax = plotyy(startyear+times/t1year,PEtoEKE_outershelf*rho0*Cp,startyear+times/t1year,hflux_eddy_vert);
xlabel('Year');
ylabel(ax(1),'Baroclinic EKE production (W/m^3)');
ylabel(ax(2),'Vertical eddy heat flux (TW)');
set(gca,'Position',[0.1 0.1 0.8 0.85]);
set(gca,'FontSize',fontsize);
set(ax(2),'FontSize',fontsize);

figure(17);
plotyy(times/t1year,-EKE_outershelf.*theta_lower_avg,times/t1year,hflux_eddy_vert);

figure(18);
plotyy(times/t1year,EKE_outershelf,times/t1year,hflux_eddy_vert);

figure(19);
plot(times/t1year,salt_lower_avg);
hold on;
plot(times/t1year,salt_upper_avg);
hold off;

figure(100);scatter(PEtoEKE_outershelf,hflux_eddy_vert)

figure(101);scatter(EKE_outershelf,hflux_eddy_vert)

figure(46);
plot(startyear+times(1:end-1)/t1year,sflux_icefront(1:end-1),'LineWidth',1.5);
hold on;
plot(startyear+times(1:end-1)/t1year,-sflux_shelfbreak(1:end-1),'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,-sflux_vert(1:end-1),'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,stend,'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,sflux_icefront(1:end-1)-sflux_shelfbreak(1:end-1)-sflux_vert(1:end-1) - stend,'k--','LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,0*times(1:end-1),'k:','LineWidth',1.5);
hold off
legend('Q_c_a_v_i_t_y','Q_s_h_e_l_f','Q_v_e_r_t','Tendency','Sum');
xlabel('Year');
ylabel('Salt flux (Gg/s)');
% set(gca,'Position',[0.05 0.1 0.9 0.85]);
set(gca,'Position',[0.1 0.15 0.8 0.75]);
set(gca,'FontSize',14);
print('./Figures/paper3/SaltBudget_onethird_notides_noeddy.png','-dpng');


figure(47);
plot(startyear+times(1:end-1)/t1year,sflux_icefront(1:end-1),'LineWidth',1.5);
hold on;
plot(startyear+times(1:end-1)/t1year,-sflux_shelfbreak(1:end-1),'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,-sflux_vert(1:end-1),'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,stend,'LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,sflux_icefront(1:end-1)-sflux_shelfbreak(1:end-1)-sflux_vert(1:end-1) - stend,'k--','LineWidth',1.5);
plot(startyear+times(1:end-1)/t1year,sflux_eddy_icefront(1:end-1),'--','LineWidth',1.5,'Color',thecolororder(1,:));
plot(startyear+times(1:end-1)/t1year,-sflux_eddy_shelfbreak(1:end-1),'--','LineWidth',1.5,'Color',thecolororder(2,:));
plot(startyear+times(1:end-1)/t1year,-sflux_eddy_vert(1:end-1),'--','LineWidth',1.5,'Color',thecolororder(3,:));
plot(startyear+times(1:end-1)/t1year,0*times(1:end-1),'k:','LineWidth',1.5);
hold off
legend('Q_c_a_v_i_t_y','Q_s_h_e_l_f','Q_v_e_r_t','Tendency','Sum','Q_c_a_v_i_t_y (eddy)','Q_s_h_e_l_f (eddy)','Q_v_e_r_t (eddy)');
xlabel('Year');
ylabel('Salt flux (Gg/s)');
% set(gca,'Position',[0.05 0.1 0.9 0.85]);
set(gca,'Position',[0.1 0.15 0.8 0.75]);
set(gca,'FontSize',14);
print('./Figures/paper3/SaltBudget_onethird_notides_eddy.png','-dpng');

figure(48);
plot(startyear+times(1:end-1)/t1year,sflux_eddy_vert(1:end-1));
hold on;
plot(startyear+times(1:end-1)/t1year,sflux_eddy_adv_vert(1:end-1));
plot(startyear+times(1:end-1)/t1year,sflux_eddy_stir_vert(1:end-1));
hold off;
title('Vertical eddy salt flux (Gg/s)');
set(gca,'Position',[0.1 0.15 0.8 0.75]);
legend('Total','Advection','Stirring');
set(gca,'FontSize',14);
xlabel('Year');
print('./Figures/paper3/SB_onethird_notides_vert.png','-dpng');


salt0 = 34.5;
figure(60);
pcolor(XC,YC,mean(w_eddy_flux.*(salt_bnd-salt0),3)*rho0*Cp);
shading interp;colorbar; 
colormap redblue; 
caxis([-100 100]);

figure(61);
pcolor(XC,YC,mean(ws_eddy_flux,3)*rho0*Cp);
shading interp;colorbar; 
colormap redblue; 
caxis([-100 100]);


figure(62);
pcolor(XC,YC,mean(ws_eddy_flux-w_eddy_flux.*(salt_bnd-salt0),3)*rho0*Cp);
shading interp;colorbar; 
colormap redblue; 
caxis([-100 100]);

figure(70);
pcolor(XC,YC,sum(theta_tavg.*DRF.*hFacC,3)./sum(DRF.*hFacC,3));
shading interp;
colorbar;
colormap jet(20);


figure(71);
pcolor(XC,YC,sum(salt_tavg.*DRF.*hFacC,3)./sum(DRF.*hFacC,3));
shading interp;
colorbar;
colormap jet(20);
