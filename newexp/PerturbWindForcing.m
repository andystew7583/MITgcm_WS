%%%
%%% PerturbWindForcing.m
%%%
%%% Modifies wind forcing files to increase winds at the front of the FRIS.
%%%

addpath ../analysis;

%%% Load base model parameters
defineGrid;

%%% "MOC grid" used to construct forcing structure that varies
%%% quasi-latitudinally
ETA = defineMOCgrid (EXF_XMC',EXF_YMC',[],[],false,false);

%%%%%%%%%% alpha
%%% choose max amount by which to perturb winds.
% alpha_max = 1.5;
% Rlat = 4; %%% Exponential decay scale 
% Rlon = 16;
% alpha_lon = -60;
% alpha_lat = -74.5;

%%% Parameters defining forcing perturbations
eta_trans_max = 4;
eta_trans_min = 3;
% alpha_u = 0.7; %%% 0.7 corresponds approximately to SSP585, see Neme et al. 2022
% alpha_v = 0.7;
% beta_u = 0.7;
% beta_v = 0.7;
alpha_u = 0.85; %%% Reduced by factor of 2 in an attempt to approximately halve buoyancy flux pert in combination with wind pert
alpha_v = 0.85;
beta_u = 0.85;
beta_v = 0.85;

%%% Matrices for applying forcing perturbations
%%% Linear variation of wind speed fraction from alpha_x for eta < eta_trans_min 
%%% to beta_x for eta > eta_trans_max, where x=u or x=v for zonal and
%%% meridional components
idx_trans = find((ETA >= eta_trans_min) & (ETA <= eta_trans_max)); %%% Transition indices
ufrac_mat = beta_u*ones(EXF_Nx,EXF_Ny);
ufrac_mat(ETA <= eta_trans_min) = alpha_u;
ufrac_mat(idx_trans) = alpha_u + (beta_u-alpha_u)*(ETA(idx_trans)-eta_trans_min)/(eta_trans_max-eta_trans_min);
vfrac_mat = beta_v*ones(EXF_Nx,EXF_Ny);
vfrac_mat(ETA <= eta_trans_min) = alpha_v;
vfrac_mat(idx_trans) = alpha_v + (beta_v-alpha_v)*(ETA(idx_trans)-eta_trans_min)/(eta_trans_max-eta_trans_min);


%%% Max componentwise wind speed
wind_max = 30;

%%% Define range of days to operate on
days_end = 3287;
days_start = 1;

%%% Load wind data
uwind = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
uwind_mod = 0*uwind;
fid = fopen(fullfile(inputfolder,zwind),'r','b');
for k=1:days_start:days_end
    uwind(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    uwind_mod(:,:,k-days_start+1) = uwind(:,:,k-days_start+1).*ufrac_mat;
end
fclose(fid);
vwind = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
vwind_mod = 0*vwind;
fid = fopen(fullfile(inputfolder,mwind),'r','b');
for k=days_start:days_end
    vwind(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    vwind_mod(:,:,k-days_start+1) = vwind(:,:,k-days_start+1).*vfrac_mat;
end
fclose(fid);

%%% Limit wind magnitudes
uwind_mod(uwind_mod>wind_max) = wind_max;
uwind_mod(uwind_mod<-wind_max) = -wind_max;
vwind_mod(vwind_mod>wind_max) = wind_max;
vwind_mod(vwind_mod<-wind_max) = -wind_max;

%%% Write modified wind forcing files
writeDataset(uwind,fullfile(inputfolder,'uwindfile_unmod.bin'),ieee,prec);
writeDataset(vwind,fullfile(inputfolder,'vwindfile_unmod.bin'),ieee,prec);
writeDataset(uwind_mod,fullfile(inputfolder,zwind),ieee,prec);
writeDataset(vwind_mod,fullfile(inputfolder,mwind),ieee,prec);