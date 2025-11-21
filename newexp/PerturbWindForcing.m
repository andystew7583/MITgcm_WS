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
alpha_eta_max = 4;
alpha_eta_min = 3;
alpha_u = 1;
alpha_v = 0.5;

%%% Matrices for applying forcing perturbations: linear variation from
%%% modified winds for eta < alpha_eta_min to unmodified winds for eta >
%%% alpha_eta_max
idx_trans = find((ETA >= alpha_eta_min) & (ETA <= alpha_eta_max)); %%% Transition indices
alpha_u_mat = ones(EXF_Nx,EXF_Ny);
alpha_u_mat(ETA <= alpha_eta_min) = alpha_u;
alpha_u_mat(idx_trans) = alpha_u + (1-alpha_u)*(ETA(idx_trans)-alpha_eta_min)/(alpha_eta_max-alpha_eta_min);
alpha_v_mat = ones(EXF_Nx,EXF_Ny);
alpha_v_mat(ETA <= alpha_eta_min) = alpha_v;
alpha_v_mat(idx_trans) = alpha_v + (1-alpha_v)*(ETA(idx_trans)-alpha_eta_min)/(alpha_eta_max-alpha_eta_min);


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
    uwind_mod(:,:,k-days_start+1) = uwind(:,:,k-days_start+1).*alpha_u_mat;
end
fclose(fid);
vwind = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
vwind_mod = 0*vwind;
fid = fopen(fullfile(inputfolder,mwind),'r','b');
for k=days_start:days_end
    vwind(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    vwind_mod(:,:,k-days_start+1) = vwind(:,:,k-days_start+1).*alpha_v_mat;
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