%%%
%%% PerturbAtmosTemp.m
%%%
%%% Modifies atmospheric temperature files.
%%%

addpath ../analysis;

%%% Load base model parameters
defineGrid;

%%% Parameters defining forcing perturbations
% deltaTemp = 5; %%% Roughly SSP585 based on https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024GL112662
deltaTemp = 2.5; %%% Reduced by factor of 2 in an attempt to ~halve buoyancy flux pert. when combined with wind pert

%%% Define range of days to operate on
days_end = 3287;
days_start = 1;

%%% Load wind data
atemp = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
atemp_mod = 0*atemp;
fid = fopen(fullfile(inputfolder,aTemp),'r','b');
for k=1:days_start:days_end
    atemp(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    atemp_mod(:,:,k-days_start+1) = atemp(:,:,k-days_start+1) + deltaTemp;
end
fclose(fid);

%%% Write modified wind forcing files
writeDataset(atemp,fullfile(inputfolder,'atempfile_unmod.bin'),ieee,prec);
writeDataset(atemp_mod,fullfile(inputfolder,aTemp),ieee,prec);