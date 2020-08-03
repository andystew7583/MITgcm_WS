%%%
%%% PerturbWindForcing.m
%%%
%%% Modifies wind forcing files to increase winds at the front of the FRIS.
%%%

%%% Load base model parameters
defineGrid;

%%%%%%%%%% alpha
%%% choose max amount by which to perturb winds.
alpha_max = 1.5;
Rlat = 4; %%% Exponential decay scale 
Rlon = 16;
alpha_lon = -60;
alpha_lat = -74.5;

%%% Max componentwise wind speed
wind_max = 30;

%%% Define range of days to operate on
days_end = 3287;
days_start = 1;

%%% Load wind data
uwind = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,zwind),'r','b');
for k=1:days_start:days_end
    uwind(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid);
vwind = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
fid = fopen(fullfile(inputfolder,mwind),'r','b');
for k=days_start:days_end
    vwind(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
end
fclose(fid);

%%% Adding Wind Perturbation and limiting magnitudes
for k=1:EXF_Nx
  for p = 1:EXF_Ny

    alpha = 1 + (alpha_max-1) * exp(-((EXF_XMC(p,k)-alpha_lon)/Rlon).^2-((EXF_YMC(p,k)-alpha_lat)/Rlat).^2);
    uwind(k,p,:) = uwind(k,p,:)*alpha;
    vwind(k,p,:) = vwind(k,p,:)*alpha;

  end
    
end

%%% Limit wind magnitudes
uwind(uwind>wind_max) = wind_max;
uwind(uwind<-wind_max) = -wind_max;
vwind(vwind>wind_max) = wind_max;
vwind(vwind<-wind_max) = -wind_max;

%%% Write modified wind forcing files
writeDataset(uwind,fullfile(inputfolder,zwind),ieee,prec);
writeDataset(vwind,fullfile(inputfolder,mwind),ieee,prec);