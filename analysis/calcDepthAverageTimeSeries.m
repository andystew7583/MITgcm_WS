%%%
%%% calcDepthAverageTimeSeries.m
%%%
%%% Computes time series of depth-integrated quantities.
%%%

%%% Read experiment data
loadexp;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Storage
tt = zeros(1,nDumps);
usq_eddy_int = zeros(Nx,Ny,nDumps);
vsq_eddy_int = zeros(Nx,Ny,nDumps);
salt_int = zeros(Nx,Ny,nDumps);
tlen = 0;

% for n=1:nDumps
for n=1:nDumps 

  tt(n) =  dumpIters(n)*deltaT/86400;  
  
  n
  
  %%% Depth-integrate zonal component of EKE
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),dumpIters(n));      
  uvelsq = rdmdsWrapper(fullfile(exppath,'/results/UVELSQ'),dumpIters(n)); 
  usq_eddy = uvelsq - uvel.^2;
  clear('uvel','uvelsq');
  usq_eddy_int(:,:,n) = sum(0.5.*usq_eddy.*DRF.*hFacW,3);
  clear('usq_eddy');
  
  %%% Depth-integrate meridional component of EKE
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));      
  vvelsq = rdmdsWrapper(fullfile(exppath,'/results/VVELSQ'),dumpIters(n)); 
  vsq_eddy = vvelsq - vvel.^2;
  clear('vvel','vvelsq');         
  vsq_eddy_int(:,:,n) = sum(0.5.*vsq_eddy.*DRF.*hFacS,3);
  clear('vsq_eddy');
  
  %%% Depth-integrate salinity
  salt = rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n));      
  salt_int(:,:,n) = sum(0.5.*salt.*DRF.*hFacC,3);
  clear('salt');
  
  tlen = tlen + 1;
  
end

%%% Write to output file
outfname = [expname,'_DepthIntegrals.mat'];
save(fullfile('products',outfname), ...
    'tt','usq_eddy_int','vsq_eddy_int','salt_int',...
    '-append','-v7.3');