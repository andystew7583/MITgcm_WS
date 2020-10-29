%%%
%%% calcSurfFluxSeries.m
%%% 
%%% Calculates time series of surface heat and salt fluxes.
%%%





%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% tmin = 19.05;
% tmax = 27.05;
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% tmin = 10.05;
% tmax = 18.05;
expname = 'hires_seq_onetwelfth_notides_RTOPO2';
tmin = 1.05;
tmax = 9.05;
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRE-DETERMINE ITERATION NUMBERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determine iteration numbers to process
itersToRead = [];
times = [];
for n=1:length(dumpIters)
 
  tyears = dumpIters(n)*deltaT/86400/365;
 
  if ((tyears >= tmin) && (tyears <= tmax))    
    itersToRead = [itersToRead dumpIters(n)];
    times = [times dumpIters(n)*deltaT];
  end
  
end
Ntime = length(itersToRead);











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HEAT/SALT FLUX CALCULATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tflux = zeros(Nx,Ny,Ntime);
sflux = zeros(Nx,Ny,Ntime);
SHIfwFlx = zeros(Nx,Ny,Ntime);
SHIhtFlx = zeros(Nx,Ny,Ntime);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))]      
  
  %%% Read fluxes
  tflux(:,:,n)  = rdmdsWrapper(fullfile(exppath,'/results/TFLUX'),itersToRead(n));
  sflux(:,:,n)  = rdmdsWrapper(fullfile(exppath,'/results/SFLUX'),itersToRead(n));
  SHIfwFlx(:,:,n)  = rdmdsWrapper(fullfile(exppath,'/results/SHIfwFlx'),itersToRead(n));
  SHIhtFlx(:,:,n)  = rdmdsWrapper(fullfile(exppath,'/results/SHIhtFlx'),itersToRead(n));
  
end

%%% Store computed data for later
outfname = [expname,'_surfFluxes'];
outfname = [outfname,'.mat'];
save(fullfile('products',outfname), ...
  'tflux','sflux','SHIfwFlx','SHIhtFlx');

