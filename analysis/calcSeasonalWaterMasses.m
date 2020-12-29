%%%
%%% calcSeasonalWaterMasses.m
%%%
%%% Computes seasonally-averaged T/S stratification.
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
expname = 'hires_seq_onetwelfth_RTOPO2';
tmin = 0.05;
tmax = 9.05;
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% 3D grid spacing matrices
DXG_3D = repmat(DXG,[1 1 Nr]);
DYG_3D = repmat(DYG,[1 1 Nr]);
DRF_3D = repmat(DRF,[Nx Ny 1]);







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











%%%%%%%%%%%%%%%%%%%%
%%% DO AVERAGING %%%
%%%%%%%%%%%%%%%%%%%%

%%% DJF
theta_djf = zeros(Nx,Ny,Nr);
salt_djf = zeros(Nx,Ny,Nr);
navg = 0;
for nyear=2:9
  for nmonth = 0:1:2
    
    n = (nyear-1)*12+nmonth

    %%% Print current time to keep track of calculation
    tyears = itersToRead(n)*deltaT/86400/365; 
    [num2str(tyears) num2str(itersToRead(n))]      

    %%% Read velocity field
    theta  = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));
    salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT'),itersToRead(n));

    if (isempty(theta) || isempty(salt))
      ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
      break;
    end
    
    theta_djf = theta_djf + theta;    
    salt_djf = salt_djf + salt;
    navg = navg + 1;
  
  end
end
theta_djf = theta_djf /navg;
salt_djf = salt_djf / navg;


%%% JJA
theta_jja = zeros(Nx,Ny,Nr);
salt_jja = zeros(Nx,Ny,Nr);
navg = 0;
for nyear=2:9
  for nmonth = 6:1:8
    
    n = (nyear-1)*12+nmonth

    %%% Print current time to keep track of calculation
    tyears = itersToRead(n)*deltaT/86400/365; 
    [num2str(tyears) num2str(itersToRead(n))]      

    %%% Read velocity field
    theta  = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));
    salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT'),itersToRead(n));

    if (isempty(theta) || isempty(salt))
      ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
      break;
    end
    
    theta_jja = theta_jja + theta;    
    salt_jja = salt_jja + salt;
    navg = navg + 1;
  
  end
end
theta_jja = theta_jja /navg;
salt_jja = salt_jja / navg;





%%% Store computed data for later
outfname = [expname,'_SeasonalStrat'];
outfname = [outfname,'.mat'];
save(fullfile('products',outfname), ...
  'theta_djf','theta_jja','salt_djf','salt_jja');
