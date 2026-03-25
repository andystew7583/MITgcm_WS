%%%
%%% plotMeanTempSeries.m
%%%
%%% Plots the instantaneous domain-mean potential temperature from 
%%% MITgcm simulations.
%%%
%%%

%%% Read experiment data
loadexp;

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(1);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
nDumps = length(dumpIters);

tt = zeros(1,nDumps);
SIvolume = NaN*ones(1,nDumps);
ptlen = 0;
for n=1:nDumps
% for n=364:364
 
n
  tt(n) =  dumpIters(n)*deltaT/86400/365;
  
  %%% Attempt to load either instantaneous velocities or their squares
  SIheff = rdmdsWrapper(fullfile(exppath,'/results/SIheff'),dumpIters(n)) ;      
  SIarea = rdmdsWrapper(fullfile(exppath,'/results/SIarea'),dumpIters(n)) ; 
  
  if (isempty(SIheff) || isempty(SIarea))   
    continue;
  end  
  
  %%% Calculate domain-mean potential temperature and salinity
  
  SIvolume(n) = sum(sum(SIheff.*DXC.*DYC));
  
  %%% Increment counter
  ptlen = ptlen + 1;
  
end

figure(31);
plot(tt(1:ptlen),SIvolume(1:ptlen),'o-');
axis tight;
xlabel('t (years)');
ylabel('Total sea ice volume (m^3)');
