%%%
%%% calcFRISMeltTimeSeries.m
%%%
%%% Calculates a time series of total FRIS melt rates for one of our
%%% experiments.
%%%

function [tt,SHImelt,SHImelt_mean,XC,YC,bathy,SHELFICEtopo] = calcFRISMeltTimeSeries (expdir,expname) 

  %%% This needs to be set to ensure we are using the correct output
  %%% frequency
  loadexp;
  diagfreq = diag_frequency(end);

  %%% Frequency of diagnostic output
  diagnum = 1; %%% Should be arbitrary because all output is provided at the same frequency
  dumpFreq = abs(diag_frequency(diagnum));
  nDumps = round(nTimeSteps*deltaT/dumpFreq);
  dumpIters = round(nIter0+(1:nDumps)*dumpFreq/deltaT);
  dumpIters = dumpIters(dumpIters > nIter0);

  %%% To store the result
  tt = zeros(1,nDumps);
  SHImelt = NaN*ones(1,nDumps);
  SHImelt_mean = zeros(Nx,Ny);
  tlen = 0;

  %%% Indices over which to integrate, i.e. defining the FRIS
  xidx = find(XC(:,1)<-29.9);
  yidx = find(YC(1,:)<-74.5);

  %%% TODO need to fix this
  for n=1:nDumps

    tt(n) =  dumpIters(n)*deltaT;
    tt(n)

    %%% Attempt to load melt ave per month
    SHIfwFlx=rdmdsWrapper(fullfile(exppath,'/results/SHIfwFlx'),dumpIters(n));  
    if (isempty(SHIfwFlx))
      continue;
    end
    
    %%% Mean local melt rate
    SHImelt_mean = SHImelt_mean + SHIfwFlx;
n
    %%% Compute area-integrated freshwater flux
    SHIfwFlx = SHIfwFlx .* RAC;
    SHImelt(n) =  sum(sum(SHIfwFlx(xidx,yidx)));


    %%% Increment counter
    tlen = tlen + 1;

  end

  SHImelt_mean = SHImelt_mean / tlen;
    
end
