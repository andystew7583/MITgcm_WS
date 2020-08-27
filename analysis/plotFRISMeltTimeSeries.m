%%%
%%% plotFRISMeltTimeSeries.m
%%%
%%% Plots a time series of total FRIS melt rates for one of our
%%% experiments.
%%%

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
SHImelt = zeros(1,nDumps);
SHImelt_mean = zeros(Nx,Ny);
tlen = 0;

%%% Indices over which to integrate, i.e. defining the FRIS
xidx = find(XC(:,1)<-29.9);
yidx = find(YC(1,:)<-74.5);

for n=1:nDumps
 
  tt(n) =  dumpIters(n)*deltaT;
  tt(n)
  
  %%% Attempt to load melt ave per month
  SHIfwFlx=rdmdsWrapper(fullfile(exppath,'/results/SHIfwFlx'),dumpIters(n));  
  if (isempty(SHIfwFlx))
    break;
  end
  
  %%% Mean local melt rate
  SHImelt_mean = SHImelt_mean + SHIfwFlx;
        
  %%% Compute area-integrated freshwater flux
  SHIfwFlx = SHIfwFlx .* RAC;
  SHImelt(n) =  sum(sum(SHIfwFlx(xidx,yidx)));
  
  
  
  %%% Increment counter
  tlen = tlen + 1;
  
end

%%% Plot the time series
figure(10);
plot(tt(1:tlen)/86400/365,SHImelt(1:tlen)/1e12*86400*365);


SHImelt_mean = SHImelt_mean / tlen;
figure(12);
pcolor(XC,YC,-SHImelt_mean/920*86400*365);
shading interp;
colorbar;
caxis([-5 5]);
colormap redblue;
