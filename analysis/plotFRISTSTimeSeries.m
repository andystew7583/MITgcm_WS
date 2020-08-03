%%%
%%% plotFRISTSTimeSeries.m
%%%
%%% Plots time series of volume-averaged FRIS cavity temperature and
%%% salinity for one of our experiments.
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
theta_avg = zeros(1,nDumps);
salt_avg = zeros(1,nDumps);
tlen = 0;

%%% Indices over which to integrate, i.e. defining the FRIS
xidx = find((XC(:,1)>-80) & (XC(:,1)<-20));
yidx = find(YC(1,:)<-75);

%%% Volume of each grid cell
VV = repmat(RAC,[1 1 Nr]) .* repmat(reshape(DRF,[1 1 Nr]),[Nx Ny 1]) .* hFacC;
FRIS_vol = sum(sum(sum(VV(xidx,yidx,:),3),2),1);

for n=1:nDumps
 
  tt(n) =  dumpIters(n)*deltaT;
  tt(n)
  
  %%% Attempt to load melt ave per month
  theta=rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));  
  salt=rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n));  
  if (isempty(theta) || isempty(salt))
    break;
  end
        
  %%% Compute volume-integrated properties
  theta = theta .* VV;
  salt = salt .* VV;
  theta_avg(n) = sum(sum(sum(theta(xidx,yidx,:),3),2),1) ./ FRIS_vol;
  salt_avg(n) = sum(sum(sum(salt(xidx,yidx,:),3),2),1) ./ FRIS_vol;
  
  %%% Increment counter
  tlen = tlen + 1;
  
end

%%% Plot the time series

figure(11);
plot(tt(1:tlen)/86400/365,theta_avg(1:tlen));

figure(12);
plot(tt(1:tlen)/86400/365,salt_avg(1:tlen));