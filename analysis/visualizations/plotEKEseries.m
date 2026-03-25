%%%
%%% plotTKE.m
%%%
%%% Plots the eddy kinetic energy output from MITgcm simulations.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Read experiment data
loadexp;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

tt = zeros(1,nDumps);
EKEavg = zeros(1,nDumps);
KElen = 0;

% for n=1:nDumps
for n=1:nDumps 

  tt(n) =  dumpIters(n)*deltaT/86400;  
  
  n
  
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),dumpIters(n));      
  uvelsq = rdmdsWrapper(fullfile(exppath,'/results/UVELSQ'),dumpIters(n)); 
  usq_eddy = uvelsq - uvel.^2;
  clear('uvel','uvelsq');
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));      
  vvelsq = rdmdsWrapper(fullfile(exppath,'/results/VVELSQ'),dumpIters(n)); 
  vsq_eddy = vvelsq - vvel.^2;
  clear('vvel','vvelsq');
        
  EKEavg(n) = sum(sum(sum(0.5.*usq_eddy.*RAW.*DRF.*hFacW))) / sum(sum(sum(RAW.*DRF.*hFacW))) ....
            + sum(sum(sum(0.5.*vsq_eddy.*RAS.*DRF.*hFacS))) / sum(sum(sum(RAS.*DRF.*hFacS)));
  
  KElen = KElen + 1;
  
end
  
figure(1);
clf;
axes('FontSize',16);
plot(tt(1:KElen)/365,EKEavg(1:KElen));
axis tight;
xlabel('t (years)');
ylabel('KE (m^2s^-^2)');


