setExpname%%%
%%% plotTKEseries.m
%%%pri
%%% Plots the instantaneous total horizontal kinetic energy output from 
%%% MITgcm simulations.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Read experiment data
loadexp;

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(end);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
nDumps = length(dumpIters);

tt = zeros(1,nDumps);
KEtot = zeros(1,nDumps);
maxSpeed = zeros(1,nDumps );
KElen = 0;

%%% To calculate time averages
uvel_mean = zeros(Nx,Ny,Nr);
vvel_mean = zeros(Nx,Ny,Nr);
KE_mean = zeros(Nx,Ny,Nr);
mean_cntr = 0;

for n=1:nDumps
% for n=364:364
 
  tt(n) =  dumpIters(n)*deltaT/86400/365;
  tt(n)
  
  %%% Attempt to load either instantaneous velocities or their squares
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n));      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n)); 
  if (isempty(uvel) || isempty(vvel))   
    break;
  end  

  %%% Calculate kinetic energy and maximum horizontal speed
  KE = 0.5*(uvel.^2+vvel.^2);             
  maxSpeed(n) = max(max(max(2*KE)));
  
  %%% Time-averages
  if (n > 66)
    uvel_mean = uvel_mean + uvel;
    vvel_mean = vvel_mean + vvel;
    KE_mean = KE_mean + KE;
    mean_cntr = mean_cntr + 1;
  end  
  
  %%% Integrate EKE over the whole domain
  KEtot(n) = 0;
  for i=1:Nx
    for j=1:Ny
      for k=1:Nr
        KEtot(n) = KEtot(n) + KE(i,j,k)*RAW(i,j)*delR(k);
      end
    end
  end
  KElen = KElen + 1;
  
end

%%% Calculate time average
uvel_mean = uvel_mean / mean_cntr;
vvel_mean = vvel_mean / mean_cntr;
KE_mean = KE_mean / mean_cntr;
EKE_mean = KE_mean - 0.5*(uvel_mean.^2+vvel_mean.^2);
EKE_mean(hFacC==0) = NaN;



figure(3);
clf;
axes('FontSize',16);
plot(tt(1:KElen),KEtot(1:KElen),'ko-');
axis tight;
xlabel('t (years)');
ylabel('KE (m^2s^-^2)');

figure(4);
clf;
axes('FontSize',16);
plot(tt(1:KElen),sqrt(maxSpeed(1:KElen)),'ko-');
axis tight;
xlabel('t (years)');
ylabel('Max speed (m^2s^-^2)');


