%%%
%%% plotTKE.m
%%%
%%% Plots the total kinetic energy output from MITgcm simulations.
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
KEtot = zeros(1,nDumps);
KElen = 0;

% for n=1:nDumps
for n=1:nDumps 

  tt(n) =  dumpIters(n)*deltaT/86400;  
  
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n));      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n));      
%   wvel = rdmdsWrapper(fullfile(exppath,'/results/WVEL_inst'),dumpIters(n));      
  uvelsq = rdmdsWrapper(fullfile(exppath,'/results/UVELSQ_inst'),dumpIters(n));
  %%%square of zonal componant of velocity
  vvelsq = rdmdsWrapper(fullfile(exppath,'/results/VVELSQ_inst'),dumpIters(n));
  %%%square of meridional componant of velocity
  wvelsq = rdmdsWrapper(fullfile(exppath,'/results/WVELSQ_inst'),dumpIters(n));
  %%%square of vertical componant of velocity
  
  if (isempty(uvelsq) || isempty(vvelsq) || isempty(wvelsq))
    break;
  end
    
%   KE = 0.5*(uvelsq + vvelsq + wvelsq - uvel.^2 - vvel.^2 - wvel.^2);
  KE = 0.5*(uvelsq + vvelsq + wvelsq);
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
  
figure(1);
clf;
axes('FontSize',16);
plot(tt(1:KElen)/365,KEtot(1:KElen));
axis tight;
xlabel('t (years)');
ylabel('KE (m^2s^-^2)');


