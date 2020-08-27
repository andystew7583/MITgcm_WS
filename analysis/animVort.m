%%%
%%% animVort.m
%%%
%%% Makes a movie of the vorticity.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Read experiment data
loadexp;

%%% Vertical grid spacing matrix
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

%%% Diagnostic indix corresponding to instantaneous velocity
diagnum = length(diag_frequency);

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(diagnum);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round(nIter0 + (1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
nDumps = length(dumpIters);

%%% Initialize movie
figure(1);
set(gcf,'Color','w');
M = moviein(nDumps);

%%% Loop through iterations
for n=12:nDumps
% for n=250:300
 
  tt(n) =  dumpIters(n)*deltaT/86400;
  tt(n)
  
  %%% Attempt to load either instantaneous velocities or their squares
%   uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n)) ;      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n)); 
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),dumpIters(n)) ;      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n)); 
  if (isempty(uvel) || isempty(vvel))   
    break;
  end
  
  %%% Plot the vorticity  
  vort = zeros(Nx,Ny);
  zlev = 20;
%   vort(:,2:Ny) = - (uvel(:,2:Ny,zlev)-uvel(:,1:Ny-1,zlev))./DYC(:,2:Ny);
%   vort(2:Nx,:) = vort(2:Nx,:) + (vvel(2:Nx,:,zlev)-vvel(1:Nx-1,:,zlev))./DXC(2:Nx,:);
  ubt = sum(uvel.*DZ.*hFacW,3) ./ sum(DZ.*hFacW,3);
  vbt = sum(vvel.*DZ.*hFacS,3) ./ sum(DZ.*hFacS,3);
  vort(:,2:Ny) = - (ubt(:,2:Ny)-ubt(:,1:Ny-1))./DYC(:,2:Ny);
  vort = vort + (vbt([2:Nx 1],:)-vbt(:,:))./DXC; 
  Omega = 2*pi*366/365/86400;
  ff = 2*Omega*sind(YG);
  pcolor(XG,YG,vort./abs(ff));
  shading interp;
%   contourf(XX/1000,YY/1000,vort/f0,-1:0.05:1,'EdgeColor','None');
  colorbar;
  colormap redblue;
%   caxis([-.5 .5]);
%   caxis([-2 2]);
  caxis([-.1 .1]);
  set(gca,'FontSize',16);
  xlabel('Longitude','interpreter','latex');
  ylabel('Latitude','interpreter','latex');
  title(['\zeta/f at t= ',num2str(round(tt(n)),'%3d'),' days']);
%   title(['\zeta/f at t= ',num2str(tt(n)/365,'%.1f'),' years']);
  M(n) = getframe(gcf);  
   
end