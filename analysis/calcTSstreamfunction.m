%%%
%%% calcTSstreamfunction.m
%%%
%%% Calculates thermohaline streamfunction in T/S space.
%%%

%%% Read experiment data
loadexp;


%%% T/S grids
dT = 0.05;
dS = 0.01;
Tmin = -3;
Tmax = 2;
Smin = 33.9;
Smax = 35;
SS = Smin:dS:Smax;
TT = Tmin:dT:Tmax;
NS = length(SS);
NT = length(TT);

%%% Experiment output iterations
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round(nIter0+(1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Determine iteraton numbers to process
readIters = [];
times = [];
for n=1:length(dumpIters)
 
  tyears = dumpIters(n)*deltaT/86400/365;
 
  if ((tyears >= tmin) && (tyears <= tmax))    
    readIters = [readIters dumpIters(n)];
    times = [times dumpIters(n)*deltaT];
  end
  
end
Ntime = length(readIters);

%%% Grids
DXG = rdmds(fullfile(resultspath,'DXG'));
DYG = rdmds(fullfile(resultspath,'DYG'));
RAC = rdmds(fullfile(resultspath,'RAC'));
DRF = rdmds(fullfile(resultspath,'DRF'));
hFacS = rdmds(fullfile(resultspath,'hFacS'));
hFacW = rdmds(fullfile(resultspath,'hFacW'));
hFacC = rdmds(fullfile(resultspath,'hFacC'));

%%% Areas of cell faces
Axz = repmat(reshape(DRF,[1 1 Nr]),[Nx Ny]).*DXG.*hFacS;
Ayz = repmat(reshape(DRF,[1 1 Nr]),[Nx Ny]).*DYG.*hFacW;
Axy = RAC;

%%% To store fluxes in TS space
% F_TS = zeros(NS-1,NT,Ntime);
psi_TS = zeros(NS,NT,Ntime);
for n=1:Ntime

  %%% To keep track of progress
  tyears = readIters(n)*deltaT/86400/365

  %%% Read model output
  theta = rdmdsWrapper(fullfile(exppath,'results','THETA'),readIters(n));
  salt = rdmdsWrapper(fullfile(exppath,'results','SALT'),readIters(n));
  uvel = rdmdsWrapper(fullfile(exppath,'results','UVEL'),readIters(n));
  vvel = rdmdsWrapper(fullfile(exppath,'results','VVEL'),readIters(n));
  wvel = rdmdsWrapper(fullfile(exppath,'results','WVEL'),readIters(n));

  %%% Make sure no vertical velocities are present at solid cell faces
  wvel(hFacC==0) = 0;

  %%% Calculate midpoint salinity and temperature
  salt_u = 0.5*(salt([1:Nx],:,:)+salt([Nx 1:Nx-1],:,:));
  salt_v = 0.5*(salt(:,[1:Ny],:)+salt(:,[Ny 1:Ny-1],:));
  salt_w = 0.5*(salt(:,:,[1:Nr])+salt(:,:,[Nr 1:Nr-1]));
%   theta_u = 0.5*(theta([1:Nx],:,:)+theta([Nx 1:Nx-1],:,:));
%   theta_v = 0.5*(theta(:,[1:Ny],:)+theta(:,[Ny 1:Ny-1],:));
%   theta_w = 0.5*(theta(:,:,[1:Nr])+theta(:,:,[Nr 1:Nr-1]));

  %%% Deal with open boundary transports by wrapping them around the boundary
  %%% artificially, closing the domain-integrated transport. Iterate around
  %%% open boundary points computing flux divergence at each point, and
  %%% balancing it with a boundary-parallel barotropic flow.
  %%% BOUNDARY-MEAN FLOW %%%
  NEtransport = sum(sum(vvel(1:Nx-1,Ny,:).*Axz(1:Nx-1,Ny,:))) ...
              + sum(sum(uvel(Nx,1:Ny-1,:).*Ayz(Nx,1:Ny-1,:)));              
  NEarea = sum(sum(Axz(1:Nx-1,Ny,:))) + sum(sum(Ayz(Nx,1:Ny-1,:)));
  Ubdry = NEtransport / NEarea;
  %%% NORTHERN BOUNDARY %%%
  for i=1:Nx-1

    %%% Compute N boundary flux divergence
    fluxconv = sum(vvel(i,Ny,:).*Axz(i,Ny,:) - Ubdry.*Axz(i,Ny,:) + uvel(i,Ny,:).*Ayz(i,Ny,:) - uvel(i+1,Ny,:).*Ayz(i+1,Ny,:), 3);
    wcarea = sum(Ayz(i+1,Ny,:),3);

    %%% Avoid NaNs
    if (wcarea > 0)

      %%% Add barotropic correction to eliminate boundary flux divergence
      uanom = fluxconv / wcarea;
      uvel(i+1,Ny,:) = uvel(i+1,Ny,:) + uanom;

      %%% Recompute flux convergence in order to recompute vertical velocity
      fluxconv = vvel(i,Ny,:).*Axz(i,Ny,:) - Ubdry.*Axz(i,Ny,:) + uvel(i,Ny,:).*Ayz(i,Ny,:) - uvel(i+1,Ny,:).*Ayz(i+1,Ny,:);
      wvel(i,Ny,2:Nr) = - cumsum(fluxconv(:,:,1:Nr-1),3) / RAC(i,Ny);

    end

  end
  %%% NORTHEAST CORNER %%%
  %%% Compute NE corner flux divergence
  fluxconv = sum(vvel(Nx,Ny,:).*Axz(Nx,Ny,:) - Ubdry.*Axz(Nx,Ny,:) + uvel(Nx,Ny,:).*Ayz(Nx,Ny,:) - Ubdry.*Ayz(Nx,Ny,:), 3);

  %%% Avoid NaNs
  wcarea = sum(Axz(Nx,Ny,:),3);
  if (wcarea > 0) 

    %%% Add barotropic correction to eliminate boundary flux divergence
    vanom = fluxconv / wcarea;
    vvel(Nx,Ny,:) = vvel(Nx,Ny,:) - vanom;

    %%% Recompute flux convergence in order to recompute vertical velocity
    fluxconv = vvel(Nx,Ny,:).*Axz(Nx,Ny,:) - Ubdry.*Axz(Nx,Ny,:) + uvel(Nx,Ny,:).*Ayz(Nx,Ny,:) - Ubdry.*Ayz(Nx,Ny,:);
    wvel(Nx,Ny,2:Nr) = - cumsum(fluxconv(:,:,1:Nr-1),3) / RAC(Nx,Ny);

  end
  %%% EASTERN BOUNDARY %%%
  for j=Ny-1:-1:1

    %%% Compute E boundary flux divergence
    fluxconv = sum(vvel(Nx,j,:).*Axz(Nx,j,:) - vvel(Nx,j+1,:).*Axz(Nx,j+1,:) + uvel(Nx,j,:).*Ayz(Nx,j,:) - Ubdry.*Ayz(Nx,j,:), 3);

    %%% Avoid NaNs
    wcarea = sum(Axz(Nx,j,:),3);
    if (wcarea > 0)

      %%% Add barotropic correction to eliminate boundary flux divergence
      vanom = fluxconv / wcarea;
      vvel(Nx,j,:) = vvel(Nx,j,:) - vanom;

      %%% Recompute flux convergence in order to recompute vertical velocity
      fluxconv = vvel(Nx,j,:).*Axz(Nx,j,:) - vvel(Nx,j+1,:).*Axz(Nx,j+1,:) + uvel(Nx,j,:).*Ayz(Nx,j,:) - Ubdry.*Ayz(Nx,j,:);
      wvel(Nx,j,2:Nr) = - cumsum(fluxconv(:,:,1:Nr-1),3) / RAC(Nx,j);

    end

  end

  %%% Transport through grid cell faces
  Tu = Ayz.*uvel;
  Tv = Axz.*vvel;
  Tw = Axy.*wvel;
  
  %%% Map vertically integrated flux convergence
  % fluxconv = zeros(Nx,Ny);
  % fluxconv(1:Nx,1:Ny) = fluxconv(1:Nx,1:Ny) + sum(Tu(1:Nx,1:Ny,:) + Tv(1:Nx,1:Ny,:), 3);
  % fluxconv(1:Nx-1,1:Ny) = fluxconv(1:Nx-1,1:Ny) - sum(Tu(2:Nx,1:Ny,:),3);
  % fluxconv(1:Nx,1:Ny-1) = fluxconv(1:Nx,1:Ny-1) - sum(Tv(1:Nx,2:Ny,:),3);
  % figure(100);
  % pcolor(XC,YC,fluxconv);
  % shading interp;
  % colorbar;

%   fluxdiv3D = zeros(Nx,Ny,Nr);
%   fluxdiv3D = fluxdiv3D - Tu - Tv + Tw;
%   fluxdiv3D(1:Nx-1,1:Ny,1:Nr) = fluxdiv3D(1:Nx-1,1:Ny,1:Nr) + Tu(2:Nx,1:Ny,:);
%   fluxdiv3D(1:Nx,1:Ny-1,1:Nr) = fluxdiv3D(1:Nx,1:Ny-1,1:Nr) + Tv(1:Nx,2:Ny,:);
%   fluxdiv3D(1:Nx,1:Ny,1:Nr-1) = fluxdiv3D(1:Nx,1:Ny,1:Nr-1) - Tw(1:Nx,1:Ny,2:Nr);
%   sum(sum(sum(fluxdiv3D(theta<-1.75))))
tstart = tic();

  %%% Loop over temperature space
  for nt = 1:NT

    nt

    %%% Create mask over all physical velocity points, leaving non-zero values 
    %%% Only where adjacent temperature points straddle the current
    %%% temperature (TT(nt)), and signed +1 (-1) if the temperature gradient
    %%% is positive (negative)
    msk_u_theta =  ((theta([1:Nx],:,:)-TT(nt)) .* (theta([Nx 1:Nx-1],:,:)-TT(nt)) < 0) ...
                .* sign(theta([1:Nx],:,:)-theta([Nx 1:Nx-1],:,:)); %%% Flow crosses temp class TT(nt), positive if dT/dx>0
    msk_v_theta =  ((theta(:,[1:Ny],:)-TT(nt)) .* (theta(:,[Ny 1:Ny-1],:)-TT(nt)) < 0) ...
                .* sign(theta(:,[1:Ny],:,:)-theta(:,[Ny 1:Ny-1],:)); %%% Flow crosses temp class TT(nt), positive if dT/dy>0
    msk_w_theta =  ((theta(:,:,[Nr 1:Nr-1])-TT(nt)) .* (theta(:,:,[1:Nr])-TT(nt)) < 0) ...
                .* sign(theta(:,:,[Nr 1:Nr-1])-theta(:,:,[1:Nr])); %%% Flow crosses temp class TT(nt), positive if dT/dz>0
              
    %%% Precompute product of transport with theta mask        
    Tu_msktheta = Tu.*msk_u_theta;    
    Tv_msktheta = Tv.*msk_v_theta;
    Tw_msktheta = Tw.*msk_w_theta;

    %%% Loop over salinity space
    for ns = 1:NS-1  

      %%% Create masks over all physical velocity points that are non-zero
      %%% only where the salinity falls within the current salinity bin
%       msk_u_salt = (salt_u>SS(ns)) & (salt_u<SS(ns+1)); %%% Salinity falls in salt class [SS(ns) SS(ns+1)]     
%       msk_v_salt = (salt_v>SS(ns)) & (salt_v<SS(ns+1));
%       msk_w_salt = (salt_w>SS(ns)) & (salt_w<SS(ns+1));
      msk_u_salt = (salt_u<SS(ns)); %%% Salinity falls in salt class [SS(ns) SS(ns+1)]     
      msk_v_salt = (salt_v<SS(ns));
      msk_w_salt = (salt_w<SS(ns));

      %%% Flux in T/S space is equal to sum of physical fluxes over all 
      %%% unmasked points     
%       F_TS(ns,nt,n) = sum(sum(sum(Tu_msktheta.*msk_u_salt + Tv_msktheta.*msk_v_salt + Tw_msktheta.*msk_w_salt)));
      psi_TS(ns,nt,n) = sum(sum(sum(Tu_msktheta.*msk_u_salt + Tv_msktheta.*msk_v_salt + Tw_msktheta.*msk_w_salt)));      
%       psi_TS(ns,nt,n) = Tu_msktheta(:)'*msk_u_salt(:) + Tv_msktheta(:)'*msk_v_salt(:) + Tw_msktheta(:)'*msk_w_salt(:);

    end
  end
  
  toc(tstart)
end


%%% Streamfunction is equal to sum of fluxes
% psi_TS(2:NS,:,:) = cumsum(F_TS,1);

%%% Store computed data for later
save([expname,'_MOC_TS.mat'],'SS','TT','times','psi_TS');

%%% Make a plot
figure(2);
[TTT,SSS]=meshgrid(TT,SS);
pcolor(SSS,TTT,mean(psi_TS,3)/1e6);
shading interp;
colormap redblue
caxis([-7 7]);
colorbar;














function psi_TS = calcStreamfunction

end
