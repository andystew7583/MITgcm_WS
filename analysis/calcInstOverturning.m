%%%
%%% calcInstOverturning.m
%%%
%%% Calculates the overturning circulation in density surfaces, similar to 
%%% that calculated using the MITgcm 'layers' package, from instantaneous
%%% model output
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
tmin = 1.01;
tmax = 7.01;

%%% Load experiment
loadexp;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Select density variable in which to compute isopycnal fluxes
densvar = 'PD0';
% densvar = 'ND1';
% densvar = 'ND2';
% densvar = 'PT';

%%% Density bins for MOC calculation  
switch (densvar)
  case 'PD0'
    % dens_levs = [30.5:1:36.5 36.6:0.1:36.8 36.9:0.01:37.4 37.42:0.02:37.6];
    dens_levs = [30.5:1:36.5 36.6:0.1:36.8 36.9:0.01:37.4 37.42:0.02:37.6] - 9.38; 
  case 'PT'
    dens_levs = [-3:.1:2];
  case 'ND1'
    dens_levs = [21:1:27 27.1:.1:27.5 27.52:.02:28.7 28.74:0.04:29.1];      
  case 'ND2'
    dens_levs = [21:1:27 27.1:.1:27.5 27.52:.02:28.7 28.74:0.04:29.1];      
end
Nd = length(dens_levs)-1;
p_ref = -rhoConst*gravity*RC(1)/1e4; %%% Reference pressure for surface-referenced potential density

%%% Define coordinate system for integrating to compute streamfunction
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
eta = -9:.1:11;
Neta = length(eta);

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
%%% For daily/12-hourly outputs
dumpStart = 1578240;
dumpStep = 86400/2/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;







%%%%%%%%%%%%%
%%% GRIDS %%%
%%%%%%%%%%%%%

%%% Create a finer vertical grid
ffac = 1;
Nrf = ffac*Nr;
delRf = zeros(1,Nrf); 
for n=1:Nr
  for m=1:ffac
    delRf((n-1)*ffac+m) = delR(n)/ffac;
  end
end
zz = - cumsum((delR + [0 delR(1:Nr-1)])/2);
zz_f = - cumsum((delRf + [0 delRf(1:Nrf-1)])/2);

%%% 3D horizontal grid spacing matrices
DXG_3D = repmat(DXG,[1 1 Nd]);
DYG_3D = repmat(DYG,[1 1 Nd]);
DZ_f = zeros(Nx,Ny,Nrf);
for k=1:Nrf
  DZ_f(:,:,k) = delRf(k);
end  


%%% Grid of actual vertical positions, accounting for partial cells
if (ffac == 1)
  
  %%% Don't need these
  kp_u = [];
  kn_u = [];
  wn_u = [];
  wp_u = [];
  kp_v = [];
  kn_v = [];
  wn_v = [];
  wp_v = [];  
  hFacW_f = [];
  hFacS_f = [];
  
else
  
  %%% Partial cell heights on fine grid
  hFacW_f = zeros(Nx,Ny,Nrf);
  for k=1:Nr
    hFacW_f(:,:,ffac*(k-1)+1:ffac*k) = hFacW(:,:,k*ones(1,ffac));              
  end
  hFacS_f = zeros(Nx,Ny,Nrf);
  for k=1:Nr
    hFacS_f(:,:,ffac*(k-1)+1:ffac*k) = hFacS(:,:,k*ones(1,ffac));              
  end
  
  ZZ_u = zeros(Nx,Ny,Nr);
  ZZ_u_f = zeros(Nx,Ny,Nrf);
  ZZ_v = zeros(Nx,Ny,Nr);
  ZZ_v_f = zeros(Nx,Ny,Nrf);
  ZZ_u(:,:,1) = - delR(1)*hFacW(:,:,1)/2;
  for k=2:Nr
    ZZ_u(:,:,k) = ZZ_u(:,:,k-1) - 0.5*delR(k-1)*hFacW(:,:,k-1) - 0.5*delR(k)*hFacW(:,:,k);
  end       
  ZZ_u_f(:,:,1) = - delRf(1)*hFacW_f(:,:,1)/2;
  for k=2:Nrf 
    ZZ_u_f(:,:,k) = ZZ_u_f(:,:,k-1) - 0.5*delRf(k-1)*hFacW_f(:,:,k-1) - 0.5*delRf(k)*hFacW_f(:,:,k);      
  end
  ZZ_v(:,:,1) = - delR(1)*hFacS(:,:,1)/2;
  for k=2:Nr
    ZZ_v(:,:,k) = ZZ_v(:,:,k-1) - 0.5*delR(k-1)*hFacS(:,:,k-1) - 0.5*delR(k)*hFacS(:,:,k);
  end       
  ZZ_v_f(:,:,1) = - delRf(1)*hFacS_f(:,:,1)/2;
  for k=2:Nrf 
    ZZ_v_f(:,:,k) = ZZ_v_f(:,:,k-1) - 0.5*delRf(k-1)*hFacS_f(:,:,k-1) - 0.5*delRf(k)*hFacS_f(:,:,k);      
  end


  %%% Matrices for vertical interpolation  
  kp_u = zeros(Nx,Ny,Nrf);
  kn_u = zeros(Nx,Ny,Nrf);
  wn_u = zeros(Nx,Ny,Nrf);
  wp_u = zeros(Nx,Ny,Nrf);
  kp_v = zeros(Nx,Ny,Nrf);
  kn_v = zeros(Nx,Ny,Nrf);
  wn_v = zeros(Nx,Ny,Nrf);
  wp_v = zeros(Nx,Ny,Nrf);
  for i=1:Nx
    for j=1:Ny

      %%% Indices of the lowest cells
      kmax_u = sum(squeeze(hFacW(i,j,:))~=0);
      kmax_u_f = ffac*kmax_u;
      kmax_v = sum(squeeze(hFacS(i,j,:))~=0);
      kmax_v_f = ffac*kmax_v;

      for k=1:Nrf

        %%% Previous and next interpolation indices
        kp_u(i,j,k) = ceil(k/ffac-0.5);
        kn_u(i,j,k) = kp_u(i,j,k) + 1;
        kp_v(i,j,k) = ceil(k/ffac-0.5);
        kn_v(i,j,k) = kp_v(i,j,k) + 1;

        %%% Fine grid cell is above highest coarse grid cell, so fine grid
        %%% gamma will just be set equal to uppermost coarse grid gamma
        if (kp_u(i,j,k) <= 0)

          kp_u(i,j,k) = 1;
          wp_u(i,j,k) = 0;
          wn_u(i,j,k) = 1;

        else

          %%% Fine grid cell is below lowest coarse grid cell, so fine grid
          %%% gamma will just be set equal to lowermost coarse grid gamma
          if (kn_u(i,j,k) > kmax_u)

            kn_u(i,j,k) = kmax_u;
            wn_u(i,j,k) = 0;
            wp_u(i,j,k) = 1;

          %%% Otherwise set weights to interpolate linearly between neighboring
          %%% coarse-grid gammas
          else

            wp_u(i,j,k) = (ZZ_u(i,j,kn_u(i,j,k))-ZZ_u_f(i,j,k))./(ZZ_u(i,j,kn_u(i,j,k))-ZZ_u(i,j,kp_u(i,j,k)));
            wn_u(i,j,k) = 1 - wp_u(i,j,k);

          end

        end

        %%% Fine grid cell is above highest coarse grid cell, so fine grid
        %%% gamma will just be set equal to uppermost coarse grid gamma
        if (kp_v(i,j,k) <= 0)

          kp_v(i,j,k) = 1;
          wp_v(i,j,k) = 0;
          wn_v(i,j,k) = 1;

        else

          %%% Fine grid cell is below lowest coarse grid cell, so fine grid
          %%% gamma will just be set equal to lowermost coarse grid gamma
          if (kn_v(i,j,k) > kmax_v)

            kn_v(i,j,k) = kmax_v;
            wn_v(i,j,k) = 0;
            wp_v(i,j,k) = 1;

          %%% Otherwise set weights to interpolate linearly between neighboring
          %%% coarse-grid gammas
          else

            wp_v(i,j,k) = (ZZ_v(i,j,kn_v(i,j,k))-ZZ_v_f(i,j,k))./(ZZ_v(i,j,kn_v(i,j,k))-ZZ_v(i,j,kp_v(i,j,k)));
            wn_v(i,j,k) = 1 - wp_v(i,j,k);

          end

        end

      end

    end
  end

  %%% Free up memory
  clear('ZZ_u','ZZ_v','ZZ_u_f','ZZ_v_f');

end









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










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ISOPYCNAL FLUX CALCULATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate time-averaged isopycnal flux, density and velocity
psi = zeros(Neta,Nd+1,Ntime);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))]      
  
  %%% Read velocity field
  uvel  = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),itersToRead(n));
  vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),itersToRead(n));              
  if (isempty(uvel) || isempty(vvel))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end

  %%% Compute/load "density" variable
  switch (densvar)

    case 'PD0'
      %%% Calculate density
      theta = rdmdsWrapper(fullfile(exppath,'/results/THETA_12hourly'),itersToRead(n));
      salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT_12hourly'),itersToRead(n)); 
      if (isempty(theta) || isempty(salt))
        ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
        break;
      end
      dens = densjmd95(salt,theta,p_ref*ones(Nx,Ny,Nr))-1000;

    case 'ND1'
     
      %%% TODO      
      dens = zeros(Nx,Ny,Nr);
      
    case 'ND2'
      
      %%% TODO
      dens = zeros(Nx,Ny,Nr);
      
    case 'PT'
      %%% 'Density' is just potential temperature
      dens = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));
      
  end  
  
  %%% Currently we can only calculate overturning using time-mean flow in
  %%% ND2 surfaces
  if (strcmp(densvar,'ND1') || strcmp(densvar,'ND2'))
    continue;
  end
  
  %%% Compute mean streamfunction
  if (ffac == 1)
    [uflux,vflux] = calcIsopFluxes (...
      uvel,vvel,dens,...
      Nx,Ny,Nr,Nrf,Nd,ffac, ...
      kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
      hFacW,hFacS,DZ_f,dens_levs);
  else
    [uflux,vflux] = calcIsopFluxes (...
      uvel,vvel,dens,...
      Nx,Ny,Nr,Nrf,Nd,ffac, ...
      kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
      hFacW_f,hFacS_f,DZ_f,dens_levs);
  end
  psi(:,:,n) = calcIsopStreamfunction(...
    uflux,vflux, ...
    Nx,Ny,Neta,Nd, ...  
    DXG_3D,DYG_3D,ETA,eta);
  clear('uflux','vflux'); 

end










%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Construct output file name
outfname = [expname,'_MOC_inst_',densvar];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];

%%% Store computed data for later
save(fullfile('products',outfname),'-v7.3', ...
  'eta','ETA','dens_levs','times', ...
  'psi');























function psi = calcIsopStreamfunction(...
  uflux,vflux, ...
  Nx,Ny,Neta,Nd, ...  
  DXG_3D,DYG_3D,ETA,eta)

  %%% Compute horizontal divergence of isopycnal fluxes
  fluxdiv = zeros(Nx,Ny,Nd);
  fluxdiv(1:Nx-1,1:Ny-1,:) = uflux(2:Nx,1:Ny-1,:) .* DYG_3D(2:Nx,1:Ny-1,:) ...
                              - uflux(1:Nx-1,1:Ny-1,:) .* DYG_3D(1:Nx-1,1:Ny-1,:) ...
                              + vflux(1:Nx-1,2:Ny,:) .* DXG_3D(1:Nx-1,2:Ny,:) ...
                              - vflux(1:Nx-1,1:Ny-1,:) .* DXG_3D(1:Nx-1,1:Ny-1,:);
                       
  %%% Integrate flux divergence across lines of constant eta (parallel to FRIS face)
  eflux = zeros(Neta,Nd);
  for m = 1:Neta
    msk = repmat(ETA<eta(m),[1 1 Nd]);
    eflux(m,:) = squeeze(sum(sum(fluxdiv.*msk,1),2));
  end

  %%% Sum fluxes to obtain streamfunction
  psi = zeros(Neta,Nd+1);
  for m=1:Nd  
    psi(:,m) = -sum(eflux(:,m:Nd),2);     
  end

end















