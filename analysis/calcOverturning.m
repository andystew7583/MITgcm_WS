%%%
%%% calcOverturning.m
%%%
%%% Calculates the overturning circulation in density surfaces, similar to 
%%% that calculated using the MITgcm 'layers' package.
%%%

%%% Options
expdir = '../experiments';
% expname = 'hires_seq_onethird_RTOPO2';
% tmin = 18.05;
% tmax = 27.05;
% expname = 'hires_seq_onesixth_RTOPO2';
% tmin = 9.05;
% tmax = 18.05;
% expname = 'hires_seq_onetwelfth_RTOPO2';
% tmin = 1.05;
% tmax = 9.05;
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
tmin = 1.01;
tmax = 7.01;
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2_SSH';
% tmin = 0.01;
% tmax = 1.01;

%%% Load experiment
loadexp;

%%% Set true to compute eddy-induced transports
calc_psi_eddy = true;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to use output from the Layers package to calculate isopycnal
%%% fluxes. N.B. if this option is selected then the density variable must
%%% be 'PD0' (surface-referenced potential density)
use_layers = true;

%%% Select density variable in which to compute isopycnal fluxes
densvar = 'PD0';
% densvar = 'ND1';
% densvar = 'ND2';
% densvar = 'PT';

%%% Density bins for MOC calculation  
if (use_layers)
  densvar = 'PD0';
  dens_levs = layers_bounds;
else
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
end
Nd = length(dens_levs)-1;
p_ref = -rhoConst*gravity*RC(1)/1e4; %%% Reference pressure for surface-referenced potential density

%%% Define coordinate system for integrating to compute streamfunction
if (use_PsiBT)

  infname = [expname,'_TSfluxes'];
  load(fullfile('products',infname),'uvel_tavg');

  %%% Calculate depth-averaged zonal velocity
  UU = sum(uvel_tavg.*repmat(DRF,[Nx Ny 1]).*hFacW,3);
  clear('uvel_tavg');
  
  %%% Calculate barotropic streamfunction
  Psi = zeros(Nx+1,Ny+1);
  Psi(2:Nx+1,2:Ny+1) = -cumsum(UU.*DYG,2);
  Psi = Psi(1:Nx,1:Ny);
  
  %%% Interpolate to cell centers
  ETA = 0.25*(Psi(1:Nx,1:Ny)+Psi([2:Nx 1],1:Ny)+Psi(1:Nx,[2:Ny 1])+Psi([2:Nx 1],[2:Ny 1]))/1e6;
  
  %%% Streamunction grid for flux calculation
  eta = -2:.1:10;
  Neta = length(eta);

else

  ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
  eta = -9:.1:11;
  Neta = length(eta);

end

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);







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
psi_mean = zeros(Neta,Nd+1,Ntime);
psi_eddy = zeros(Neta,Nd+1,Ntime);
dens_tavg = zeros(Nx,Ny,Nr);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))]      
  
  %%% Read velocity field
  uvel  = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),itersToRead(n));
  vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),itersToRead(n));              
  if (isempty(uvel) || isempty(vvel))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end

  %%% Compute/load "density" variable
  switch (densvar)

    case 'PD0'
      %%% Calculate density
      theta = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));
      salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT'),itersToRead(n)); 
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
  dens_tavg = dens_tavg + dens/Ntime; %%% Add to time average
  
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
  psi_mean(:,:,n) = calcIsopStreamfunction(...
    uflux,vflux, ...
    Nx,Ny,Neta,Nd, ...  
    DXG_3D,DYG_3D,ETA,eta);
  clear('uflux','vflux');
  
  %%% Calculate eddy streamfunction, if required
  if (calc_psi_eddy)
   
    if (use_layers)
    
      %%% Read isopycnal fluxes
      uflux  = rdmdsWrapper(fullfile(exppath,'/results/LaUH1RHO'),itersToRead(n));
      vflux  = rdmdsWrapper(fullfile(exppath,'/results/LaVH1RHO'),itersToRead(n));              
      if (isempty(uflux) || isempty(vflux))
        ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
        break;
      end

      %%% Compute total isopycnal streamfunction
      psi_tot = calcIsopStreamfunction(...
        uflux,vflux, ...
        Nx,Ny,Neta,Nd, ...  
        DXG_3D,DYG_3D,ETA,eta);
      
      %%% Eddy component is difference between total and mean
      psi_eddy(:,:,n) = psi_tot - psi_mean(:,:,n);
      
    else

      %%% Load heat and salt fluxes
      wvel  = rdmdsWrapper(fullfile(exppath,'/results/WVEL'),itersToRead(n));      
      uvelth = rdmdsWrapper(fullfile(exppath,'/results/UVELTH'),itersToRead(n));
      vvelth = rdmdsWrapper(fullfile(exppath,'/results/VVELTH'),itersToRead(n));
      wvelth = rdmdsWrapper(fullfile(exppath,'/results/WVELTH'),itersToRead(n));
      uvelslt = rdmdsWrapper(fullfile(exppath,'/results/UVELSLT'),itersToRead(n));
      vvelslt = rdmdsWrapper(fullfile(exppath,'/results/VVELSLT'),itersToRead(n));
      wvelslt = rdmdsWrapper(fullfile(exppath,'/results/WVELSLT'),itersToRead(n));

      %%% Compute eddy-induced velocities
      [u_eddy,v_eddy,w_eddy] = calcTRMvelocity (...
        uvel,vvel,wvel,theta,salt, ...
        uvelth,vvelth,wvelth, ...
        uvelslt,vvelslt,wvelslt, ...
        hFacC,hFacW,hFacS, ...
        DXG,DYG,RAC,DXC,DYC, ...
        DRF,DRC,RC,RF,...
        rhoConst,gravity);
      clear('uvel','vvel','wvel','theta','salt','uvelth','vvelth','wvelth','uvelslt','vvelslt','wvelslt','w_eddy');

      %%% Compute eddy-induced streamfunction
      if (ffac == 1)
        [uflux,vflux] = calcIsopFluxes (...
          u_eddy,v_eddy,dens,...
          Nx,Ny,Nr,Nrf,Nd,ffac, ...
          kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
          hFacW,hFacS,DZ_f,dens_levs);
      else
        [uflux,vflux] = calcIsopFluxes (...
          u_eddy,v_eddy,dens,...
          Nx,Ny,Nr,Nrf,Nd,ffac, ...
          kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
          hFacW_f,hFacS_f,DZ_f,dens_levs);
      end
      clear('u_eddy','v_eddy','dens');
      psi_eddy(:,:,n) = calcIsopStreamfunction(...
        uflux,vflux, ...
        Nx,Ny,Neta,Nd, ...  
        DXG_3D,DYG_3D,ETA,eta);
      clear('uflux','vflux');
    
    end
    
  end


end

%%% Currently we can only calculate overturning using time-mean flow in
%%% ND1 surfaces
if (strcmp(densvar,'ND1'))
  densfname = fullfile('products',[expname,'_ND1.mat']);
  load(densfname,'gg_ref','pt_ref','ss_ref');
  theta_tavg = pt_ref;
  clear('pt_ref');  
  salt_tavg = ss_ref;
  clear('ss_ref');
  dens_tavg = gg_ref;
  clear('gg_ref');
  dens_tavg(isnan(dens_tavg)) = 0;
end

%%% Read time-averaged variables
uvel_tavg = readIters(exppath,'UVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
vvel_tavg = readIters(exppath,'VVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);

%%% Partition mean streamfunction in to standing and fluctuating components
if (ffac == 1)
  [uflux,vflux] = calcIsopFluxes (...
    uvel_tavg,vvel_tavg,dens_tavg,...
    Nx,Ny,Nr,Nrf,Nd,ffac, ...
    kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
    hFacW,hFacS,DZ_f,dens_levs);
else
  [uflux,vflux] = calcIsopFluxes (...
    uvel_tavg,vvel_tavg,dens_tavg,...
    Nx,Ny,Nr,Nrf,Nd,ffac, ...
    kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
    hFacW_f,hFacS_f,DZ_f,dens_levs);
end
psi_mean_stand = calcIsopStreamfunction(...
  uflux,vflux, ...
  Nx,Ny,Neta,Nd, ...  
  DXG_3D,DYG_3D,ETA,eta);
psi_mean_fluc = mean(psi_mean,3)-psi_mean_stand;

%%% If we're using ND1 as the density variable, estimate the eddy
%%% streamfunction using the multi-annual mean fluxes
if (calc_psi_eddy && strcmp(densvar,'ND1'))
      
  %%% Load time-mean fluxes and properties
  wvel_tavg = readIters(exppath,'WVEL',dumpIters,deltaT,tmin*86400*365,tmax*86400*365,Nx,Ny,Nr);
  uvelth_tavg = readIters(exppath,'UVELTH',dumpIters,deltaT,tmin*86400*365,tmax*86400*365,Nx,Ny,Nr);
  vvelth_tavg = readIters(exppath,'VVELTH',dumpIters,deltaT,tmin*86400*365,tmax*86400*365,Nx,Ny,Nr);
  wvelth_tavg = readIters(exppath,'WVELTH',dumpIters,deltaT,tmin*86400*365,tmax*86400*365,Nx,Ny,Nr);
  uvelslt_tavg = readIters(exppath,'UVELSLT',dumpIters,deltaT,tmin*86400*365,tmax*86400*365,Nx,Ny,Nr);
  vvelslt_tavg = readIters(exppath,'VVELSLT',dumpIters,deltaT,tmin*86400*365,tmax*86400*365,Nx,Ny,Nr);
  wvelslt_tavg = readIters(exppath,'WVELSLT',dumpIters,deltaT,tmin*86400*365,tmax*86400*365,Nx,Ny,Nr);

  %%% Compute eddy-induced velocities
  [u_eddy,v_eddy,w_eddy] = calcTRMvelocity (...
    uvel_tavg,vvel_tavg,wvel_tavg,theta_tavg,salt_tavg, ...
    uvelth_tavg,vvelth_tavg,wvelth_tavg, ...
    uvelslt_tavg,vvelslt_tavg,wvelslt_tavg, ...
    hFacC,hFacW,hFacS, ...
    DXG,DYG,RAC,DXC,DYC, ...
    DRF,DRC,RC,RF,...
    rhoConst,gravity);

  %%% Compute eddy-induced streamfunction
  if (ffac == 1)
    [uflux,vflux] = calcIsopFluxes (...
      u_eddy,v_eddy,dens_tavg,...
      Nx,Ny,Nr,Nrf,Nd,ffac, ...
      kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
      hFacW,hFacS,DZ_f,dens_levs);
  else
    [uflux,vflux] = calcIsopFluxes (...
      u_eddy,v_eddy,dens_tavg,...
      Nx,Ny,Nr,Nrf,Nd,ffac, ...
      kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
      hFacW_f,hFacS_f,DZ_f,dens_levs);
  end
  
  %%% Assign to 'fluctuating' component of the streamfunction as this
  %%% includes all temporal fluctuations from the multi-annual mean
  psi_mean_fluc = calcIsopStreamfunction(...
    uflux,vflux, ...
    Nx,Ny,Neta,Nd, ...  
    DXG_3D,DYG_3D,ETA,eta);
  
end










%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Construct output file name
outfname = [expname,'_MOC_',densvar];
if (calc_psi_eddy)
  if (use_layers)
    estr = '_layers';
  else
    estr = '_TRM';
  end
else
  estr = '_noeddy';
end
outfname = [outfname,estr];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else
  if (deform_cavity)
    outfname = [outfname,'_deform'];
  end
end
outfname = [outfname,'.mat'];

%%% Store computed data for later
save(fullfile('products',outfname),'-v7.3', ...
  'eta','ETA','dens_levs','times', ...
  'psi_mean','psi_mean_stand','psi_mean_fluc', ...
  'uvel_tavg','vvel_tavg','dens_tavg');

if (calc_psi_eddy)
  save(fullfile('products',outfname),'psi_eddy', ...
  ... %'psi_eddy_stand','psi_eddy_fluc',
  '-append');
end






















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















