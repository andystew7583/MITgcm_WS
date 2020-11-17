%%%
%%% calcAABWcirculation.m
%%%
%%% Calculates horizontal circulation in isopycnal layers.
%%%

%%% Options
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% tmin = 19.05;
% tmax = 27.05;
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% tmin = 10.05;
% tmax = 18.05;
expname = 'hires_seq_onetwelfth_RTOPO2';
tmin = 1.05;
tmax = 9.05;

%%% Load experiment
loadexp;

%%% Select density variable in which to compute isopycnal fluxes
densvar = 'PD0';
dens_levs = layers_bounds;
Nd = length(dens_levs)-1;
p_ref = -rhoConst*gravity*RC(1)/1e4; %%% Reference pressure for surface-referenced potential density

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

%%% Partial cell heights on fine grid
hFacW_f = zeros(Nx,Ny,Nrf);
for k=1:Nr
  hFacW_f(:,:,ffac*(k-1)+1:ffac*k) = hFacW(:,:,k*ones(1,ffac));              
end
hFacS_f = zeros(Nx,Ny,Nrf);
for k=1:Nr
  hFacS_f(:,:,ffac*(k-1)+1:ffac*k) = hFacS(:,:,k*ones(1,ffac));              
end

%%% 3D horizontal grid spacing matrices
DXG_3D = repmat(DXG,[1 1 Nd]);
DYG_3D = repmat(DYG,[1 1 Nd]);

%%% Grid of actual vertical positions, accounting for partial cells
ZZ_u = zeros(Nx,Ny,Nr);
ZZ_u_f = zeros(Nx,Ny,Nrf);
ZZ_v = zeros(Nx,Ny,Nr);
ZZ_v_f = zeros(Nx,Ny,Nrf);
DZ = zeros(Nx,Ny,Nr);
DZ_f = zeros(Nx,Ny,Nrf);
PP = zeros(Nx,Ny,Nr);
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
for k=1:Nr
  DZ(:,:,k) = delR(k);
end   
for k=1:Nrf
  DZ_f(:,:,k) = delRf(k);
end   
for k=1:Nr
  PP(:,:,k) = -delR(k);
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

%%% Calculate time-averaged isopycnal fluxes and layer thicknesses
uflux_tavg = zeros(Nx,Ny,Nd);
vflux_tavg = zeros(Nx,Ny,Nd);
uthic_tavg = zeros(Nx,Ny,Nd);
vthic_tavg = zeros(Nx,Ny,Nd);
uvel_tavg = zeros(Nx,Ny,Nr);
vvel_tavg = zeros(Nx,Ny,Nr);
dens_tavg = zeros(Nx,Ny,Nr);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))] 
  
  %%% Read isopycnal fluxes
  uflux = rdmdsWrapper(fullfile(exppath,'/results/LaUH1RHO'),itersToRead(n));
  vflux = rdmdsWrapper(fullfile(exppath,'/results/LaVH1RHO'),itersToRead(n)); 
  uthic = rdmdsWrapper(fullfile(exppath,'/results/LaHw1RHO'),itersToRead(n));
  vthic = rdmdsWrapper(fullfile(exppath,'/results/LaHs1RHO'),itersToRead(n));
  if (isempty(uflux) || isempty(vflux) || isempty(uthic) || isempty(vthic))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end
  %%% Add to time averages
  uflux_tavg = uflux_tavg + uflux/Ntime;
  vflux_tavg = vflux_tavg + vflux/Ntime;
  uthic_tavg = uthic_tavg + uthic/Ntime;
  vthic_tavg = vthic_tavg + vthic/Ntime;
  
  %%% Read velocity field
  uvel  = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),itersToRead(n));
  vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),itersToRead(n));              
  if (isempty(uvel) || isempty(vvel))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end  
  %%% Add to time averages
  uvel_tavg = uvel_tavg + uvel/Ntime;
  vvel_tavg = vvel_tavg + vvel/Ntime;

  %%% Calculate density
  theta = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));
  salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT'),itersToRead(n)); 
  if (isempty(theta) || isempty(salt))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end
  dens = densjmd95(salt,theta,p_ref*ones(Nx,Ny,Nr))-1000;      
  dens_tavg = dens_tavg + dens/Ntime; %%% Add to time average      

end

%%% Compute time-mean component of isopycnal volume flux
[uflux_mean,vflux_mean] = calcIsopFluxes (...
  uvel_tavg,vvel_tavg,dens_tavg,...
  Nx,Ny,Nr,Nrf,Nd,ffac, ...
  kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
  hFacW_f,hFacS_f,DZ_f,dens_levs);





%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Construct output file name
outfname = [expname,'_AABWcirc_.mat'];

%%% Store computed data for later
save(fullfile('products',outfname), ...
  'dens_levs','times', ...
  'uflux_tavg','vflux_tavg','uthic_tavg','vthic_tavg','uflux_mean','vflux_mean');