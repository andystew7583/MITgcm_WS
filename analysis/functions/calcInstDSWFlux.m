%%%
%%% calcInstDSWFlux.m
%%%
%%% Calculates the flux of dense water across the mouth of the Filchner
%%% trough from instantaneous model output.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
tmin = 1.01;
tmax = 7.01;



%%% Load experiment
loadexp;

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

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
%%% For daily/12-hourly outputs
dumpStart = 1578240;
dumpStep = 86400/2/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;

%%% Indices over which to compute fluxes
ft_lat = -74.54;
ft_w_lon = -35.3;
ft_e_lon = -30.6;
ft_yidx = find(YC(1,:)<ft_lat,1,'last');
ft_w_xidx = find(XC(:,1)>ft_w_lon,1,'first');
ft_e_xidx = find(XC(:,1)>ft_e_lon,1,'first');
Nx_ft = length(ft_w_xidx:ft_e_xidx);




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

%%% Crop matrices to the region of interest
DXG_3D = DXG_3D(ft_w_xidx:ft_e_xidx,ft_yidx,:);
DYG_3D = DYG_3D(ft_w_xidx:ft_e_xidx,ft_yidx,:);
DZ_f = DZ_f(ft_w_xidx:ft_e_xidx,ft_yidx,:);



%%% Grid of actual vertical positions, accounting for partial cells
if (ffac == 1)
  
  %%% Don't need these  
  kp_v = [];
  kn_v = [];
  wn_v = [];
  wp_v = [];    
  hFacS_f = hFacS(ft_w_xidx:ft_e_xidx,ft_yidx,:);
  
else
  
  %%% Partial cell heights on fine grid
  hFacS_f = zeros(Nx,Ny,Nrf);
  for k=1:Nr
    hFacS_f(:,:,ffac*(k-1)+1:ffac*k) = hFacS(:,:,k*ones(1,ffac));              
  end
  
  ZZ_v = zeros(Nx,Ny,Nr);
  ZZ_v_f = zeros(Nx,Ny,Nrf);
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
  
  
  %%% Extract only the parts of the matrices that we need for DSW flux
  %%% calculation
  wn_v = wn_v(ft_w_xidx:ft_e_xidx,ft_yidx,:);
  wp_v = wp_v(ft_w_xidx:ft_e_xidx,ft_yidx,:);
  kn_v = kn_v(ft_w_xidx:ft_e_xidx,ft_yidx,:);
  kp_v = kp_v(ft_w_xidx:ft_e_xidx,ft_yidx,:);
  hFacS_f = hFacS_f(ft_w_xidx:ft_e_xidx,ft_yidx,:);
  ZZ_v_f = ZZ_v_f(ft_w_xidx:ft_e_xidx,ft_yidx,:);
  ZZ_v = ZZ_v(ft_w_xidx:ft_e_xidx,ft_yidx,:);

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
T_ft = zeros(Nd+1,Ntime);
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
  
  %%% Calculate density on v points via straightforward linear
  %%% interpolation
  vvel = vvel(ft_w_xidx:ft_e_xidx,ft_yidx,:);
  dens_v = 0.5* (dens(ft_w_xidx:ft_e_xidx,ft_yidx,:) + dens(ft_w_xidx:ft_e_xidx,ft_yidx-1,:)); 
  
  %%% Interpolate u, v and density onto a finer vertical grid        
  vvel_f = zeros(Nx,Ny,Nrf);
  dens_v_f = NaN*zeros(Nx,Ny,Nrf);
  if (ffac == 1)

    %%% Shortcut if fine grid resolution = coarse grid resolution    
    vvel_f = vvel;        

    dens_v_f = dens_v;

  else   

    %%% Velocity uniform throughout each coarse grid cell to preserve
    %%% mass conservation
    for k=1:Nr      
      vvel_f(:,:,ffac*(k-1)+1:ffac*k) = vvel(:,:,k*ones(1,ffac));          
    end

    %%% Linearly interpolate density
    for i=1:Nx
      for j=3:Ny-1 %%% Restrict to wet grid cells          
        dens_v_f(i,j,:) = wp_v(i,j,:).*dens_v(i,j,squeeze(kp_v(i,j,:))) + wn_v(i,j,:).*dens_v(i,j,squeeze(kn_v(i,j,:)));
      end
    end

  end      
  
  %%% Free up memory
  clear('dens_u','dens_v');

  %%% Calculate fluxes within density surfaces
  vdz = vvel_f.*hFacS_f.*DZ_f;  
  vflux = zeros(Nx_ft,1,Nd);  
  vflux(:,:,Nd) = vflux(:,:,Nd) + sum(vdz.*(dens_v_f>dens_levs(Nd)),3);
  vflux(:,:,1) = vflux(:,:,1) + sum(vdz.*(dens_v_f<=dens_levs(2)),3);
  for m=2:Nd-1  
    vflux(:,:,m) = vflux(:,:,m) + sum(vdz.*((dens_v_f>dens_levs(m)) & (dens_v_f<=dens_levs(m+1))),3);
  end   
  vflux = sum(vflux.*DXG_3D,1); %%% Multiply by grid cell widths to get transports through each density layer in m^3/s
    
  %%% Sum fluxes to obtain transports
  for m=1:Nd  
    T_ft(m,n) = sum(vflux(:,:,m:Nd),3);     
  end

end










%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Construct output file name
outfname = [expname,'_FTtrans_inst_',densvar];
outfname = [outfname,'.mat'];

%%% Store computed data for later
save(fullfile('products',outfname),'-v7.3', ...
  'ft_lat','ft_w_lon','ft_e_lon', ...
  'ft_yidx','ft_w_xidx','ft_e_xidx',... 
  'dens_levs','times','T_ft');























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