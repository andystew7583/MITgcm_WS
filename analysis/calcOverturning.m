%%%
%%% calcOverturning.m
%%%
%%% Calculates the overturning circulation, calculated using the MITgcm 
%%% 'layers' package.
%%%

%%% Load experiment
loadexp;

%%% Density bins for MOC calculation  
dens_levs = [30.5:1:36.5 36.6:0.1:36.8 36.9:0.01:37.4 37.42:0.02:37.6];
Nd = length(dens_levs)-1;
p_ref = 2000;

%%% Define coordinate system for integrating to compute streamfunction
lat1 = -74.8;
lon1 = -61; 
lat2 = -78.35;
lon2 = -36.7;
phic = atan((lat2-lat1)/(lon2-lon1));
CHI = (XC-lon1)*cos(phic) + (YC-lat1)*sin(phic);
ETA = (YC-lat1)*cos(phic) - (XC-lon1)*sin(phic);
pcolor(XC,YC,ETA)
eta = -9:.1:11;
Neta = length(eta);

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
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










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ISOPYCNAL FLUX CALCULATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate time-averaged isopycnal flux, density and velocity
psi = zeros(Neta,Nd+1,Ntime);
uvel_tavg = zeros(Nx,Ny,Nr);
vvel_tavg = zeros(Nx,Ny,Nr);
dens_tavg = zeros(Nx,Ny,Nr);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = readIters(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(dumpIters(n))]      
  
  %%% Read velocity field
  uvel  = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),dumpIters(n));
  vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));              
  if (isempty(uvel) || isempty(vvel))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end

  %%% TODO calculate or read density

  %%% Calculate density
  theta = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));
  salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n));    
  dens = densjmd95(salt,theta,p_ref*ones(Nx,Ny,Nr))-1000;


  %%% Compute streamfunction
  psi(:,:,n) = calcIsopStreamfunction(...
    uvel,vvel,dens,...
    Nx,Ny,Nr,Nrf,Neta,Nd,ffac, ...
    kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
    DXG_3D,DYG_3D,hFacW_f,hFacS_f,DZ_f,dens_levs,ETA,eta);   
  

  %%% Add to time averages
  uvel_tavg = uvel_tavg + uvel;
  vvel_tavg = vvel_tavg + vvel;
  dens_tavg = dens_tavg + dens;

end

%%% Divide by Ntime to get time mean
uvel_tavg = uvel_tavg / Ntime;
vvel_tavg = vvel_tavg / Ntime;
dens_tavg = dens_tavg / Ntime;


%%% Compute ``mean'' streamfunction
psi_mean = calcIsopStreamfunction(...
    uvel_tavg,vvel_tavg,dens_tavg,...
    Nx,Ny,Nr,Nrf,Neta,Nd,ffac, ...
    kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
    DXG_3D,DYG_3D,hFacW_f,hFacS_f,DZ_f,dens_levs,ETA,eta);

psi_eddy = psi-psi_mean;








%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Store computed data for later
save([expname,'_MOC_dens.mat'],'eta','ETA','CHI','lon1','lat1','lon2','lat2', ...
                                'dens_levs','times','psi','psi_mean','psi_eddy');




figure(1);
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,mean(psi,3)/1e6);
% pcolor(EE,DD,mean(psi_mean,3)/1e6);
shading interp;
caxis([-7 7]);
set(gca,'YDir','reverse');
set(gca,'YLim',[36.8 37.4]);
colormap redblue(56);
colorbar;














%%%
%%% Functionalized to avoid repeated code
%%%
function psi = calcIsopStreamfunction(...
  uvel,vvel,dens,...
  Nx,Ny,Nr,Nrf,Neta,Nd,ffac, ...
  kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
  DXG_3D,DYG_3D,hFacW_f,hFacS_f,DZ_f,dens_levs,ETA,eta)

  %%% Calculate density on u and v points via straightforward linear
  %%% interpolation
  dens_u = 0.5* (dens(1:Nx,:,:) + dens([2:Nx 1],:,:));
  dens_v = 0.5* (dens(:,1:Ny,:) + dens(:,[2:Ny 1],:)); 
  
  %%% Interpolate u, v and density onto a finer vertical grid      
  uvel_f = zeros(Nx,Ny,Nrf);
  dens_u_f = NaN*zeros(Nx,Ny,Nrf);
  vvel_f = zeros(Nx,Ny,Nrf);
  dens_v_f = NaN*zeros(Nx,Ny,Nrf);
  if (ffac == 1)

    %%% Shortcut if fine grid resolution = coarse grid resolution
    uvel_f = uvel;        
    dens_u_f = dens_u;
    vvel_f = vvel;        
    dens_v_f = dens_v;

  else   

    %%% Velocity uniform throughout each coarse grid cell to preserve
    %%% mass conservation
    for k=1:Nr
      uvel_f(:,:,ffac*(k-1)+1:ffac*k) = uvel(:,:,k*ones(1,ffac)); 
      vvel_f(:,:,ffac*(k-1)+1:ffac*k) = vvel(:,:,k*ones(1,ffac));          
    end

    %%% Linearly interpolate density
    for i=1:Nx
      for j=3:Ny-1 %%% Restrict to wet grid cells  
        dens_u_f(i,j,:) = wp_u(i,j,:).*dens_u(i,j,squeeze(kp_u(i,j,:))) + wn_u(i,j,:).*dens_u(i,j,squeeze(kn_u(i,j,:)));
        dens_v_f(i,j,:) = wp_v(i,j,:).*dens_v(i,j,squeeze(kp_v(i,j,:))) + wn_v(i,j,:).*dens_v(i,j,squeeze(kn_v(i,j,:)));
      end
    end

  end            

  %%% Calculate fluxes within density surfaces
  udz = uvel_f.*hFacW_f.*DZ_f;
  vdz = vvel_f.*hFacS_f.*DZ_f;
  uflux = zeros(Nx,Ny,Nd);
  vflux = zeros(Nx,Ny,Nd);
  uflux(:,:,Nd) = uflux(:,:,Nd) + sum(udz.*(dens_u_f>dens_levs(Nd)),3);
  uflux(:,:,1) = uflux(:,:,1) + sum(udz.*(dens_u_f<=dens_levs(2)),3);
  vflux(:,:,Nd) = vflux(:,:,Nd) + sum(vdz.*(dens_v_f>dens_levs(Nd)),3);
  vflux(:,:,1) = vflux(:,:,1) + sum(vdz.*(dens_v_f<=dens_levs(2)),3);
  for m=2:Nd-1
    uflux(:,:,m) = uflux(:,:,m) + sum(udz.*((dens_u_f>dens_levs(m)) & (dens_u_f<=dens_levs(m+1))),3);
    vflux(:,:,m) = vflux(:,:,m) + sum(vdz.*((dens_v_f>dens_levs(m)) & (dens_v_f<=dens_levs(m+1))),3);
  end   

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