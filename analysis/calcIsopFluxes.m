%%%
%%% calcIsopFluxes.m
%%%
%%% Computes vertically averaged volume fluxes in a series of density
%%% surfaces.
%%%
function [uflux,vflux] = calcIsopFluxes (...
  uvel,vvel,dens,...
  Nx,Ny,Nr,Nrf,Nd,ffac, ...
  kp_u,kn_u,wp_u,wn_u,kp_v,kn_v,wp_v,wn_v, ...
  hFacW_f,hFacS_f,DZ_f,dens_levs)

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
  
end