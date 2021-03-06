%%%
%%% calcTRMvelocity
%%%
%%% Computes the eddy-induced velocity at each point via the TRM
%%% streamfunction. All inputs are standard MITgcm outputs. N.B. the
%%% code expects RC and RF to have dimensions of [1 1 Nr] and [1 1 Nr+1]
%%% respectively.
%%%
function [u_eddy,v_eddy,w_eddy] = calcTRMvelocity (...
  uvel,vvel,wvel,theta,salt, ...
  uvelth,vvelth,wvelth, ...
  uvelslt,vvelslt,wvelslt, ...
  hFacC,hFacW,hFacS, ...
  DXG,DYG,RAC,DXC,DYC, ...
  DRF,DRC,RC,RF,...
  rhoConst,gravity)

  %%% Reference stratification to regularize the TRM where stratification is weak
  %%% N.B. This differs from the actual stratification N^2 by a factor of g
  dbuoy_dz_ref = 1e-8;

  %%% Grid sizes
  Nx = size(hFacC,1);
  Ny = size(hFacC,2);
  Nr = size(hFacC,3);

  %%% Remove dry grid cells. Should ensure that streamfunction only gets
  %%% calculated at points surrouned by wet cells
  salt(hFacC==0) = NaN;
  theta(hFacC==0) = NaN;

  %%% Calculate midpoint salinity and temperature
  salt_u = 0.5*(salt([1:Nx],:,:)+salt([Nx 1:Nx-1],:,:));
  salt_v = 0.5*(salt(:,[1:Ny],:)+salt(:,[Ny 1:Ny-1],:));
  salt_w = NaN*ones(Nx,Ny,Nr+1);
  salt_w(:,:,2:Nr) = 0.5*(salt(:,:,1:Nr-1)+salt(:,:,2:Nr));
  theta_u = 0.5*(theta([1:Nx],:,:)+theta([Nx 1:Nx-1],:,:));
  theta_v = 0.5*(theta(:,[1:Ny],:)+theta(:,[Ny 1:Ny-1],:));
  theta_w = NaN*ones(Nx,Ny,Nr+1);
  theta_w(:,:,2:Nr) = 0.5*(theta(:,:,1:Nr-1)+theta(:,:,2:Nr));
 
  %%% Compute thermal expansion and haline contraction coefficients
  press_c = -rhoConst*gravity*repmat(RC,[Nx Ny 1])/1e4; %%% N.B. Units in dbar
%   press_c = -rhoConst*gravity*repmat(RC(1),[Nx Ny Nr])/1e4; %%% N.B. Units in dbar
%   press_w = -rhoConst*gravity*repmat(RC(1),[Nx Ny Nr+1])/1e4;
  [alpha_u,beta_u] = calcAlphaBeta(salt_u,theta_u,press_c);
  [alpha_v,beta_v] = calcAlphaBeta(salt_v,theta_v,press_c);
  clear('press_c');
  press_w = -rhoConst*gravity*repmat(RF,[Nx Ny 1])/1e4;
  [alpha_w,beta_w] = calcAlphaBeta(salt_w,theta_w,press_w);
  clear('press_w');
 
  %%% Compute eddy heat and salt fluxes on cell faces
  uvelslt_eddy = uvelslt - uvel .* salt_u;
  vvelslt_eddy = vvelslt - vvel .* salt_v;
  wvelslt_eddy = NaN*ones(Nx,Ny,Nr+1);
  wvelslt_eddy(:,:,1:Nr) = wvelslt - wvel .* salt_w(:,:,1:Nr);
  uvelth_eddy = uvelth - uvel .* theta_u;
  vvelth_eddy = vvelth - vvel .* theta_v;
  wvelth_eddy = NaN*ones(Nx,Ny,Nr+1);
  wvelth_eddy(:,:,1:Nr) = wvelth - wvel .* theta_w(:,:,1:Nr);
  clear('salt_u','salt_v','salt_w','theta_u','theta_v','theta_w');
  
  %%% Compute eddy 'buoyancy' fluxes on cell faces
  uvelbuoy_eddy_u = alpha_u.*uvelth_eddy - beta_u.*uvelslt_eddy;
  vvelbuoy_eddy_v = alpha_v.*vvelth_eddy - beta_v.*vvelslt_eddy;
  wvelbuoy_eddy_w = alpha_w.*wvelth_eddy - beta_w.*wvelslt_eddy;
  clear('uvelth_eddy','vvelth_eddy','wvelth_eddy','uvelslt_eddy','vvelslt_eddy','wvelslt_eddy');
  
  %%% Interpolate eddy 'buoyancy' fluxes to cell corners
  uvelbuoy_eddy_uw = NaN*ones(Nx,Ny,Nr+1);
  uvelbuoy_eddy_uw(:,:,2:Nr) = 0.5*(uvelbuoy_eddy_u(:,:,1:Nr-1)+uvelbuoy_eddy_u(:,:,2:Nr));
  vvelbuoy_eddy_vw = NaN*ones(Nx,Ny,Nr+1);
  vvelbuoy_eddy_vw(:,:,2:Nr) = 0.5*(vvelbuoy_eddy_v(:,:,1:Nr-1)+vvelbuoy_eddy_v(:,:,2:Nr));
  wvelbuoy_eddy_uw = 0.5*(wvelbuoy_eddy_w([1:Nx],:,:)+wvelbuoy_eddy_w([Nx 1:Nx-1],:,:));
  wvelbuoy_eddy_vw = 0.5*(wvelbuoy_eddy_w(:,[1:Ny],:)+wvelbuoy_eddy_w(:,[Ny 1:Ny-1],:));
  clear('uvelbuoy_eddy_u','vvelbuoy_eddy_v','wvelbuoy_eddy_w');
  
  %%% Compute mean temperature and salinity gradients on cell faces
  DXC_3D = repmat(DXC,[1 1 Nr]);
  DYC_3D = repmat(DYC,[1 1 Nr]);
  DRC_3D = repmat(reshape(DRC,[1 1 Nr+1]),[Nx Ny 1]);
  dsalt_dx_u = (salt([1:Nx],:,:)-salt([Nx 1:Nx-1],:,:)) ./ DXC_3D;
  dsalt_dy_v = (salt(:,[1:Ny],:)-salt(:,[Ny 1:Ny-1],:)) ./ DYC_3D;  
  dsalt_dz_w = NaN*ones(Nx,Ny,Nr+1);
  dsalt_dz_w(:,:,2:Nr) = -diff(salt,1,3) ./ DRC_3D(:,:,2:Nr);
  dtheta_dx_u = (theta([1:Nx],:,:)-theta([Nx 1:Nx-1],:,:)) ./ DXC_3D;
  dtheta_dy_v = (theta(:,[1:Ny],:)-theta(:,[Ny 1:Ny-1],:)) ./ DYC_3D;  
  dtheta_dz_w = NaN*ones(Nx,Ny,Nr+1);
  dtheta_dz_w(:,:,2:Nr) = -diff(theta,1,3) ./ DRC_3D(:,:,2:Nr);
  clear('DXC_3D','DYC_3D','DRC_3D');
  
  %%% Compute mean 'buoyancy' gradients on cell faces
  dbuoy_dx_u = alpha_u.*dtheta_dx_u - beta_u.*dsalt_dx_u;
  dbuoy_dy_v = alpha_v.*dtheta_dy_v - beta_v.*dsalt_dy_v;
  dbuoy_dz_w = alpha_w.*dtheta_dz_w - beta_w.*dsalt_dz_w;
  clear('dtheta_dx_u','dtheta_dy_v','dtheta_dz_w','dsalt_dx_u','dsalt_dy_v','dsalt_dz_w');
  clear('alpha_u','alpha_v','alpha_w','beta_u','beta_v','beta_w');
  
  %%% Interpolate mean 'buoyancy' gradients to cell corners
  dbuoy_dx_uw = NaN*ones(Nx,Ny,Nr+1);
  dbuoy_dx_uw(:,:,2:Nr) = 0.5*(dbuoy_dx_u(:,:,1:Nr-1)+dbuoy_dx_u(:,:,2:Nr));
  dbuoy_dy_vw = NaN*ones(Nx,Ny,Nr+1);
  dbuoy_dy_vw(:,:,2:Nr) = 0.5*(dbuoy_dy_v(:,:,1:Nr-1)+dbuoy_dy_v(:,:,2:Nr));
  dbuoy_dz_uw = 0.5*(dbuoy_dz_w([1:Nx],:,:)+dbuoy_dz_w([Nx 1:Nx-1],:,:));
  dbuoy_dz_vw = 0.5*(dbuoy_dz_w(:,[1:Ny],:)+dbuoy_dz_w(:,[Ny 1:Ny-1],:));
  clear('dbuoy_dx_u','dbuoy_dy_v','dbuoy_dz_w');
  
  %%% Compute components of TRM streamfunction
  PsiX = (uvelbuoy_eddy_uw .* dbuoy_dz_uw - wvelbuoy_eddy_uw .* dbuoy_dx_uw) ./ (dbuoy_dz_ref.^2 + dbuoy_dx_uw.^2 + dbuoy_dz_uw.^2);
  PsiY = (vvelbuoy_eddy_vw .* dbuoy_dz_vw - wvelbuoy_eddy_vw .* dbuoy_dy_vw) ./ (dbuoy_dz_ref.^2 + dbuoy_dy_vw.^2 + dbuoy_dz_vw.^2);
%   PsiX = (uvelbuoy_eddy_uw) ./ sqrt(dbuoy_dz_ref.^2 + dbuoy_dz_uw.^2);
%   PsiY = (vvelbuoy_eddy_vw) ./ sqrt(dbuoy_dz_ref.^2 + dbuoy_dz_vw.^2);
  clear('uvelbuoy_eddy_uw','vvelbuoy_eddy_vw','wvelbuoy_eddy_uw','wvelbuoy_eddy_vw', ...
    'dbuoy_dz_uw','dbuoy_dz_vw','dbuoy_dx_uw','dbuoy_dy_vw');
  
  %%% NaNs should correspond to land points
  PsiX(isnan(PsiX)) = 0;
  PsiY(isnan(PsiY)) = 0;
  
  %%% Compute eddy velocities from streamfunction
  DXG_3D = repmat(DXG,[1 1 Nr]);
  DYG_3D = repmat(DYG,[1 1 Nr]);
  RAC_3D = repmat(RAC,[1 1 Nr]); 
  DRF_3D = repmat(reshape(DRF,[1 1 Nr]),[Nx Ny 1]);
  u_eddy = diff(PsiX,1,3) ./ (DRF_3D .* hFacW); %%% N.B. this is -dPsiX/dz
  u_eddy(hFacW==0) = 0;
  v_eddy = diff(PsiY,1,3) ./ (DRF_3D .* hFacS);
  v_eddy(hFacS==0) = 0;
  w_eddy = ((PsiX([2:Nx 1],:,1:Nr) - PsiX(1:Nx,:,1:Nr)) .* DYG_3D + (PsiY(:,[2:Ny 1],1:Nr) - PsiY(:,1:Ny,1:Nr)) .* DXG_3D) ./ RAC_3D;
  clear('PsiX','PsiY','DXG_3D','DYG_3D','RAC_3D','DRF_3D');

end