%%%
%%% calcHeatFunctionEddyDecomp.m
%%% 
%%% Calculates total heat (and salt) functions in quasi-latitude or
%%% streamline coordinates. This function decomposes the "eddy" component
%%% of the streamfunction into advective and stirring components by
%%% computing the temporal residual mean velocity.
%%%





%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% tmin = 19.05;
% tmax = 27.05;
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% tmin = 10.05;
% tmax = 18.05;
% expname = 'hires_seq_onetwelfth_RTOPO2';
% tmin = 1.05;
% tmax = 9.05;
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
tmin = 1.05;
tmax = 7.05;
loadexp;

%%% Reference surface freezing temperature
theta0 = -1.9;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to use depth-averaged temperature as the coordinate system
use_meanT = false;

%%% Set true to use grounding line coordinate
gl_coord = true;

%%% Index of the upper grid cell face dividing the upper and lower portions
%%% of the water column
if (gl_coord)
  zidx_icefront = 15;
else
  zidx_icefront = 25;
end

%%% Define coordinate system for integrating to compute heatfunction
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
  
  if (use_meanT)

    %%% Load time-mean temperature
    outfname = [expname,'_TSfluxes.mat'];
    load(fullfile('./products',outfname),'theta_tavg');
    
    %%% Coordinate
    ETA = sum(theta_tavg.*DRF.*hFacC,3)./sum(DRF.*hFacC,3);
    eta = -2.6:0.025:0;
    Neta = length(eta);

  else

    ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity,gl_coord);
    eta = -9:.05:11;
    Neta = length(eta);

  end

end

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);









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
%%% HEAT/SALT FLUX CALCULATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psiT_eddy_adv = zeros(Neta,Nr+1,Ntime);
psiT_eddy_stir = zeros(Neta,Nr+1,Ntime);
psiS_eddy_adv = zeros(Neta,Nr+1,Ntime);
psiS_eddy_stir = zeros(Neta,Nr+1,Ntime);
psi_eddy = zeros(Neta,Nr+1,Ntime);
w_eddy_flux = zeros(Nx,Ny,Ntime);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))]      
   
  %%% Reference stratification to regularize the TRM where stratification is weak
  %%% N.B. This differs from the actual stratification N^2 by a factor of g
  dbuoy_dz_ref = 1e-8;
%   dbuoy_dz_ref = 1e-9;

  %%% Grid sizes
  Nx = size(hFacC,1);
  Ny = size(hFacC,2);
  Nr = size(hFacC,3);
  
  %%% Load temperature and salinity fields
  theta  = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));
  salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT'),itersToRead(n));

  %%% Remove dry grid cells. Should ensure that streamfunction only gets
  %%% calculated at points surrouned by wet cells
  salt(hFacC==0) = NaN;
  theta(hFacC==0) = NaN;

  %%% ZONAL FLUXES %%%

  %%% Calculate midpoint salinity and temperature
  salt_u = 0.5*(salt([1:Nx],:,:)+salt([Nx 1:Nx-1],:,:));
  theta_u = 0.5*(theta([1:Nx],:,:)+theta([Nx 1:Nx-1],:,:));

  %%% Load u-vel and u-fluxes
  uvelslt_eddy  = rdmdsWrapper(fullfile(exppath,'/results/UVELSLT'),itersToRead(n));
  uvelth_eddy  = rdmdsWrapper(fullfile(exppath,'/results/UVELTH'),itersToRead(n));  
  uvel  = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),itersToRead(n));

  %%% Compute u-eddy fluxes
  uvelslt_eddy = uvelslt_eddy - uvel .* salt_u;
  uvelth_eddy = uvelth_eddy - uvel .* theta_u;
  clear('uvel');

  %%% Compute thermal expansion and haline contraction coefficients
  press_c = -rhoConst*gravity*repmat(RC,[Nx Ny 1])/1e4; %%% N.B. Units in dbar
  [alpha_u,beta_u] = calcAlphaBeta(salt_u,theta_u,press_c);
  clear('salt_u','theta_u','press_c');

  %%% Compute eddy buoyancy flux
  uvelbuoy_eddy_u = alpha_u.*uvelth_eddy - beta_u.*uvelslt_eddy;
  clear('uvelth_eddy','uvelslt_eddy');

  %%% Compute buoyancy gradient
  dsalt_dx_u = (salt([1:Nx],:,:)-salt([Nx 1:Nx-1],:,:)) ./ repmat(DXC,[1 1 Nr]);
  dtheta_dx_u = (theta([1:Nx],:,:)-theta([Nx 1:Nx-1],:,:)) ./ repmat(DXC,[1 1 Nr]);
  dbuoy_dx_u = alpha_u.*dtheta_dx_u - beta_u.*dsalt_dx_u;
  clear('dsalt_dx_u','dtheta_dx_u','alpha_u','beta_u');

  %%% MERIDIONAL FLUXES %%%

  %%% Calculate midpoint salinity and temperature
  salt_v = 0.5*(salt(:,[1:Ny],:)+salt(:,[Ny 1:Ny-1],:));
  theta_v = 0.5*(theta(:,[1:Ny],:)+theta(:,[Ny 1:Ny-1],:));

  %%% Load v-vel and v-fluxes
  vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),itersToRead(n));
  vvelth_eddy  = rdmdsWrapper(fullfile(exppath,'/results/VVELTH'),itersToRead(n));
  vvelslt_eddy  = rdmdsWrapper(fullfile(exppath,'/results/VVELSLT'),itersToRead(n));

  %%% Compute v-eddy fluxes
  vvelslt_eddy = vvelslt_eddy - vvel .* salt_v;
  vvelth_eddy = vvelth_eddy - vvel .* theta_v;
  clear('vvel');

  %%% Compute thermal expansion and haline contraction coefficients
  press_c = -rhoConst*gravity*repmat(RC,[Nx Ny 1])/1e4; %%% N.B. Units in dbar
  [alpha_v,beta_v] = calcAlphaBeta(salt_v,theta_v,press_c);
  clear('salt_v','theta_v','press_c');

  %%% Compute eddy buoyancy flux
  vvelbuoy_eddy_v = alpha_v.*vvelth_eddy - beta_v.*vvelslt_eddy;
  clear('vvelth_eddy','vvelslt_eddy');

  %%% Compute buoyancy gradient
  dsalt_dy_v = (salt(:,[1:Ny],:)-salt(:,[Ny 1:Ny-1],:)) ./ repmat(DYC,[1 1 Nr]);  
  dtheta_dy_v = (theta(:,[1:Ny],:)-theta(:,[Ny 1:Ny-1],:)) ./ repmat(DYC,[1 1 Nr]);
  dbuoy_dy_v = alpha_v.*dtheta_dy_v - beta_v.*dsalt_dy_v;
  clear('dsalt_dy_v','dtheta_dy_v','alpha_v','beta_v');

  %%% VERTICAL COMPONENT %%%

  %%% Calculate midpoint salinity and temperature
  salt_w = NaN*ones(Nx,Ny,Nr+1);
  salt_w(:,:,2:Nr) = 0.5*(salt(:,:,1:Nr-1)+salt(:,:,2:Nr));
  theta_w = NaN*ones(Nx,Ny,Nr+1);
  theta_w(:,:,2:Nr) = 0.5*(theta(:,:,1:Nr-1)+theta(:,:,2:Nr));
 
  %%% Load w-vel and w-fluxes
  %%% Compute w-eddy fluxes
  wvel  = rdmdsWrapper(fullfile(exppath,'/results/WVEL'),itersToRead(n));    
  wvelslt = rdmdsWrapper(fullfile(exppath,'/results/WVELSLT'),itersToRead(n));  
  wvelslt_eddy = NaN*ones(Nx,Ny,Nr+1);
  wvelslt_eddy(:,:,1:Nr) = wvelslt - wvel .* salt_w(:,:,1:Nr);
  clear('wvelslt')
  wvelth = rdmdsWrapper(fullfile(exppath,'/results/WVELTH'),itersToRead(n));
  wvelth_eddy = NaN*ones(Nx,Ny,Nr+1);
  wvelth_eddy(:,:,1:Nr) = wvelth - wvel .* theta_w(:,:,1:Nr);
  clear('wvel','wvelth');

  %%% Compute thermal expansion and haline contraction coefficients
  press_w = -rhoConst*gravity*repmat(RF,[Nx Ny 1])/1e4;
  [alpha_w,beta_w] = calcAlphaBeta(salt_w,theta_w,press_w);
  clear('press_w','theta_w','salt_w');

  %%% Compute eddy buoyancy flux
  wvelbuoy_eddy_w = alpha_w.*wvelth_eddy - beta_w.*wvelslt_eddy;
  clear('wvelth_eddy','wvelslt_eddy');

  %%% Compute buoyancy gradient
  DRC_3D = repmat(reshape(DRC,[1 1 Nr+1]),[Nx Ny 1]);
  dsalt_dz_w = NaN*ones(Nx,Ny,Nr+1);
  dsalt_dz_w(:,:,2:Nr) = -diff(salt,1,3) ./ DRC_3D(:,:,2:Nr);    
  dtheta_dz_w = NaN*ones(Nx,Ny,Nr+1);
  dtheta_dz_w(:,:,2:Nr) = -diff(theta,1,3) ./ DRC_3D(:,:,2:Nr);
  clear('theta','salt');
  dbuoy_dz_w = alpha_w.*dtheta_dz_w - beta_w.*dsalt_dz_w;
  clear('DRC_3D','dsalt_dz_w','dtheta_dz_w','alpha_w','beta_w');

  %%% STREAMFUNCTION %%%

  %%% Interpolate eddy 'buoyancy' fluxes to cell corners
  uvelbuoy_eddy_uw = NaN*ones(Nx,Ny,Nr+1);
  uvelbuoy_eddy_uw(:,:,2:Nr) = 0.5*(uvelbuoy_eddy_u(:,:,1:Nr-1)+uvelbuoy_eddy_u(:,:,2:Nr));
  clear('uvelbuoy_eddy_u');
  vvelbuoy_eddy_vw = NaN*ones(Nx,Ny,Nr+1);
  vvelbuoy_eddy_vw(:,:,2:Nr) = 0.5*(vvelbuoy_eddy_v(:,:,1:Nr-1)+vvelbuoy_eddy_v(:,:,2:Nr));
  clear('vvelbuoy_eddy_v');
  wvelbuoy_eddy_uw = 0.5*(wvelbuoy_eddy_w([1:Nx],:,:)+wvelbuoy_eddy_w([Nx 1:Nx-1],:,:));
  wvelbuoy_eddy_vw = 0.5*(wvelbuoy_eddy_w(:,[1:Ny],:)+wvelbuoy_eddy_w(:,[Ny 1:Ny-1],:));
  clear('wvelbuoy_eddy_w');
    
  %%% Interpolate mean 'buoyancy' gradients to cell corners
  dbuoy_dx_uw = NaN*ones(Nx,Ny,Nr+1);
  dbuoy_dx_uw(:,:,2:Nr) = 0.5*(dbuoy_dx_u(:,:,1:Nr-1)+dbuoy_dx_u(:,:,2:Nr));
  clear('dbuoy_dx_u');
  dbuoy_dy_vw = NaN*ones(Nx,Ny,Nr+1);
  dbuoy_dy_vw(:,:,2:Nr) = 0.5*(dbuoy_dy_v(:,:,1:Nr-1)+dbuoy_dy_v(:,:,2:Nr));
  clear('dbuoy_dy_v');
  dbuoy_dz_uw = 0.5*(dbuoy_dz_w([1:Nx],:,:)+dbuoy_dz_w([Nx 1:Nx-1],:,:));
  dbuoy_dz_vw = 0.5*(dbuoy_dz_w(:,[1:Ny],:)+dbuoy_dz_w(:,[Ny 1:Ny-1],:));
  clear('dbuoy_dz_w');
  
  %%% Compute components of TRM streamfunction
  PsiX = (uvelbuoy_eddy_uw .* dbuoy_dz_uw - wvelbuoy_eddy_uw .* dbuoy_dx_uw) ./ (dbuoy_dz_ref.^2 + dbuoy_dx_uw.^2 + dbuoy_dz_uw.^2);
  PsiY = (vvelbuoy_eddy_vw .* dbuoy_dz_vw - wvelbuoy_eddy_vw .* dbuoy_dy_vw) ./ (dbuoy_dz_ref.^2 + dbuoy_dy_vw.^2 + dbuoy_dz_vw.^2);
  clear('uvelbuoy_eddy_uw','vvelbuoy_eddy_vw','wvelbuoy_eddy_uw','wvelbuoy_eddy_vw', ...
    'dbuoy_dz_uw','dbuoy_dz_vw','dbuoy_dx_uw','dbuoy_dy_vw');
  
  %%% NaNs should correspond to land points
  PsiX(isnan(PsiX)) = 0;
  PsiY(isnan(PsiY)) = 0;
  
  %%% Compute eddy velocities from streamfunction
  RAC_3D = repmat(RAC,[1 1 Nr]); 
  DXG_3D = repmat(DXG,[1 1 Nr]);
  DYG_3D = repmat(DYG,[1 1 Nr]);
  w_eddy = ((PsiX([2:Nx 1],:,1:Nr) - PsiX(1:Nx,:,1:Nr)) .* DYG_3D + (PsiY(:,[2:Ny 1],1:Nr) - PsiY(:,1:Ny,1:Nr)) .* DXG_3D) ./ RAC_3D;
  clear('DXG_3D','DYG_3D','RAC_3D');
  DRF_3D = repmat(reshape(DRF,[1 1 Nr]),[Nx Ny 1]);
  u_eddy = diff(PsiX,1,3) ./ (DRF_3D .* hFacW); %%% N.B. this is -dPsiX/dz
  u_eddy(hFacW==0) = 0;
  v_eddy = diff(PsiY,1,3) ./ (DRF_3D .* hFacS);  
  v_eddy(hFacS==0) = 0;
  clear('PsiX','PsiY','DRF_3D');

  %%% Store eddy velocity at the interface between upper/lower shelf
  %%% waters
  w_eddy_flux(:,:,n) = w_eddy(:,:,zidx_icefront);
  clear('w_eddy');

  %%% 3D grid spacing matrices
  DXG_3D = repmat(DXG,[1 1 Nr]);
  DYG_3D = repmat(DYG,[1 1 Nr]);
  DRF_3D = repmat(DRF,[Nx Ny 1]);

  %%% Load temperature and salinity fields
  theta  = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));
  salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT'),itersToRead(n));
  
  %%% Decompose eddy heat and salt fluxes into advective and stirring components

  [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
    u_eddy,v_eddy,theta,0*u_eddy,0*v_eddy, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psiT_eddy_adv(:,:,n) = eflux_mean;  
  
  [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
    u_eddy,v_eddy,salt,0*u_eddy,0*v_eddy, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psiS_eddy_adv(:,:,n) = eflux_mean;  

  %%% Eddy-induced overturning in depth space
  eflux = calcQuasiLatFluxes (...
    u_eddy,v_eddy, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psi_eddy(:,:,n) = eflux;

  
  %%% Clear memory
  clear('u_eddy','v_eddy','w_eddy','salt','theta');

    
  
end

%%% Store computed data for later
outfname = [expname,'_HeatFunctionEddyDecomp'];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else 
  if (use_meanT)
    outfname = [outfname,'_meanT'];
  else 
    if (deform_cavity)
      outfname = [outfname,'_deform'];
    elseif (gl_coord)
      outfname = [outfname,'_GLcoord'];
    end
  end
end
outfname = [outfname,'.mat'];
save(fullfile('products',outfname), ...
  'psiT_eddy_adv',...
  'psiS_eddy_adv',...
  'w_eddy_flux', ...
  'psi_eddy', ...
  '-v7.3');





%%%
%%% Convenience function to compute mean and eddy fluxes
%%%
function [trflux_tot,trflux_mean,trflux_eddy] = calcMeanEddyFluxes (...
  uvel,vvel,tracer,uveltr,vveltr, ...
  Nx,Ny,Nr,Neta, ...  
  DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta)

  %%% Calculate midpoint tracer
  tracer_u = 0.5*(tracer([1:Nx],:,:)+tracer([Nx 1:Nx-1],:,:));
  tracer_v = 0.5*(tracer(:,[1:Ny],:)+tracer(:,[Ny 1:Ny-1],:));
  
  %%% Compute eddy heat and salt fluxes on cell faces
  uveltr_mean = uvel.*tracer_u;
  vveltr_mean = vvel.*tracer_v;
  uveltr_eddy = uveltr - uveltr_mean;
  vveltr_eddy = vveltr - vveltr_mean;
  
  %%% Calculate fluxes in quasi-latitude coordinates
  trflux_tot = calcQuasiLatFluxes (...
    uveltr,vveltr, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  trflux_mean = calcQuasiLatFluxes (...
    uveltr_mean,vveltr_mean, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  trflux_eddy = calcQuasiLatFluxes (...
    uveltr_eddy,vveltr_eddy, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);

end








%%%
%%% Convenience function to comute fluxes in quasi-latitude space
%%%
function eflux = calcQuasiLatFluxes (...
  uflux,vflux, ...
  Nx,Ny,Nr,Neta, ...  
  DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta)

  %%% Integrate fluxes verticall and horizontally over each cell face
  uflux_yzint = zeros(Nx,Ny,Nr+1);
  vflux_xzint = zeros(Nx,Ny,Nr+1);
  uflux_yzint(:,:,1:Nr) = cumsum(uflux .* DYG_3D .* DRF_3D .* hFacW,3,'reverse');
  vflux_xzint(:,:,1:Nr) = cumsum(vflux .* DXG_3D .* DRF_3D .* hFacS,3,'reverse');

  %%% Compute horizontal divergence of isopycnal fluxes
  fluxdiv = zeros(Nx,Ny,Nr+1);
  fluxdiv(1:Nx-1,1:Ny-1,:) = uflux_yzint(2:Nx,1:Ny-1,:) ...
                          - uflux_yzint(1:Nx-1,1:Ny-1,:) ...
                          + vflux_xzint(1:Nx-1,2:Ny,:) ...
                          - vflux_xzint(1:Nx-1,1:Ny-1,:);
                       
  %%% Integrate flux divergence across lines of constant eta 
  eflux = zeros(Neta,Nr+1,1);
  for m = 1:Neta
    msk = repmat(ETA<eta(m),[1 1 Nr+1]);
    eflux(m,:) = squeeze(sum(sum(fluxdiv.*msk,1),2));
  end  

end




