%%%
%%% calcHeatFunction.m
%%% 
%%% Calculates total heat (and salt) functions in quasi-latitude or streamline coordinates.
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

%%% Index of the upper grid cell face dividing the upper and lower portions
%%% of the water column
zidx_icefront = 25;

%%% Reference surface freezing temperature
theta0 = -1.9;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to decompose eddy fluxes
calc_eddy_decomp = true;

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

%%% 3D grid spacing matrices
DXG_3D = repmat(DXG,[1 1 Nr]);
DYG_3D = repmat(DYG,[1 1 Nr]);
DRF_3D = repmat(DRF,[Nx Ny 1]);







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

psiT_pos_mean = zeros(Neta,Nr+1,Ntime);
psiT_neg_mean = zeros(Neta,Nr+1,Ntime);
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
  
  %%% Read velocity field
  uvel  = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),itersToRead(n));
  vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),itersToRead(n));  
  theta  = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));
  salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT'),itersToRead(n));
  uvelth  = rdmdsWrapper(fullfile(exppath,'/results/UVELTH'),itersToRead(n));
  vvelth  = rdmdsWrapper(fullfile(exppath,'/results/VVELTH'),itersToRead(n));
  uvelslt  = rdmdsWrapper(fullfile(exppath,'/results/UVELSLT'),itersToRead(n));
  vvelslt  = rdmdsWrapper(fullfile(exppath,'/results/VVELSLT'),itersToRead(n));
  if (isempty(uvel) || isempty(vvel) || isempty(theta) || isempty(salt) ...
      || isempty(uvelth) || isempty(vvelth) || isempty(uvelslt) || isempty(vvelslt) )
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end

  %%% Load heat and salt fluxes         
  wvel  = rdmdsWrapper(fullfile(exppath,'/results/WVEL'),itersToRead(n));  
  wvelth = rdmdsWrapper(fullfile(exppath,'/results/WVELTH'),itersToRead(n));
  wvelslt = rdmdsWrapper(fullfile(exppath,'/results/WVELSLT'),itersToRead(n));

   
  %%% Reference stratification to regularize the TRM where stratification is weak
  %%% N.B. This differs from the actual stratification N^2 by a factor of g
  dbuoy_dz_ref = 1e-8;
%   dbuoy_dz_ref = 1e-9;

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


    
  %%% Clear memory
  clear('wvel','uvelth','vvelth','wvelth','uvelslt','vvelslt','wvelslt');
  
  %%% Compute mean and eddy fluxes in quasi-latitude coordinates
  [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
    u_eddy,v_eddy,theta,0*u_eddy,0*v_eddy, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psiT_eddy_adv(:,:,n) = eflux_mean;
  psiT_eddy_stir(:,:,n) = psiT_eddy(:,:,n) - psiT_eddy_adv(:,:,n);
  [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
    u_eddy,v_eddy,salt,0*u_eddy,0*v_eddy, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psiS_eddy_adv(:,:,n) = eflux_mean;
  psiS_eddy_stir(:,:,n) = psiS_eddy(:,:,n) - psiS_eddy_adv(:,:,n);

  %%% Store eddy velocity at the interface between upper/lower shelf
  %%% waters
  w_eddy_flux(:,:,n) = w_eddy(:,:,zidx_icefront);

  %%% Eulerian-mean overturning
  eflux = calcQuasiLatFluxes (...
    u_eddy,v_eddy, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psi_eddy(:,:,n) = eflux;

  
  %%% Clear memory
  clear('u_eddy','v_eddy','w_eddy','salt');
    

  %%% Positive and negative temperature anomalies
  theta_pos = theta - theta0;
  theta_neg = min(theta_pos,0);
  theta_pos = max(theta_pos,0);  

  %%% Mean flux of positive temperature anomaly
  [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
    uvel,vvel,theta_pos,uvel,vvel, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psiT_pos_mean(:,:,n) = eflux_mean;
  
  %%% Mean flux of negative temperature anomaly
  [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
    uvel,vvel,theta_neg,uvel,vvel, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psiT_neg_mean(:,:,n) = eflux_mean;
  
  %%% Eulerian-mean overturning
  eflux = calcQuasiLatFluxes (...
    uvel,vvel, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psi_tot(:,:,n) = eflux;

  %%% Clear memory
  clear('uvel','vvel','theta','theta_pos','theta_neg');
    
  
end

%%% Read time-averaged variables
uvel_tavg = readIters(exppath,'UVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
vvel_tavg = readIters(exppath,'VVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
theta_tavg = readIters(exppath,'THETA',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
salt_tavg = readIters(exppath,'SALT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
uvelth_tavg = readIters(exppath,'UVELTH',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
vvelth_tavg = readIters(exppath,'VVELTH',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
uvelslt_tavg = readIters(exppath,'UVELSLT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
vvelslt_tavg = readIters(exppath,'VVELSLT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);

%%% Compute multi-annual mean mean and eddy fluxes in quasi-latitude coordinates
%%% N.B. Total multi-annual mean flux = flux_stand + flux_fluc + flux_eddy
[eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
  uvel_tavg,vvel_tavg,theta_tavg,uvelth_tavg,vvelth_tavg, ...
  Nx,Ny,Nr,Neta, ...  
  DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
psiT_stand = eflux_mean;
psiT_fluc = eflux_eddy - mean(psiT_eddy,3);
[eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
  uvel_tavg,vvel_tavg,salt_tavg,uvelslt_tavg,vvelslt_tavg, ...
  Nx,Ny,Nr,Neta, ...  
  DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
psiS_stand = eflux_mean;
psiS_fluc = eflux_eddy - mean(psiS_eddy,2);


%%% Store computed data for later
outfname = [expname,'_HeatFunction'];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else
  if (deform_cavity)
    outfname = [outfname,'_deform'];
  end
end
outfname = [outfname,'.mat'];
save(fullfile('products',outfname), ...
  'eta','ETA','times', ...
  'psiT_tot','psiT_mean','psiT_stand','psiT_fluc','psiT_eddy',...
  'psiS_tot','psiS_mean','psiS_stand','psiS_fluc','psiS_eddy',...
  'psiT_pos_mean','psiT_neg_mean','psi_tot','psi_eddy', ...
  'uvel_tavg','vvel_tavg','theta_tavg','salt_tavg', ... 
  'uvelth_tavg','vvelth_tavg','uvelslt_tavg','vvelslt_tavg','-v7.3');
if (calc_eddy_decomp)
  save(fullfile('products',outfname), ...
  'psiT_eddy_adv','psiT_eddy_stir',...
  'psiS_eddy_adv','psiS_eddy_stir',...
  'w_eddy_flux', ...
  '-append','-v7.3');
end
clear('uvel_tavg','vvel_tavg','theta_tavg','salt_tavg','uvelth_tavg','vvelth_tavg','uvelslt_tavg','vvelslt_tavg');





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










% uvel_tavg = readIters(exppath,'UVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% vvel_tavg = readIters(exppath,'VVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% theta_tavg = readIters(exppath,'THETA',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% salt_tavg = readIters(exppath,'SALT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% uvelth_tavg = readIters(exppath,'UVELTH',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% vvelth_tavg = readIters(exppath,'VVELTH',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% uvelslt_tavg = readIters(exppath,'UVELSLT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% vvelslt_tavg = readIters(exppath,'VVELSLT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% wvelth_tavg = readIters(exppath,'WVELTH',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% wvelslt_tavg = readIters(exppath,'WVELSLT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% wvel_tavg = readIters(exppath,'WVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
% %%% Compute eddy-induced velocities
%     [u_eddy,v_eddy,w_eddy] = calcTRMvelocity (...
%       uvel_tavg,vvel_tavg,wvel_tavg,theta_tavg,salt_tavg, ...
%       uvelth_tavg,vvelth_tavg,wvelth_tavg, ...
%       uvelslt_tavg,vvelslt_tavg,wvelslt_tavg, ...
%       hFacC,hFacW,hFacS, ...
%       DXG,DYG,RAC,DXC,DYC, ...
%       DRF,DRC,RC,RF,...
%       rhoConst,gravity);
% w_eddy_flux = w_eddy(:,:,zidx_icefront);
% clear('uvelth_tavg','vvelth_tavg','uvelslt_tavg','vvelslt_tavg','wvelth_tavg','wvelslt_tavg');
% %%% Positive and negative temperature anomalies
%   theta_pos = theta_tavg - theta0;
%   theta_neg = min(theta_pos,0);
%   theta_pos = max(theta_pos,0);
% %%% Mean flux of positive temperature anomaly
%   [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
%     uvel_tavg,vvel_tavg,theta_pos,uvel_tavg,vvel_tavg, ...
%     Nx,Ny,Nr,Neta, ...  
%     DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
%   psiT_pos_mean = eflux_mean;
% %%% Mean flux of negative temperature anomaly
%   [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
%     uvel_tavg,vvel_tavg,theta_neg,uvel_tavg,vvel_tavg, ...
%     Nx,Ny,Nr,Neta, ...  
%     DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
%   psiT_neg_mean = eflux_mean;
% save('./products/tmp.mat','w_eddy_flux','psiT_pos_mean','psiT_neg_mean');














% 
% setExpname;
% expname = 'hires_seq_onetwentyfourth_RTOPO2';
% loadexp;
% addpath ../utils/matlab/
% 
% iter = 46203;
% 
% Eta = rdmds(fullfile(basedir,expname,'results','Eta'),iter);
% figure(1);pcolor(XC,YC,Eta);shading interp;colorbar
% 
% Heff = rdmds(fullfile(basedir,expname,'results','HEFF'),iter);
% Heff(hFacC(:,:,1)==0) = NaN;
% figure(2);pcolor(XC,YC,Heff);shading interp;colorbar
% 
% 
% KPPdiffKzT = rdmds(fullfile(basedir,expname,'results','KPPdiffKzT'),iter);
% KPPdiffKzT(hFacC==0) = NaN;
% [i,j] = find(KPPdiffKzT==max(KPPdiffKzT(:)));
% [j,k] = find(squeeze(KPPdiffKzT(i,:,:))==max(KPPdiffKzT(:)));
% figure(3);pcolor(XC,YC,log10(KPPdiffKzT(:,:,k)));shading flat;colorbar;caxis([-2 3]);
% [ZZ,YY] = meshgrid(RC,YC(i,:));
% figure(4);pcolor(YY,ZZ,log10(squeeze(KPPdiffKzT(i,:,:))));shading flat;colorbar;caxis([-2 3]);
% 
% 
% Qnet = rdmds(fullfile(basedir,expname,'results','Qnet'),iter);
% Qnet(hFacC(:,:,1)==0) = NaN;
% figure(5);pcolor(XC,YC,Qnet);shading flat;colorbar
% [i,j]=find(abs(Qnet)==max(abs(Qnet(:))));
% 
% Vwind = rdmds(fullfile(basedir,expname,'results','VWIND'),iter);
% Vwind(hFacC(:,:,1)==0) = NaN;
% figure(6);pcolor(XC,YC,Vwind);shading interp;colorbar
% 
% Uwind = rdmds(fullfile(basedir,expname,'results','UWIND'),iter);
% Uwind(hFacC(:,:,1)==0) = NaN;
% figure(7);pcolor(XC,YC,Uwind);shading interp;colorbar
% 
% Area = rdmds(fullfile(basedir,expname,'results','AREA'),iter);
% Area(hFacC(:,:,1)==0) = NaN;
% figure(8);pcolor(XC,YC,Area);shading interp;colorbar
% 
% 
% T = rdmds(fullfile(basedir,expname,'results','T'),iter);
% T(hFacC==0) = NaN;
% figure(9);pcolor(XC,YC,T(:,:,1));shading interp;colorbar
% 
% S = rdmds(fullfile(basedir,expname,'results','S'),iter);
% S(hFacC==0) = NaN;
% figure(10);pcolor(XC,YC,S(:,:,1));shading interp;colorbar
% 
% W = rdmds(fullfile(basedir,expname,'results','W'),iter);
% figure(11);pcolor(XC,YC,W(:,:,25));shading interp;colorbar
% colormap redblue
% [iw,jw] = find(W==max(abs(W(:))));
% [jw,kw] = find(squeeze(W(iw,:,:))==max(abs(W(:))));
% 
% U = rdmds(fullfile(basedir,expname,'results','U'),iter);
% figure(12);pcolor(XC,YC,U(:,:,1));shading interp;colorbar
% colormap redblue
% 
% V = rdmds(fullfile(basedir,expname,'results','V'),iter);
% figure(13);pcolor(XC,YC,V(:,:,1));shading interp;colorbar
% colormap redblue
% 
% EmPmR = rdmds(fullfile(basedir,expname,'results','EmPmR'),iter);
% EmPmR(hFacC(:,:,1)==0) = NaN;
% figure(14);pcolor(XC,YC,EmPmR);shading flat;colorbar
% 
% 
% Uice = rdmds(fullfile(basedir,expname,'results','UICE'),iter);
% figure(15);pcolor(XC,YC,Uice);shading flat;colorbar; colormap redblue; caxis([-4 4]);
% 
% 
% Vice = rdmds(fullfile(basedir,expname,'results','VICE'),iter);
% figure(16);pcolor(XC,YC,Vice);shading flat;colorbar; colormap redblue; caxis([-4 4]);
% 
% 
% figure(17);pcolor(XC,YC,Uice-U(:,:,1));shading flat;colorbar; colormap redblue; caxis([-1 1]);
% figure(18);pcolor(XC,YC,Vice-V(:,:,1));shading flat;colorbar; colormap redblue; caxis([-4 4]);
% 
% EXFqnet = rdmds(fullfile(basedir,expname,'results','EXFqnet'),iter);
% EXFqnet(hFacC(:,:,1)==0) = NaN;
% figure(19);pcolor(XC,YC,EXFqnet);shading flat;colorbar
% 
% oceaTaux = rdmds(fullfile(basedir,expname,'results','oceTAUX'),iter);
% oceaTaux(hFacC(:,:,1)==0) = NaN;
% figure(20);pcolor(XC,YC,oceaTaux);shading flat;colorbar
% 
% SItflux = rdmds(fullfile(basedir,expname,'results','SItflux'),iter);
% SItflux(hFacC(:,:,1)==0) = NaN;
% figure(21);pcolor(XC,YC,SItflux);shading flat;colorbar
% 
% TFLUX = rdmds(fullfile(basedir,expname,'results','TFLUX'),iter);
% TFLUX(hFacC(:,:,1)==0) = NaN;
% figure(22);pcolor(XC,YC,TFLUX);shading flat;colorbar
% 
% SIhl = rdmds(fullfile(basedir,expname,'results','SIhl'),iter);
% SIhl(hFacC(:,:,1)==0) = NaN;
% figure(23);pcolor(XC,YC,SIhl);shading flat;colorbar
% 
% SItices = rdmds(fullfile(basedir,expname,'results','SItices'),iter);
% SItices(hFacC(:,:,1)==0) = NaN;
% figure(24);pcolor(XC,YC,SItices);shading flat;colorbar
% 
% EXFhl = rdmds(fullfile(basedir,expname,'results','EXFhl'),iter);
% EXFhl(hFacC(:,:,1)==0) = NaN;
% figure(25);pcolor(XC,YC,EXFhl);shading flat;colorbar
% 
% EXFhs = rdmds(fullfile(basedir,expname,'results','EXFhs'),iter);
% EXFhs(hFacC(:,:,1)==0) = NaN;
% figure(26);pcolor(XC,YC,EXFhs);shading flat;colorbar
% 
% EXFtaux = rdmds(fullfile(basedir,expname,'results','EXFtaux'),iter);
% EXFtaux(hFacC(:,:,1)==0) = NaN;
% figure(27);pcolor(XC,YC,EXFtaux);shading flat;colorbar
