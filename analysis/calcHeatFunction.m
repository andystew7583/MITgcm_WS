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
% tmax = 1.1;
tmax = 7.05;
loadexp;

%%% Reference surface freezing temperature
theta0 = -1.9;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use grounding line coordinate
gl_coord = true;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to use depth-averaged temperature as the coordinate system
use_meanT = false;

%%% Set true to decompose eddy fluxes
calc_eddy_decomp = false;

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

psiT_tot = zeros(Neta,Nr+1,Ntime);
psiT_mean = zeros(Neta,Nr+1,Ntime);
psiT_pos_mean = zeros(Neta,Nr+1,Ntime);
psiT_neg_mean = zeros(Neta,Nr+1,Ntime);
psiT_eddy = zeros(Neta,Nr+1,Ntime);
psiT_eddy_adv = zeros(Neta,Nr+1,Ntime);
psiT_eddy_stir = zeros(Neta,Nr+1,Ntime);
psiS_tot = zeros(Neta,Nr+1,Ntime);
psiS_mean = zeros(Neta,Nr+1,Ntime);
psiS_eddy = zeros(Neta,Nr+1,Ntime);
psiS_eddy_adv = zeros(Neta,Nr+1,Ntime);
psiS_eddy_stir = zeros(Neta,Nr+1,Ntime);
psi_tot = zeros(Neta,Nr+1,Ntime);
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
  
  %%% Compute mean and eddy fluxes in quasi-latitude coordinates
  [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
    uvel,vvel,theta,uvelth,vvelth, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psiT_tot(:,:,n) = eflux_tot;
  psiT_mean(:,:,n) = eflux_mean;
  psiT_eddy(:,:,n) = eflux_eddy;

  [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
    uvel,vvel,salt,uvelslt,vvelslt, ...
    Nx,Ny,Nr,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  psiS_tot(:,:,n) = eflux_tot;
  psiS_mean(:,:,n) = eflux_mean;
  psiS_eddy(:,:,n) = eflux_eddy;

  if (calc_eddy_decomp)
    
    %%% Load heat and salt fluxes         
    wvel  = rdmdsWrapper(fullfile(exppath,'/results/WVEL'),itersToRead(n));  
    wvelth = rdmdsWrapper(fullfile(exppath,'/results/WVELTH'),itersToRead(n));
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
    
  else
    
    %%% Clear memory
    clear('salt','uvelth','vvelth','uvelslt','vvelslt');
    
  end

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