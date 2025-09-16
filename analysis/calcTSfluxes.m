%%%
%%% calcTSfluxes.m
%%% 
%%% Calculates total heat and salt fluxes in quasi-latitude coordinates.
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
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
% tmin = 1.05;
% tmax = 9.05;
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
tmin = 1.05;
tmax = 7.05;
loadexp;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use grounding line coordinate
gl_coord = true;

%%% Set true to decompose eddy fluxes
calc_eddy_decomp = false;

%%% Define coordinate system for integrating to compute streamfunction
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity,gl_coord);
eta = -9:.05:11;
Neta = length(eta);

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

thflux_tot = zeros(Neta,Ntime);
thflux_mean = zeros(Neta,Ntime);
thflux_eddy = zeros(Neta,Ntime);
thflux_eddy_adv = zeros(Neta,Ntime);
thflux_eddy_stir = zeros(Neta,Ntime);
sltflux_tot = zeros(Neta,Ntime);
sltflux_mean = zeros(Neta,Ntime);
sltflux_eddy = zeros(Neta,Ntime);
sltflux_eddy_adv = zeros(Neta,Ntime);
sltflux_eddy_stir = zeros(Neta,Ntime);
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
    Nx,Ny,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  thflux_tot(:,n) = eflux_tot;
  thflux_mean(:,n) = eflux_mean;
  thflux_eddy(:,n) = eflux_eddy;
  [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
    uvel,vvel,salt,uvelslt,vvelslt, ...
    Nx,Ny,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  sltflux_tot(:,n) = eflux_tot;
  sltflux_mean(:,n) = eflux_mean;
  sltflux_eddy(:,n) = eflux_eddy;

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
    clear('uvel','vvel','wvel','uvelth','vvelth','wvelth','uvelslt','vvelslt','wvelslt','w_eddy');
    
    %%% Compute mean and eddy fluxes in quasi-latitude coordinates
    [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
      u_eddy,v_eddy,theta,0*u_eddy,0*v_eddy, ...
      Nx,Ny,Neta, ...  
      DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
    thflux_eddy_adv(:,n) = eflux_mean;
    thflux_eddy_stir(:,n) = thflux_eddy(:,n) - thflux_eddy_adv(:,n);
    [eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
      u_eddy,v_eddy,salt,0*u_eddy,0*v_eddy, ...
      Nx,Ny,Neta, ...  
      DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
    sltflux_eddy_adv(:,n) = eflux_mean;
    sltflux_eddy_stir(:,n) = sltflux_eddy(:,n) - sltflux_eddy_adv(:,n);
    
    clear('u_eddy','v_eddy','theta','salt');
    
  else
    
    %%% Just clear memory
    clear('uvel','vvel','theta','salt','uvelth','vvelth','uvelslt','vvelslt');
    
  end
  
  
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
tflux_tavg = readIters(exppath,'TFLUX',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,1);
sflux_tavg = readIters(exppath,'SFLUX',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,1);
SHIfwFlx_tavg = readIters(exppath,'SHIfwFlx',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,1);
SHIhtFlx_tavg = readIters(exppath,'SHIhtFlx',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,1);

%%% Compute multi-annual mean mean and eddy fluxes in quasi-latitude coordinates
%%% N.B. Total multi-annual mean flux = flux_stand + flux_fluc + flux_eddy
[eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
  uvel_tavg,vvel_tavg,theta_tavg,uvelth_tavg,vvelth_tavg, ...
  Nx,Ny,Neta, ...  
  DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
thflux_stand = eflux_mean;
thflux_fluc = eflux_eddy - mean(thflux_eddy,2);
[eflux_tot,eflux_mean,eflux_eddy] = calcMeanEddyFluxes (...
  uvel_tavg,vvel_tavg,salt_tavg,uvelslt_tavg,vvelslt_tavg, ...
  Nx,Ny,Neta, ...  
  DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
sltflux_stand = eflux_mean;
sltflux_fluc = eflux_eddy - mean(sltflux_eddy,2);


%%% Store computed data for later
outfname = [expname,'_TSfluxes'];
if (deform_cavity)
  outfname = [outfname,'_deform'];
elseif (gl_coord)
  outfname = [outfname,'_GLcoord'];
end
outfname = [outfname,'.mat'];
save(fullfile('products',outfname), ...
  'eta','ETA','times', ...
  'thflux_stand','thflux_fluc','thflux_eddy',...
  'sltflux_stand','sltflux_fluc','sltflux_eddy',...
  'uvel_tavg','vvel_tavg','theta_tavg','salt_tavg', ... 
  'uvelth_tavg','vvelth_tavg','uvelslt_tavg','vvelslt_tavg', ...
  'tflux_tavg','sflux_tavg','SHIfwFlx_tavg','SHIhtFlx_tavg','-v7.3');
clear('uvel_tavg','vvel_tavg','theta_tavg','salt_tavg','uvelth_tavg','vvelth_tavg','uvelslt_tavg','vvelslt_tavg');

%%% Add vertical fluxes
wvel_tavg = readIters(exppath,'WVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
wvelth_tavg = readIters(exppath,'WVELTH',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
wvelslt_tavg = readIters(exppath,'WVELSLT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
save(fullfile('products',outfname), ...
  'wvel_tavg','wvelth_tavg','wvelslt_tavg', ...
  '-append','-v7.3');

%%% Add eddy decomposition, if applicable
if (calc_eddy_decomp)
  save(fullfile('products',outfname), ...
    'thflux_eddy_adv','thflux_eddy_stir', ...
    'sltflux_eddy_adv','sltflux_eddy_stir', ...
    '-append','-v7.3');
end








%%%
%%% Convenience function to compute mean and eddy fluxes
%%%
function [trflux_tot,trflux_mean,trflux_eddy] = calcMeanEddyFluxes (...
  uvel,vvel,tracer,uveltr,vveltr, ...
  Nx,Ny,Neta, ...  
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
    Nx,Ny,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  trflux_mean = calcQuasiLatFluxes (...
    uveltr_mean,vveltr_mean, ...
    Nx,Ny,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);
  trflux_eddy = calcQuasiLatFluxes (...
    uveltr_eddy,vveltr_eddy, ...
    Nx,Ny,Neta, ...  
    DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta);

end








%%%
%%% Convenience function to comute fluxes in quasi-latitude space
%%%
function eflux = calcQuasiLatFluxes (...
  uflux,vflux, ...
  Nx,Ny,Neta, ...  
  DXG_3D,DYG_3D,DRF_3D,hFacW,hFacS,ETA,eta)

  %%% Integrate fluxes verticall and horizontally over each cell face
  uflux_yzint = sum(uflux .* DYG_3D .* DRF_3D .* hFacW,3);
  vflux_xzint = sum(vflux .* DXG_3D .* DRF_3D .* hFacS,3);

  %%% Compute horizontal divergence of isopycnal fluxes
  fluxdiv = zeros(Nx,Ny);
  fluxdiv(1:Nx-1,1:Ny-1) = uflux_yzint(2:Nx,1:Ny-1) ...
                              - uflux_yzint(1:Nx-1,1:Ny-1) ...
                              + vflux_xzint(1:Nx-1,2:Ny) ...
                              - vflux_xzint(1:Nx-1,1:Ny-1);
                       
  %%% Integrate flux divergence across lines of constant eta 
  eflux = zeros(Neta,1);
  for m = 1:Neta
    msk = ETA<eta(m);
    eflux(m) = squeeze(sum(sum(fluxdiv.*msk,1),2));
  end  

end