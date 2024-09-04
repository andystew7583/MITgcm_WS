%%%
%%% calcPosNegHeatFunction.m
%%% 
%%% Calculates total heat (and salt) functions in quasi-latitude or
%%% streamline coordinates. This script specifically computes the
%%% 'positive' and 'negative' sub-components of the mean component of the heatfunction.
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

%%% Set true to use depth-averaged temperature as the coordinate system
use_meanT = true;

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

    ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
    eta = -9:.1:11;
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

psiT_pos_mean = zeros(Neta,Nr+1,Ntime);
psiT_neg_mean = zeros(Neta,Nr+1,Ntime);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))]      
  
  %%% Read velocity field
  uvel  = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),itersToRead(n));
  vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),itersToRead(n));  
  theta  = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));  
  if (isempty(uvel) || isempty(vvel) || isempty(theta))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
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
 
  %%% Clear memory
  clear('uvel','vvel','theta','theta_pos','theta_neg');
      
end

%%% Store computed data for later
outfname = [expname,'_PosNegHeatFunction'];
if (use_PsiBT)
  outfname = [outfname,'_PsiBT'];
else 
  if (use_meanT)
    outfname = [outfname,'_meanT'];
  else 
    if (deform_cavity)
      outfname = [outfname,'_deform'];
    end
  end
end
outfname = [outfname,'.mat'];
save(fullfile('products',outfname), ...
  'eta','ETA','times', ...  
  'psiT_pos_mean','psiT_neg_mean','-v7.3');





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
