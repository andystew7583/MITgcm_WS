%%%
%%% calcMeanOverturning.m
%%%
%%% Calculates the Eulerian-mean overturning circulation in quasi-latitude 
%%% coordinates.
%%%

%%% Options
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% tmin = 18.05;
% tmax = 19.05;
% expname = 'hires_seq_onesixth_RTOPO2';
% tmin = 9.05*86400*365;
% tmax = 18.05*86400*365;
expname = 'hires_seq_onetwelfth_notides_RTOPO2';
tmin = 1.05;
tmax = 2.05;

%%% Load experiment
loadexp;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Define coordinate system for integrating to compute streamfunction
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
eta = -9:.1:11;
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
%%% STREAMFUNCTION CALCULATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate time-averaged isopycnal flux, density and velocity
psi_EM = zeros(Neta,Nr+1,Ntime);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))]      
  
  %%% Read velocity field
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),itersToRead(n));
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),itersToRead(n));              
  if (isempty(uvel) || isempty(vvel))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end

  %%% Integrate fluxes verticall and horizontally over each cell face
  uflux_yzint = uvel .* DYG_3D .* DRF_3D .* hFacW;
  vflux_xzint = vvel .* DXG_3D .* DRF_3D .* hFacS;

  %%% Compute horizontal divergence of isopycnal fluxes
  fluxdiv = zeros(Nx,Ny,Nr);
  fluxdiv(1:Nx-1,1:Ny-1,:) = uflux_yzint(2:Nx,1:Ny-1,:) ...
                              - uflux_yzint(1:Nx-1,1:Ny-1,:) ...
                              + vflux_xzint(1:Nx-1,2:Ny,:) ...
                              - vflux_xzint(1:Nx-1,1:Ny-1,:);
                       
  %%% Integrate flux divergence across lines of constant eta (parallel to FRIS face)
  eflux = zeros(Neta,Nr);
  for m = 1:Neta
    msk = repmat(ETA<eta(m),[1 1 Nr]);
    eflux(m,:) = squeeze(sum(sum(fluxdiv.*msk,1),2));
  end

  %%% Sum fluxes to obtain streamfunction
  for k=1:Nr  
    psi_EM(:,k,n) = -sum(eflux(:,k:Nr),2);     
  end
 

end









%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Construct output file name
outfname = [expname,'_EMOC'];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];

%%% Store computed data for later
save(fullfile('products',outfname), ...
  'eta','ETA','RF','times', ...
  'psi_EM');













