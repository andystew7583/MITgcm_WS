%%%
%%% calcOverturning_totOnly.m
%%%
%%% Calculates the overturning circulation in density surfaces, similar to 
%%% that calculated using the MITgcm 'layers' package. This script only
%%% calculates the full MOC from the LAYERS package, with no decomposition
%%% into mean/eddy components.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2_SSH';
tmin = 0.01;
tmax = 1.01;

%%% Load experiment
loadexp;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use output from the Layers package to calculate isopycnal
%%% fluxes. N.B. if this option is selected then the density variable must
%%% be 'PD0' (surface-referenced potential density)
use_layers = true;

%%% Define density variable
densvar = 'PD0';
dens_levs = layers_bounds;
Nd = length(dens_levs)-1;
p_ref = -rhoConst*gravity*RC(1)/1e4; %%% Reference pressure for surface-referenced potential density

%%% MOC coordinate
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
eta = -9:.1:11;
Neta = length(eta);


%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(end));
if (dumpFreq == 2592000)  
  startDate = num2str(startDate_1);
  startYear = str2num(startDate(1:4));
  startMonth = str2num(startDate(5:6));
  startday = str2num(startDate(7:8));
  startDate = datenum([startDate(1:4),'-',startDate(5:6),'-',startDate(7:8)]);
  endDate = startDate+endTime/86400;
  currentYear = startYear;
  currentMonth = startMonth;
  currentDate = startDate;
  dumpIters = [];
  while (true)
    currentMonth = currentMonth + 1;
    if (currentMonth > 12)
      currentMonth = 1;
      currentYear = currentYear + 1;      
    end
    currentDate = datenum([num2str(currentYear),'-',num2str(currentMonth),'-01']);
    if (currentDate > endDate)
      break;
    end
    dumpIters = [dumpIters ceil((currentDate-startDate)*86400/deltaT)];
  end
  nDumps = length(dumpIters);
else
  nDumps = round(endTime/dumpFreq);
  dumpIters = round((1:nDumps)*dumpFreq/deltaT);
  dumpIters = dumpIters(dumpIters > nIter0);
  nDumps = length(dumpIters);
end







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
psi_tot = zeros(Neta,Nd+1,Ntime);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))]      
    
  %%% Read isopycnal fluxes
  uflux  = rdmdsWrapper(fullfile(exppath,'/results/LaUH1RHO'),itersToRead(n));
  vflux  = rdmdsWrapper(fullfile(exppath,'/results/LaVH1RHO'),itersToRead(n));              
  if (isempty(uflux) || isempty(vflux))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end

  %%% Compute total isopycnal streamfunction
  psi_tot(:,:,n) = calcIsopStreamfunction(...
    uflux,vflux, ...
    Nx,Ny,Neta,Nd, ...  
    DXG_3D,DYG_3D,ETA,eta);
      
end





%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Construct output file name
outfname = [expname,'_MOCtotOnly_',densvar];
estr = '_layers';
outfname = [outfname,estr];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];

%%% Store computed data for later
save(fullfile('products',outfname),'-v7.3', ...
  'eta','ETA','dens_levs','times', ...
  'psi_tot');























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















