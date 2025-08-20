%%%
%%% calcFilchnerDSWFlux.m
%%%
%%% Calculates the flux of dense water across the mouth of the Filchner
%%% trough from monthly-averaged isopycnal fluxes.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_nest_onethirtieth_notides_RTOPO2';
tmin = 1.01;
tmax = 7.01;



%%% Load experiment
loadexp;

%%% Select density variable in which to compute isopycnal fluxes
densvar = 'PD0';
% densvar = 'ND1';
% densvar = 'ND2';
% densvar = 'PT';

%%% Density bins for MOC calculation  
switch (densvar)
  case 'PD0'
    % dens_levs = [30.5:1:36.5 36.6:0.1:36.8 36.9:0.01:37.4 37.42:0.02:37.6];
    dens_levs = [30.5:1:36.5 36.6:0.1:36.8 36.9:0.01:37.4 37.42:0.02:37.6] - 9.38; 
  case 'PT'
    dens_levs = [-3:.1:2];
  case 'ND1'
    dens_levs = [21:1:27 27.1:.1:27.5 27.52:.02:28.7 28.74:0.04:29.1];      
  case 'ND2'
    dens_levs = [21:1:27 27.1:.1:27.5 27.52:.02:28.7 28.74:0.04:29.1];      
end
Nd = length(dens_levs)-1;
p_ref = -rhoConst*gravity*RC(1)/1e4; %%% Reference pressure for surface-referenced potential density

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(end));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Indices over which to compute fluxes
ft_lat = -74.54;
ft_w_lon = -35.3;
ft_e_lon = -30.6;
ft_yidx = find(YC(1,:)<ft_lat,1,'last');
ft_w_xidx = find(XC(:,1)>ft_w_lon,1,'first');
ft_e_xidx = find(XC(:,1)>ft_e_lon,1,'first');
Nx_ft = length(ft_w_xidx:ft_e_xidx);


%%% GRIDS %%%

DXG_3D = repmat(DXG,[1 1 Nd]);
DXG_3D = DXG_3D(ft_w_xidx:ft_e_xidx,ft_yidx);








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
T_ft = zeros(Nd+1,Ntime);
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  [num2str(tyears) num2str(itersToRead(n))]      
  
  %%% Read monthly-averaged isopycnal fluxes
  vflux  = rdmdsWrapper(fullfile(exppath,'/results/LaVH1RHO'),itersToRead(n));              
  if (isempty(vflux))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
    break;
  end

  %%% Restrict attention to the Filchner Trough
  vflux = vflux(ft_w_xidx:ft_e_xidx,ft_yidx,:);
  vflux = sum(vflux.*DXG_3D,1); %%% Multiply by grid cell widths to get transports through each density layer in m^3/s
    
  %%% Sum fluxes to obtain transports
  for m=1:Nd  
    T_ft(m,n) = sum(vflux(:,:,m:Nd),3);     
  end

end










%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Construct output file name
outfname = [expname,'_FTtrans_',densvar];
outfname = [outfname,'.mat'];

%%% Store computed data for later
save(fullfile('products',outfname),'-v7.3', ...
  'ft_lat','ft_w_lon','ft_e_lon', ...
  'ft_yidx','ft_w_xidx','ft_e_xidx',... 
  'dens_levs','times','T_ft');





















