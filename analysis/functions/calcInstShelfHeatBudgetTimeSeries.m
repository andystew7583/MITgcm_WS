%%%
%%% calcInstShelfHeatBudgetTimeSeries.m
%%%
%%% Computes time series of quantities related to continental shelf heat/salt budgets.
%%%

%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;

%%% Index of the upper grid cell face dividing the upper and lower portions
%%% of the water column
zidx_icefront = 25;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to use depth-averaged temperature as the coordinate system
use_meanT = false;

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
  
  %%% Interpolate tot cell centers
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

    %%% Bounds and horizontal mask for heat budget volume
    eta_icefront = -1.1;
    eta_shelfbreak = 3.5;             

  end

end

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
%%% For daily/12-hourly outputs
dumpStart = 1578240;
dumpStep = 86400/2/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;
itersToRead = dumpIters;
times = dumpIters*deltaT;
Ntime = length(itersToRead);


%%% Calendar information
days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays = [0 cumsum(days_per_month)];
Nmonths = length(days_per_month);





%%% Storage
tt = zeros(1,Nmonths);
hflux_icefront = zeros(1,Nmonths);
sflux_icefront = zeros(1,Nmonths);
mflux_icefront = zeros(1,Nmonths);
hflux_shelfbreak = zeros(1,Nmonths);
sflux_shelfbreak = zeros(1,Nmonths);
mflux_shelfbreak = zeros(1,Nmonths);
wt = zeros(1,Nmonths);
ws = zeros(1,Nmonths);
wm = zeros(1,Nmonths);
ttend = zeros(1,Nmonths);
stend = zeros(1,Nmonths);
tlen = zeros(1,Nmonths);

%%% Loop over iterations
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  tt(n) = itersToRead(n)*deltaT;
  disp([n num2str(tyears) num2str(itersToRead(n))])
    
  %%% Load instantaneous fields  
  salt = rdmdsWrapper(fullfile(exppath,'/results/SALT_12hourly'),itersToRead(n));      
  theta = rdmdsWrapper(fullfile(exppath,'/results/THETA_12hourly'),itersToRead(n));

  %%% Number of days into the year
  nday = (n-1)/2;



  %%% TENDENCIES %%%

  %%% If this snapshot corresponds to the start or end of a month, add its
  %%% contribution to the monthly tendency of volume-integrated T/S
  month_end_match = find(nday == cumdays);
  if (~isempty(month_end_match))
    msk = (ETA > eta_icefront) & (ETA < eta_shelfbreak);
    theta_zint = sum(theta(:,:,zidx_icefront+1:end).*DRF(:,:,zidx_icefront+1:end).*hFacC(:,:,zidx_icefront+1:end),3);
    theta_int = sum(sum(theta_zint.*RAC.*msk));
    salt_zint = sum(salt(:,:,zidx_icefront+1:end).*DRF(:,:,zidx_icefront+1:end).*hFacC(:,:,zidx_icefront+1:end),3);
    salt_int = sum(sum(salt_zint.*RAC.*msk));
    if (month_end_match > 1)
      ttend(month_end_match-1) = ttend(month_end_match-1) + theta_int / (days_per_month(month_end_match-1)*86400);
      stend(month_end_match-1) = stend(month_end_match-1) + salt_int / (days_per_month(month_end_match-1)*86400);
    end
    if (month_end_match < Nmonths+1)
      ttend(month_end_match) = ttend(month_end_match) - theta_int / (days_per_month(month_end_match)*86400);
      stend(month_end_match) = stend(month_end_match) - salt_int / (days_per_month(month_end_match)*86400);
    end    
  end



  %%% HORIZONTAL FLUXES %%%

  %%% Load horizontal velocity field
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),itersToRead(n));      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),itersToRead(n));
    
  %%% Compute depth-integrated heat fluxes at shelf break and ice front
  uvelth = uvel .* 0.5.*(theta([1:Nx],:,:)+theta([Nx 1:Nx-1],:,:));  
  uflux_yzint = sum(uvelth(:,:,zidx_icefront:end) .* DYG .* DRF(:,:,zidx_icefront:end) .* hFacW(:,:,zidx_icefront:end),3);
  clear('uvelth');
  vvelth = vvel .* 0.5.*(theta(:,[1:Ny],:)+theta(:,[Ny 1:Ny-1],:));   
  vflux_xzint = sum(vvelth(:,:,zidx_icefront:end) .* DXG .* DRF(:,:,zidx_icefront:end) .* hFacS(:,:,zidx_icefront:end),3);
  clear('vvelth');    
  fluxdiv = zeros(Nx,Ny);
  fluxdiv(1:Nx-1,1:Ny-1) = uflux_yzint(2:Nx,1:Ny-1) ...
                          - uflux_yzint(1:Nx-1,1:Ny-1) ...
                          + vflux_xzint(1:Nx-1,2:Ny) ...
                          - vflux_xzint(1:Nx-1,1:Ny-1);                        
  msk = ETA<eta_icefront;
  hflux_i = squeeze(sum(sum(fluxdiv.*msk,1),2));  
  msk = ETA<eta_shelfbreak;
  hflux_s = squeeze(sum(sum(fluxdiv.*msk,1),2));  

  %%% Compute depth-integrated salt fluxes at shelf break and ice front
  uvelslt = uvel .* 0.5.*(salt([1:Nx],:,:)+salt([Nx 1:Nx-1],:,:));  
  uflux_yzint = sum(uvelslt(:,:,zidx_icefront:end) .* DYG .* DRF(:,:,zidx_icefront:end) .* hFacW(:,:,zidx_icefront:end),3);
  clear('uvelslt');
  vvelslt = vvel .* 0.5.*(salt(:,[1:Ny],:)+salt(:,[Ny 1:Ny-1],:));   
  vflux_xzint = sum(vvelslt(:,:,zidx_icefront:end) .* DXG .* DRF(:,:,zidx_icefront:end) .* hFacS(:,:,zidx_icefront:end),3);
  clear('vvelslt');  
  fluxdiv = zeros(Nx,Ny);
  fluxdiv(1:Nx-1,1:Ny-1) = uflux_yzint(2:Nx,1:Ny-1) ...
                          - uflux_yzint(1:Nx-1,1:Ny-1) ...
                          + vflux_xzint(1:Nx-1,2:Ny) ...
                          - vflux_xzint(1:Nx-1,1:Ny-1);                        
  msk = ETA<eta_icefront;
  sflux_i = squeeze(sum(sum(fluxdiv.*msk,1),2));  
  msk = ETA<eta_shelfbreak;
  sflux_s = squeeze(sum(sum(fluxdiv.*msk,1),2));  

  %%% Compute depth-integrated volume fluxes at shelf break and ice front  
  uflux_yzint = sum(uvel(:,:,zidx_icefront:end) .* DYG .* DRF(:,:,zidx_icefront:end) .* hFacW(:,:,zidx_icefront:end),3);  
  vflux_xzint = sum(vvel(:,:,zidx_icefront:end) .* DXG .* DRF(:,:,zidx_icefront:end) .* hFacS(:,:,zidx_icefront:end),3);  
  fluxdiv = zeros(Nx,Ny);
  fluxdiv(1:Nx-1,1:Ny-1) = uflux_yzint(2:Nx,1:Ny-1) ...
                          - uflux_yzint(1:Nx-1,1:Ny-1) ...
                          + vflux_xzint(1:Nx-1,2:Ny) ...
                          - vflux_xzint(1:Nx-1,1:Ny-1);                        
  msk = ETA<eta_icefront;
  mflux_i = squeeze(sum(sum(fluxdiv.*msk,1),2));  
  msk = ETA<eta_shelfbreak;
  mflux_s = squeeze(sum(sum(fluxdiv.*msk,1),2));  


  %%% Clear memory
  clear('uvel','vvel');



  %%% VERTICAL FLUXES %%%

  %%% Load vertical velocity
  wvel = rdmdsWrapper(fullfile(exppath,'/results/WVEL_12hourly'),itersToRead(n));

  %%% Compute area-integrated vertical heat flux
  wt_tmp = wvel(:,:,zidx_icefront) .* 0.5.*(theta(:,:,zidx_icefront)+theta(:,:,zidx_icefront-1));
  ws_tmp = wvel(:,:,zidx_icefront) .* 0.5.*(salt(:,:,zidx_icefront)+salt(:,:,zidx_icefront-1));
  msk = (ETA > eta_icefront) & (ETA < eta_shelfbreak);
  wt_tmp = sum(sum(wt_tmp.*RAC.*msk));
  ws_tmp = sum(sum(ws_tmp.*RAC.*msk));
  wm_tmp = sum(sum(wvel(:,:,zidx_icefront).*RAC.*msk));

  clear('wvel','theta','salt');

  

  %%% ASSIGNMENT TO MONTHLY AVERAGES

  %%% Assign fluxes to monthly averages based on centered averaging in time
  nmonth = find(nday>=cumdays,1,'last');
  if (nday == cumdays(nmonth))
    if (nmonth < Nmonths+1)
      hflux_shelfbreak(nmonth) = hflux_shelfbreak(nmonth) + hflux_s/2;
      hflux_icefront(nmonth) = hflux_icefront(nmonth) + hflux_i/2;
      sflux_shelfbreak(nmonth) = sflux_shelfbreak(nmonth) + sflux_s/2;
      sflux_icefront(nmonth) = sflux_icefront(nmonth) + sflux_i/2;
      mflux_shelfbreak(nmonth) = mflux_shelfbreak(nmonth) + mflux_s/2;
      mflux_icefront(nmonth) = mflux_icefront(nmonth) + mflux_i/2;
      wt(nmonth) = wt(nmonth) + wt_tmp/2;
      ws(nmonth) = ws(nmonth) + ws_tmp/2;
      wm(nmonth) = wm(nmonth) + wm_tmp/2;
      tlen(nmonth) = tlen(nmonth) + 0.5;
    end
    if (nmonth > 1)
      hflux_shelfbreak(nmonth-1) = hflux_shelfbreak(nmonth-1) + hflux_s/2;
      hflux_icefront(nmonth-1) = hflux_icefront(nmonth-1) + hflux_i/2;
      sflux_shelfbreak(nmonth-1) = sflux_shelfbreak(nmonth-1) + sflux_s/2;
      sflux_icefront(nmonth-1) = sflux_icefront(nmonth-1) + sflux_i/2;
      mflux_shelfbreak(nmonth-1) = mflux_shelfbreak(nmonth-1) + mflux_s/2;
      mflux_icefront(nmonth-1) = mflux_icefront(nmonth-1) + mflux_i/2;
      wt(nmonth-1) = wt(nmonth-1) + wt_tmp/2;
      ws(nmonth-1) = ws(nmonth-1) + ws_tmp/2;
      wm(nmonth-1) = wm(nmonth-1) + wm_tmp/2;
      tlen(nmonth-1) = tlen(nmonth-1) + 0.5;
    end
  else
    hflux_shelfbreak(nmonth) = hflux_shelfbreak(nmonth) + hflux_s;
    hflux_icefront(nmonth) = hflux_icefront(nmonth) + hflux_i;
    sflux_shelfbreak(nmonth) = sflux_shelfbreak(nmonth) + sflux_s;
    sflux_icefront(nmonth) = sflux_icefront(nmonth) + sflux_i;
    mflux_shelfbreak(nmonth) = mflux_shelfbreak(nmonth) + mflux_s;
    mflux_icefront(nmonth) = mflux_icefront(nmonth) + mflux_i;
    wt(nmonth) = wt(nmonth) + wt_tmp;
    ws(nmonth) = ws(nmonth) + ws_tmp;
    wm(nmonth) = wm(nmonth) + wm_tmp;
    tlen(nmonth) = tlen(nmonth) + 1;
  end
   
end

%%% Finish time averaging
hflux_icefront = hflux_icefront ./ tlen;
hflux_shelfbreak = hflux_shelfbreak ./ tlen;
sflux_icefront = sflux_icefront ./ tlen;
sflux_shelfbreak = sflux_shelfbreak ./ tlen;
mflux_icefront = mflux_icefront ./ tlen;
mflux_shelfbreak = mflux_shelfbreak ./ tlen;
wt = wt ./ tlen;
ws = ws ./ tlen;
wm = wm ./ tlen;

%%% Write to output file
outfname = [expname,'_InstShelfHeatBudget.mat'];
save(fullfile('products',outfname), ...
    'tt', ...
    'hflux_icefront','hflux_shelfbreak',...
    'sflux_icefront','sflux_shelfbreak',...
    'mflux_icefront','mflux_shelfbreak',...
    'wt','ws','wm', ...
    'ttend','stend', ...
    '-v7.3');
