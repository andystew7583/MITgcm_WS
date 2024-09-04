%%%
%%% calcShelfHeatBudgetTimeSeries.m
%%%
%%% Computes time series of quantities related to continental shelf heat/salt budgets.
%%%

%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
expname = 'hires_seq_onethird_RTOPO2';
tmin = 19.05;
tmax = 27.05;
% expname = 'hires_seq_onesixth_RTOPO2';
% tmin = 10.05;
% tmax = 18.05;
% expname = 'hires_seq_onetwelfth_RTOPO2';
% tmin = 1.05;
% tmax = 9.05;
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% tmin = 1.05;
% tmax = 7.05;
loadexp;

%%% Index of the upper grid cell face dividing the upper and lower portions
%%% of the water column
zidx_icefront = 25;

%%% Reference surface freezing temperature
theta0 = -1.9;


%%% Required for vertical integrals
DRF3D = repmat(DRF,[Nx Ny 1]);

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
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








%%% Storage
tt = zeros(1,Ntime);
usq_eddy_int = zeros(Nx,Ny,Ntime);
vsq_eddy_int = zeros(Nx,Ny,Ntime);
salt_int_lower = zeros(Nx,Ny,Ntime);
salt_int_upper = zeros(Nx,Ny,Ntime);
theta_int_lower = zeros(Nx,Ny,Ntime);
theta_int_upper = zeros(Nx,Ny,Ntime);
theta_pos_int_lower = zeros(Nx,Ny,Ntime);
theta_pos_int_upper = zeros(Nx,Ny,Ntime);
theta_neg_int_lower = zeros(Nx,Ny,Ntime);
theta_neg_int_upper = zeros(Nx,Ny,Ntime);
salt_bnd = zeros(Nx,Ny,Ntime);
theta_bnd = zeros(Nx,Ny,Ntime);
theta_pos_bnd = zeros(Nx,Ny,Ntime);
theta_neg_bnd = zeros(Nx,Ny,Ntime);
PEtoEKE_int = zeros(Nx,Ny,Ntime);
MKEtoEKE_int = zeros(Nx,Ny,Ntime);
wt_tot_flux = zeros(Nx,Ny,Ntime);
ws_tot_flux = zeros(Nx,Ny,Ntime);
wt_mean_flux = zeros(Nx,Ny,Ntime);
ws_mean_flux = zeros(Nx,Ny,Ntime);
wt_eddy_flux = zeros(Nx,Ny,Ntime);
ws_eddy_flux = zeros(Nx,Ny,Ntime);
w_flux = zeros(Nx,Ny,Ntime);
tlen = 0;

%%% Loop over iterations
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  tt(n) = itersToRead(n)*deltaT;
  [num2str(tyears) num2str(itersToRead(n))] 



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% VELOCITY-RELATED DIAGNOSTICS %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Depth-integrate zonal component of EKE
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),itersToRead(n));      
  uvelsq = rdmdsWrapper(fullfile(exppath,'/results/UVELSQ'),itersToRead(n)); 
  uvelsq_eddy = uvelsq - uvel.^2;  
  clear('uvelsq');
  usq_eddy_int(:,:,n) = sum(0.5.*uvelsq_eddy.*DRF.*hFacW,3);  
  
  %%% Depth-integrate meridional component of EKE
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),itersToRead(n));      
  vvelsq = rdmdsWrapper(fullfile(exppath,'/results/VVELSQ'),itersToRead(n)); 
  vvelsq_eddy = vvelsq - vvel.^2;
  clear('vvelsq');         
  vsq_eddy_int(:,:,n) = sum(0.5.*vvelsq_eddy.*DRF.*hFacS,3);  

  %%% Calculate MKE->EKE
  
  %%% Remove land points
  uvel(hFacW==0) = NaN;
  vvel(hFacS==0) = NaN;
  
  %%% Calculate components involving u'^2 and v'^2
  MKEtoEKE = - uvelsq_eddy.*(uvel([2:Nx 1],:,:)-uvel(:,:,:))./DXG ...
             - vvelsq_eddy.*(vvel(:,[2:Ny 1],:)-vvel(:,:,:))./DYG;
  clear('uvelsq_eddy','vvelsq_eddy');
  
  %%% Add components involving u'v'
  uvvel = rdmdsWrapper(fullfile(exppath,'/results/UV_VEL_Z'),itersToRead(n));
  uvvel_mean = 0.5.*(vvel(1:Nx,1:Ny,:)+vvel([Nx 1:Nx-1],1:Ny,:)) ...
          .* 0.5.*(uvel(1:Nx,1:Ny,:)+uvel(1:Nx,[Ny 1:Ny-1],:));
  uvvel_eddy = uvvel - uvvel_mean;
  clear('uvvel','uvvel_mean');
  MKEtoEKE = MKEtoEKE ...
             - uvvel_eddy.*(uvel(:,1:Ny,:)-uvel(:,[Ny 1:Ny-1],:))./DYC ...
             - uvvel_eddy.*(vvel(1:Nx,:,:)-vvel([Nx 1:Nx-1],:,:))./DXC;
  clear('uvvel_eddy');

  %%% Load vertical velocity
  wvel = rdmdsWrapper(fullfile(exppath,'/results/WVEL'),itersToRead(n));

  %%% Vertical transport across boundary
  w_flux(:,:,n) = wvel(:,:,zidx_icefront);
  
  %%% Add component involving u'w'
  uwvel = rdmdsWrapper(fullfile(exppath,'/results/WU_VEL'),itersToRead(n));
  uwvel_mean = 0.5.*(uvel(:,:,1:Nr)+uvel(:,:,[Nr 1:Nr-1])) ...
             .* 0.5.*(wvel(1:Nx,:,:)+wvel([Nx 1:Nx-1],:,:));
  uwvel_eddy = uwvel - uwvel_mean;
  clear('uwvel','uwvel_mean');
  MKEtoEKE = MKEtoEKE ...        
            - uwvel_eddy.*(uvel(:,:,[Nr 1:Nr-1])-uvel(:,:,1:Nr))./repmat(DRC(1:Nr),[Nx Ny 1]);
  clear('uwvel_eddy','uvel')
  
  %%% Add component involving v'w'
  vwvel = rdmdsWrapper(fullfile(exppath,'/results/WV_VEL'),itersToRead(n));
  vwvel_mean = 0.5.*(vvel(:,:,1:Nr)+vvel(:,:,[Nr 1:Nr-1])) ...
             .* 0.5.*(wvel(:,1:Ny,:)+wvel(:,[Ny 1:Ny-1],:));        
  vwvel_eddy = vwvel - vwvel_mean;
  clear('vwvel','vwvel_mean');    
  MKEtoEKE = MKEtoEKE ...        
            - vwvel_eddy.*(vvel(:,:,[Nr 1:Nr-1])-vvel(:,:,1:Nr))./repmat(DRC(1:Nr),[Nx Ny 1]);
  clear('vwvel_eddy','vvel');
  
  %%% Integrate vertically
  MKEtoEKE_tmp = nansum(MKEtoEKE.*DRF3D.*hFacC,3);
  MKEtoEKE_tmp(isinf(MKEtoEKE_tmp)) = 0;  
  MKEtoEKE_int(:,:,n) = MKEtoEKE_tmp;

  %%% Free up memory
  clear('MKEtoEKE');



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% SALINITY-RELATED DIAGNOSTICS %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Load salinity
  salt = rdmdsWrapper(fullfile(exppath,'/results/SALT'),itersToRead(n));      

  %%% Compute depth-integrated salinity diagnostics
  salt_int_lower(:,:,n) = sum(salt(:,:,zidx_icefront:end).*DRF(:,:,zidx_icefront:end).*hFacC(:,:,zidx_icefront:end),3);
  salt_int_upper(:,:,n) = sum(salt(:,:,1:zidx_icefront-1).*DRF(:,:,1:zidx_icefront-1).*hFacC(:,:,1:zidx_icefront-1),3);  

  %%% Salinity on w-points
  salt_w = 0*ones(Nx,Ny,Nr);
  salt_w(:,:,2:Nr) = 0.5*(salt(:,:,1:Nr-1)+salt(:,:,2:Nr));
  salt_w(:,:,1) = salt(:,:,1);

  %%% Salinity at the boundary
  salt_bnd(:,:,n) = salt_w(:,:,zidx_icefront);
    
  %%% Free up memory
  clear('salt');
  
  %%% Vertical salt flux
  wvelslt = rdmdsWrapper(fullfile(exppath,'/results/WVELSLT'),itersToRead(n));
  
  %%% Mean and eddy components
  wvelslt_mean = wvel.*salt_w;
  wvelslt_eddy = wvelslt - wvelslt_mean;

  %%% Salt fluxes across boundary
  ws_tot_flux(:,:,n) = wvelslt(:,:,zidx_icefront);
  ws_mean_flux(:,:,n) = wvelslt_mean(:,:,zidx_icefront);
  ws_eddy_flux(:,:,n) = wvelslt_eddy(:,:,zidx_icefront);
 
  %%% Free up memory 
  clear('wvelslt','wvelslt_mean');



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% TEMPERATURE-RELATED DIAGNOSTICS %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Temperature on w-points
  theta = rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));

  %%% Positive and negative temperature anomalies
  theta_pos = theta - theta0;
  theta_neg = min(theta_pos,0);
  theta_pos = max(theta_pos,0);  

  %%% Compute depth-integrated temperature diagnostics
  theta_int_lower(:,:,n) = sum(theta(:,:,zidx_icefront:end).*DRF(:,:,zidx_icefront:end).*hFacC(:,:,zidx_icefront:end),3);
  theta_int_upper(:,:,n) = sum(theta(:,:,1:zidx_icefront-1).*DRF(:,:,1:zidx_icefront-1).*hFacC(:,:,1:zidx_icefront-1),3);
  theta_pos_int_lower(:,:,n) = sum(theta_pos(:,:,zidx_icefront:end).*DRF(:,:,zidx_icefront:end).*hFacC(:,:,zidx_icefront:end),3);
  theta_pos_int_upper(:,:,n) = sum(theta_pos(:,:,1:zidx_icefront-1).*DRF(:,:,1:zidx_icefront-1).*hFacC(:,:,1:zidx_icefront-1),3);
  theta_neg_int_lower(:,:,n) = sum(theta_neg(:,:,zidx_icefront:end).*DRF(:,:,zidx_icefront:end).*hFacC(:,:,zidx_icefront:end),3);
  theta_neg_int_upper(:,:,n) = sum(theta_neg(:,:,1:zidx_icefront-1).*DRF(:,:,1:zidx_icefront-1).*hFacC(:,:,1:zidx_icefront-1),3);
  clear('theta_pos','theta_neg');

  %%% Temperature on w-points
  theta_w = 0*ones(Nx,Ny,Nr);
  theta_w(:,:,2:Nr) = 0.5*(theta(:,:,1:Nr-1)+theta(:,:,2:Nr));
  theta_w(:,:,1) = theta(:,:,1);

  %%% Positive and negative temperature anomalies
  theta_pos_w = theta_w - theta0;
  theta_neg_w = min(theta_pos_w,0);
  theta_pos_w = max(theta_pos_w,0); 

  %%% Temperature at the boundary
  theta_bnd(:,:,n) = theta_w(:,:,zidx_icefront);
  theta_pos_bnd(:,:,n) = theta_pos_w(:,:,zidx_icefront);
  theta_neg_bnd(:,:,n) = theta_neg_w(:,:,zidx_icefront);

  %%% Free up memory
  clear('theta','theta_pos_w','theta_neg_w');
  
  %%% Vertical eddy heat flux
  wvelth = rdmdsWrapper(fullfile(exppath,'/results/WVELTH'),itersToRead(n));
  wvelth_mean = wvel.*theta_w;
  wvelth_eddy = wvelth - wvelth_mean;

  %%% Heat fluxes across boundary
  wt_tot_flux(:,:,n) = wvelth(:,:,zidx_icefront);
  wt_mean_flux(:,:,n) = wvelth_mean(:,:,zidx_icefront);
  wt_eddy_flux(:,:,n) = wvelth_eddy(:,:,zidx_icefront);

  %%% Free up memeory
  clear('wvelth','wvelth_mean','wvel');



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% BUOYANCY-RELATED DIAGNOSTICS %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Compute thermal expansion and haline contraction coefficients
  press_w = -rhoConst*gravity*repmat(RF(1:Nr),[Nx Ny 1])/1e4;
  [alpha_w,beta_w] = calcAlphaBeta(salt_w,theta_w,press_w);

  %%% Free up memory
  clear('salt_w','theta_w','press_w');
  
  %%% Compute baroclinic energy production
  PEtoEKE = gravity*(alpha_w.*wvelth_eddy - beta_w.*wvelslt_eddy);

  %%% Free up memory
  clear('alpha_w','beta_w','wvelth_eddy','wvelslt_eddy');
  
  %%% Integrate vertically
  PEtoEKE_tmp = sum(PEtoEKE.*DRF3D.*hFacC,3);
  PEtoEKE_tmp(isinf(PEtoEKE_tmp)) = 0;
  PEtoEKE_int(:,:,n) = PEtoEKE_tmp;

  %%% Free up mmeory
  clear('PEtoEKE');



  tlen = tlen + 1;
  
end


%%% Write to output file
outfname = [expname,'_ShelfHeatBudget.mat'];
save(fullfile('products',outfname), ...
    'tt','usq_eddy_int','vsq_eddy_int',...
    'salt_int_lower','salt_int_upper',...
    'theta_int_lower','theta_int_upper',...
    'theta_pos_int_lower','theta_pos_int_upper',...
    'theta_neg_int_lower','theta_neg_int_upper',...
    'salt_bnd','theta_bnd','theta_pos_bnd','theta_neg_bnd',...
    'PEtoEKE_int','MKEtoEKE_int', ...
    'wt_tot_flux','wt_tot_flux', ...
    'ws_tot_flux','ws_tot_flux', ...
    'wt_mean_flux','wt_mean_flux', ...
    'ws_mean_flux','ws_mean_flux', ...
    'wt_eddy_flux','wt_eddy_flux', ...
    'ws_eddy_flux','ws_eddy_flux', ...
    'w_flux', ...
    '-v7.3');
