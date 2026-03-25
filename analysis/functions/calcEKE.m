%%%
%%% calcEKE.m
%%%
%%% Calculates EKE and EKE production.
%%%

%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_RTOPO2';
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
expname = 'hires_nest_onethirtieth_notides_RTOPO2';
loadexp;

%%% Time frame in model years
% tmin = 19.05*86400*365;
% tmax = 27.05*86400*365;
% tmin = 10.05*86400*365;
% tmax = 18.05*86400*365;
% tmin = 1.05*86400*365;
% tmax = 9.05*86400*365;
tmin = 0.05;
tmax = 7.05;

%%% Required for vertical integrals
DRF3D = repmat(DRF,[Nx Ny 1]);
Hocean = sum(DRF3D.*hFacC,3);





%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%

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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATE PRODUCTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Storage
tt = zeros(1,Ntime);
EKE_surf = zeros(Nx,Ny,Ntime);
EKE_50m = zeros(Nx,Ny,Ntime);
EKE_zavg = zeros(Nx,Ny,Ntime);
PEtoEKE_zavg = zeros(Nx,Ny,Ntime);
MKEtoEKE_zavg = zeros(Nx,Ny,Ntime);
tlen = 0;

%%% Loop over iterations
for n=1:Ntime
 
  %%% Print current time to keep track of calculation
  tyears = itersToRead(n)*deltaT/86400/365; 
  tt(n) = itersToRead(n)*deltaT;
  [num2str(tyears) num2str(itersToRead(n))] 





  %%% Calculate EKE
  
  %%% Time-average required fields
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),itersToRead(n));   
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),itersToRead(n));   
  uvelsq = rdmdsWrapper(fullfile(exppath,'/results/UVELSQ'),itersToRead(n));   
  vvelsq = rdmdsWrapper(fullfile(exppath,'/results/VVELSQ'),itersToRead(n));   
  
  %%% Calculate EKE
  uvelsq_eddy = uvelsq - uvel.^2;
  vvelsq_eddy = vvelsq - vvel.^2;
  EKE = 0.5 * ( 0.5*(uvelsq_eddy(1:Nx,1:Ny,:)+uvelsq_eddy([2:Nx 1],1:Ny,:)) ...
              + 0.5*(vvelsq_eddy(1:Nx,1:Ny,:)+vvelsq_eddy(1:Nx,[2:Ny 1],:)) );
  EKE_zavg(:,:,n) = sum(EKE.*DRF3D.*hFacC,3) ./ Hocean;
  EKE_50m(:,:,n) = sum(EKE(:,:,1:9).*DRF3D(:,:,1:9).*hFacC(:,:,1:9),3) ./ sum(DRF3D(:,:,1:9).*hFacC(:,:,1:9),3);
  EKE_surf(:,:,n) = EKE(:,:,1);
  clear('uvelsq','vvelsq');
  
  
  
  
  %%% Calculate MKE->EKE
  
  %%% Time-average required fields
  wvel =  rdmdsWrapper(fullfile(exppath,'/results/WVEL'),itersToRead(n));   
  
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
  
  %%% Add component involving u'v'
  uwvel =  rdmdsWrapper(fullfile(exppath,'/results/WU_VEL'),itersToRead(n));   
  uwvel_mean = 0.5.*(uvel(:,:,1:Nr)+uvel(:,:,[Nr 1:Nr-1])) ...
             .* 0.5.*(wvel(1:Nx,:,:)+wvel([Nx 1:Nx-1],:,:));
  uwvel_eddy = uwvel - uwvel_mean;
  clear('uwvel','uwvel_mean');
  MKEtoEKE = MKEtoEKE ...        
            - uwvel_eddy.*(uvel(:,:,[Nr 1:Nr-1])-uvel(:,:,1:Nr))./repmat(DRC(1:Nr),[Nx Ny 1]);
  clear('uwvel_eddy','uvel')
  
  %%% Add component involving v'w'
  vwvel =  rdmdsWrapper(fullfile(exppath,'/results/WV_VEL'),itersToRead(n));   
  vwvel_mean = 0.5.*(vvel(:,:,1:Nr)+vvel(:,:,[Nr 1:Nr-1])) ...
             .* 0.5.*(wvel(:,1:Ny,:)+wvel(:,[Ny 1:Ny-1],:));        
  vwvel_eddy = vwvel - vwvel_mean;
  clear('vwvel','vwvel_mean');    
  MKEtoEKE = MKEtoEKE ...        
            - vwvel_eddy.*(vvel(:,:,[Nr 1:Nr-1])-vvel(:,:,1:Nr))./repmat(DRC(1:Nr),[Nx Ny 1]);
  clear('vwvel_eddy','vvel');
  
  %%% Integrate vertically
  MKEtoEKE_zavg(:,:,n) = nansum(MKEtoEKE.*DRF3D.*hFacC,3) ./ Hocean;
  
  
  
  
  
  
  
  %%% Calculate vertical eddy buoyancy flux
  
  %%% Salinity on w-points
  salt = rdmdsWrapper(fullfile(exppath,'/results/SALT'),itersToRead(n));   
  salt_w = 0*ones(Nx,Ny,Nr);
  salt_w(:,:,2:Nr) = 0.5*(salt(:,:,1:Nr-1)+salt(:,:,2:Nr));
  salt_w(:,:,1) = salt(:,:,1);
  clear('salt');
  
  %%% Vertical eddy salt flux
  wvelslt =  rdmdsWrapper(fullfile(exppath,'/results/WVELSLT'),itersToRead(n));   
  wvelslt_eddy = 0*ones(Nx,Ny,Nr);
  wvelslt_eddy(:,:,1:Nr) = wvelslt - wvel .* salt_w(:,:,1:Nr);
  clear('wvelslt');
  
  %%% Temperature on w-points
  theta =  rdmdsWrapper(fullfile(exppath,'/results/THETA'),itersToRead(n));   
  theta_w = 0*ones(Nx,Ny,Nr);
  theta_w(:,:,2:Nr) = 0.5*(theta(:,:,1:Nr-1)+theta(:,:,2:Nr));
  theta_w(:,:,1) = theta(:,:,1);
  clear('theta');
  
  %%% Vertical eddy heat flux
  wvelth =  rdmdsWrapper(fullfile(exppath,'/results/WVELTH'),itersToRead(n));   
  wvelth_eddy = 0*ones(Nx,Ny,Nr);
  wvelth_eddy(:,:,1:Nr) = wvelth - wvel .* theta_w(:,:,1:Nr); 
  clear('wvelth','wvel');
  
  %%% Compute thermal expansion and haline contraction coefficients
  press_w = -rhoConst*gravity*repmat(RF(1:Nr),[Nx Ny 1])/1e4;
  [alpha_w,beta_w] = calcAlphaBeta(salt_w,theta_w,press_w);
  clear('salt_w','theta_w','press_w');
  
  %%% Compute baroclinic energy production
  PEtoEKE = gravity*(alpha_w.*wvelth_eddy - beta_w.*wvelslt_eddy);
  clear('alpha_w','beta_w','wvelth_eddy','wvelslt_eddy');
  
  %%% Integrate vertically
  PEtoEKE_zavg(:,:,n) = sum(PEtoEKE.*DRF3D.*hFacC,3) ./ Hocean;
 



end

%%% Catch infinite values occurring where water column thickness is zero
EKE_zavg(isinf(EKE_zavg)) = 0;
MKEtoEKE_zavg(isinf(MKEtoEKE_zavg)) = 0;
PEtoEKE_zavg(isinf(PEtoEKE_zavg)) = 0;




%%%%%%%%%%%%%%%%%%
%%%%% OUTPUT %%%%%
%%%%%%%%%%%%%%%%%%

outfname = [expname,'_EKE.mat'];
save(fullfile('products',outfname),'EKE_surf','EKE_zavg','EKE_50m','PEtoEKE_zavg','MKEtoEKE_zavg','Hocean','-v7.3');

