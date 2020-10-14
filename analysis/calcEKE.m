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
expname = 'hires_seq_onesixth_RTOPO2';
loadexp;

%%% Time frame over which to average thermodynamic variables to create
%%% climatology
% tmin = 18.05*86400*365;
% tmax = 27.05*86400*365;
tmin = 9.05*86400*365;
tmax = 18.05*86400*365;








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

%%% Time-average required fields
uvel = readIters(exppath,'UVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
vvel = readIters(exppath,'VVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
wvel = readIters(exppath,'WVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
salt = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
theta = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
wvelslt = readIters(exppath,'WVELSLT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
wvelth = readIters(exppath,'WVELTH',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
uvelsq = readIters(exppath,'UVELSQ',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
vvelsq = readIters(exppath,'VVELSQ',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATE PRODUCTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DRF3D = repmat(DRF,[Nx Ny 1]);

%%% Calculate EKE
uvelsq_eddy = uvelsq - uvel.^2;
vvelsq_eddy = vvelsq - vvel.^2;
EKE = 0.5 * ( 0.5*(uvelsq_eddy(1:Nx,1:Ny,:)+uvelsq_eddy([2:Nx 1],1:Ny,:)) ...
            + 0.5*(vvelsq_eddy(1:Nx,1:Ny,:)+vvelsq_eddy(1:Nx,[2:Ny 1],:)) );
EKE_zavg = sum(EKE.*DRF3D.*hFacC,3) ./ sum(DRF3D.*hFacC,3);
EKE_zavg(isinf(EKE_zavg)) = 0;

%%% Calculate vertical eddy buoyancy flux
salt_w = 0*ones(Nx,Ny,Nr);
salt_w(:,:,2:Nr) = 0.5*(salt(:,:,1:Nr-1)+salt(:,:,2:Nr));
salt_w(:,:,1) = salt(:,:,1);
theta_w = 0*ones(Nx,Ny,Nr);
theta_w(:,:,2:Nr) = 0.5*(theta(:,:,1:Nr-1)+theta(:,:,2:Nr));
theta_w(:,:,1) = theta(:,:,1);
press_w = -rhoConst*gravity*repmat(RF(1:Nr),[Nx Ny 1])/1e4;
[alpha_w,beta_w] = calcAlphaBeta(salt_w,theta_w,press_w);
wvelslt_eddy = 0*ones(Nx,Ny,Nr);
wvelslt_eddy(:,:,1:Nr) = wvelslt - wvel .* salt_w(:,:,1:Nr);
wvelth_eddy = 0*ones(Nx,Ny,Nr);
wvelth_eddy(:,:,1:Nr) = wvelth - wvel .* theta_w(:,:,1:Nr); 
PEtoEKE = gravity*(alpha_w.*wvelth_eddy - beta_w.*wvelslt_eddy);
PEtoEKE_zavg = sum(PEtoEKE.*DRF3D.*hFacC,3) ./ sum(DRF3D.*hFacC,3);
PEtoEKE_zavg(isinf(PEtoEKE_zavg)) = 0;












%%%%%%%%%%%%%%%%%%
%%%%% OUTPUT %%%%%
%%%%%%%%%%%%%%%%%%

outfname = [expname,'_EKE.mat'];
save(fullfile('products',outfname),'EKE','EKE_zavg','PEtoEKE','PEtoEKE_zavg');

