%%%
%%% calcTSstreamfunction.m
%%%
%%% Calculates thermohaline streamfunction in T/S space.
%%%

%%% Read experiment data
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_RTOPO2';
expname = 'hires_seq_onetwelfth_RTOPO2';
expdir = '/data3/MITgcm_WS/experiments';
loadexp;
% tmin = 18.05;
% tmax = 27.05;
% tmin = 9.05;
% tmax = 18.05;
tmin = 1.05;
tmax = 9.05;

%%% Set true to calculate streamfunction by integrating w.r.t. S as well as
%%% w.r.t. T
calc_intS = false;

%%% Set true to compute eddy-induced transports
calc_psi_eddy = true; 

%%% T/S grids
dT = 0.05;
dS = 0.01;
Tmin = -3;
Tmax = 2;
Smin = 33.9;
Smax = 35;
SS = Smin:dS:Smax;
TT = Tmin:dT:Tmax;
NS = length(SS);
NT = length(TT);

%%% Experiment output iterations
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Determine iteraton numbers to process
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

%%% Grids
DXG = rdmds(fullfile(resultspath,'DXG'));
DYG = rdmds(fullfile(resultspath,'DYG'));
RAC = rdmds(fullfile(resultspath,'RAC'));
DRF = rdmds(fullfile(resultspath,'DRF'));
hFacS = rdmds(fullfile(resultspath,'hFacS'));
hFacW = rdmds(fullfile(resultspath,'hFacW'));
hFacC = rdmds(fullfile(resultspath,'hFacC'));

%%% Areas of cell faces
Axz = repmat(reshape(DRF,[1 1 Nr]),[Nx Ny]).*DXG.*hFacS;
Ayz = repmat(reshape(DRF,[1 1 Nr]),[Nx Ny]).*DYG.*hFacW;
Axy = RAC;

%%% To store fluxes in TS space
psi_TS_mean_intT = zeros(NS,NT,Ntime);
psi_TS_mean_intS = zeros(NS,NT,Ntime);
psi_TS_eddy_intT = zeros(NS,NT,Ntime);
psi_TS_eddy_intS = zeros(NS,NT,Ntime);
uvel_tavg = zeros(Nx,Ny,Nr);
vvel_tavg = zeros(Nx,Ny,Nr);
wvel_tavg = zeros(Nx,Ny,Nr);
salt_tavg = zeros(Nx,Ny,Nr);
theta_tavg = zeros(Nx,Ny,Nr);
for n=1:Ntime

  %%% To keep track of progress
  tyears = itersToRead(n)*deltaT/86400/365

  %%% Read model output
  theta = rdmdsWrapper(fullfile(exppath,'results','THETA'),itersToRead(n));
  salt = rdmdsWrapper(fullfile(exppath,'results','SALT'),itersToRead(n));
  uvel = rdmdsWrapper(fullfile(exppath,'results','UVEL'),itersToRead(n));
  vvel = rdmdsWrapper(fullfile(exppath,'results','VVEL'),itersToRead(n));
  wvel = rdmdsWrapper(fullfile(exppath,'results','WVEL'),itersToRead(n));

  %%% Compute mean streamfunction
  tstart = tic();
  psi_TS_mean_intT(:,:,n) = calcStreamfunction_intT (...
    uvel,vvel,wvel,salt,theta, ...
    Nx,Ny,Nr,NS,NT,SS,TT,Ayz,Axz,Axy);
  disp(toc(tstart))
  
  if (calc_intS)
    psi_TS_mean_intS(:,:,n) = calcStreamfunction_intS (...
      uvel,vvel,wvel,salt,theta, ...
      Nx,Ny,Nr,NS,NT,SS,TT,Ayz,Axz,Axy);
  end
  
  %%% Calculate eddy streamfunction, if required
  if (calc_psi_eddy)

    %%% Load heat and salt fluxes     
    uvelth = rdmdsWrapper(fullfile(exppath,'/results/UVELTH'),itersToRead(n));
    vvelth = rdmdsWrapper(fullfile(exppath,'/results/VVELTH'),itersToRead(n));
    wvelth = rdmdsWrapper(fullfile(exppath,'/results/WVELTH'),itersToRead(n));
    uvelslt = rdmdsWrapper(fullfile(exppath,'/results/UVELSLT'),itersToRead(n));
    vvelslt = rdmdsWrapper(fullfile(exppath,'/results/VVELSLT'),itersToRead(n));
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
    clear('uvelth','vvelth','wvelth','uvelslt','vvelslt','wvelslt');

    %%% Compute streamfunction
    psi_TS_eddy_intT(:,:,n) = calcStreamfunction_intT (...
      u_eddy,v_eddy,w_eddy,salt,theta, ...
      Nx,Ny,Nr,NS,NT,SS,TT,Ayz,Axz,Axy);
    if (calc_intS)
      psi_TS_eddy_intS(:,:,n) = calcStreamfunction_intS (...
        u_eddy,v_eddy,w_eddy,salt,theta, ...
        Nx,Ny,Nr,NS,NT,SS,TT,Ayz,Axz,Axy);
    end


  end
  
  %%% Add to time averages
  uvel_tavg = uvel_tavg + uvel;
  vvel_tavg = vvel_tavg + vvel;
  wvel_tavg = wvel_tavg + wvel;
  salt_tavg = salt_tavg + salt;
  theta_tavg = theta_tavg + theta;
 
end

%%% Divide by Ntime to get time mean
uvel_tavg = uvel_tavg / Ntime;
vvel_tavg = vvel_tavg / Ntime;
wvel_tavg = wvel_tavg / Ntime;
salt_tavg = salt_tavg / Ntime;
theta_tavg = theta_tavg / Ntime;

%%% Make sure no vertical velocities are present at solid cell faces
wvel(hFacC==0) = 0;

%%% Compute "standing" component of streamfunction
psi_TS_mean_stand_intT = calcStreamfunction_intT (...
    uvel_tavg,vvel_tavg,wvel_tavg,salt_tavg,theta_tavg, ...
    Nx,Ny,Nr,NS,NT,SS,TT,Ayz,Axz,Axy);
if (calc_intS)
  psi_TS_mean_stand_intS = calcStreamfunction_intS (...
      uvel_tavg,vvel_tavg,wvel_tavg,salt_tavg,theta_tavg, ...
      Nx,Ny,Nr,NS,NT,SS,TT,Ayz,Axz,Axy);
else
  psi_TS_mean_stand_intS = zeros(NS,NT);
end
  
%%% Compute "fluctuating" component of streamfunction
psi_TS_mean_fluc_intT = mean(psi_TS_mean_intT,3) - psi_TS_mean_stand_intT;
psi_TS_mean_fluc_intS = mean(psi_TS_mean_intS,3) - psi_TS_mean_stand_intS;








%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Construct output file name
outfname = [expname,'_MOC_TS'];
if (calc_psi_eddy)  
  estr = '_TRM';  
else
  estr = '_noeddy';
end
outfname = [outfname,estr];
outfname = [outfname,'.mat'];

%%% Store computed data for later
save(fullfile('products',outfname),'SS','TT','times', ...
  'psi_TS_mean_intT','psi_TS_mean_stand_intT','psi_TS_mean_fluc_intT', ...
  'psi_TS_mean_intS','psi_TS_mean_stand_intS','psi_TS_mean_fluc_intS' ...
);

if (calc_psi_eddy)
  save(fullfile('products',outfname), ...
  'psi_TS_eddy_intT','psi_TS_eddy_intS', ...
  '-append');
end









%%%
%%% Deal with open boundary transports by wrapping them around the boundary
%%% artificially, closing the domain-integrated transport. Iterate around
%%% open boundary points computing flux divergence at each point, and
%%% balancing it with a boundary-parallel barotropic flow.
%%%
function [Tu,Tv,Tw] = wrapBdryTransport (...
  uvel,vvel,wvel, ...
  Nx,Ny,Nr,Ayz,Axz,Axy)
  
  %%% BOUNDARY-MEAN FLOW %%%
  
  NEtransport = sum(sum(vvel(1:Nx-1,Ny,:).*Axz(1:Nx-1,Ny,:))) ...
              + sum(sum(uvel(Nx,1:Ny-1,:).*Ayz(Nx,1:Ny-1,:)));              
  NEarea = sum(sum(Axz(1:Nx-1,Ny,:))) + sum(sum(Ayz(Nx,1:Ny-1,:)));
  Ubdry = NEtransport / NEarea;
  
  %%% NORTHERN BOUNDARY %%%
  
  for i=1:Nx-1

    %%% Compute N boundary flux divergence
    fluxconv = sum(vvel(i,Ny,:).*Axz(i,Ny,:) - Ubdry.*Axz(i,Ny,:) + uvel(i,Ny,:).*Ayz(i,Ny,:) - uvel(i+1,Ny,:).*Ayz(i+1,Ny,:), 3);
    wcarea = sum(Ayz(i+1,Ny,:),3);

    %%% Avoid NaNs
    if (wcarea > 0)

      %%% Add barotropic correction to eliminate boundary flux divergence
      uanom = fluxconv / wcarea;
      uvel(i+1,Ny,:) = uvel(i+1,Ny,:) + uanom;

      %%% Recompute flux convergence in order to recompute vertical velocity
      fluxconv = vvel(i,Ny,:).*Axz(i,Ny,:) - Ubdry.*Axz(i,Ny,:) + uvel(i,Ny,:).*Ayz(i,Ny,:) - uvel(i+1,Ny,:).*Ayz(i+1,Ny,:);
      wvel(i,Ny,2:Nr) = - cumsum(fluxconv(:,:,1:Nr-1),3) / Axy(i,Ny);

    end

  end
  
  %%% NORTHEAST CORNER %%%
  
  %%% Compute NE corner flux divergence
  fluxconv = sum(vvel(Nx,Ny,:).*Axz(Nx,Ny,:) - Ubdry.*Axz(Nx,Ny,:) + uvel(Nx,Ny,:).*Ayz(Nx,Ny,:) - Ubdry.*Ayz(Nx,Ny,:), 3);

  %%% Avoid NaNs
  wcarea = sum(Axz(Nx,Ny,:),3);
  if (wcarea > 0) 

    %%% Add barotropic correction to eliminate boundary flux divergence
    vanom = fluxconv / wcarea;
    vvel(Nx,Ny,:) = vvel(Nx,Ny,:) - vanom;

    %%% Recompute flux convergence in order to recompute vertical velocity
    fluxconv = vvel(Nx,Ny,:).*Axz(Nx,Ny,:) - Ubdry.*Axz(Nx,Ny,:) + uvel(Nx,Ny,:).*Ayz(Nx,Ny,:) - Ubdry.*Ayz(Nx,Ny,:);
    wvel(Nx,Ny,2:Nr) = - cumsum(fluxconv(:,:,1:Nr-1),3) / Axy(Nx,Ny);

  end
  
  %%% EASTERN BOUNDARY %%%
  
  for j=Ny-1:-1:1

    %%% Compute E boundary flux divergence
    fluxconv = sum(vvel(Nx,j,:).*Axz(Nx,j,:) - vvel(Nx,j+1,:).*Axz(Nx,j+1,:) + uvel(Nx,j,:).*Ayz(Nx,j,:) - Ubdry.*Ayz(Nx,j,:), 3);

    %%% Avoid NaNs
    wcarea = sum(Axz(Nx,j,:),3);
    if (wcarea > 0)

      %%% Add barotropic correction to eliminate boundary flux divergence
      vanom = fluxconv / wcarea;
      vvel(Nx,j,:) = vvel(Nx,j,:) - vanom;

      %%% Recompute flux convergence in order to recompute vertical velocity
      fluxconv = vvel(Nx,j,:).*Axz(Nx,j,:) - vvel(Nx,j+1,:).*Axz(Nx,j+1,:) + uvel(Nx,j,:).*Ayz(Nx,j,:) - Ubdry.*Ayz(Nx,j,:);
      wvel(Nx,j,2:Nr) = - cumsum(fluxconv(:,:,1:Nr-1),3) / Axy(Nx,j);

    end

  end

  %%% Transport through grid cell faces
  Tu = Ayz.*uvel;
  Tv = Axz.*vvel;
  Tw = Axy.*wvel;
  
   
  %%% DEBUG CODE
  %%% Map vertically integrated flux convergence
  % fluxconv = zeros(Nx,Ny);
  % fluxconv(1:Nx,1:Ny) = fluxconv(1:Nx,1:Ny) + sum(Tu(1:Nx,1:Ny,:) + Tv(1:Nx,1:Ny,:), 3);
  % fluxconv(1:Nx-1,1:Ny) = fluxconv(1:Nx-1,1:Ny) - sum(Tu(2:Nx,1:Ny,:),3);
  % fluxconv(1:Nx,1:Ny-1) = fluxconv(1:Nx,1:Ny-1) - sum(Tv(1:Nx,2:Ny,:),3);
  % figure(100);
  % pcolor(XC,YC,fluxconv);
  % shading interp;
  % colorbar;

  %%% DEBUG CODE
%   fluxdiv3D = zeros(Nx,Ny,Nr);
%   fluxdiv3D = fluxdiv3D - Tu - Tv + Tw;
%   fluxdiv3D(1:Nx-1,1:Ny,1:Nr) = fluxdiv3D(1:Nx-1,1:Ny,1:Nr) + Tu(2:Nx,1:Ny,:);
%   fluxdiv3D(1:Nx,1:Ny-1,1:Nr) = fluxdiv3D(1:Nx,1:Ny-1,1:Nr) + Tv(1:Nx,2:Ny,:);
%   fluxdiv3D(1:Nx,1:Ny,1:Nr-1) = fluxdiv3D(1:Nx,1:Ny,1:Nr-1) - Tw(1:Nx,1:Ny,2:Nr);
%   sum(sum(sum(fluxdiv3D(theta<-1.75))))

end












%%%
%%% Computes TS streamfunction by integrating w.r.t. salinity at constant
%%% temperature.
%%%
function psi_TS = calcStreamfunction_intS (...
  uvel,vvel,wvel,salt,theta, ...
  Nx,Ny,Nr,NS,NT,SS,TT,Ayz,Axz,Axy)


  %%% Calculate midpoint salinity 
  salt_u = 0.5*(salt([1:Nx],:,:)+salt([Nx 1:Nx-1],:,:));
  salt_v = 0.5*(salt(:,[1:Ny],:)+salt(:,[Ny 1:Ny-1],:));
  salt_w = 0.5*(salt(:,:,[1:Nr])+salt(:,:,[Nr 1:Nr-1]));

  %%% Wrap open boundary transports to close whole-domain transport
  [Tu,Tv,Tw] = wrapBdryTransport(uvel,vvel,wvel,Nx,Ny,Nr,Ayz,Axz,Axy);

  %%% Loop over temperature space
  psi_TS = zeros(NS,NT);
  for nt = 1:NT

    nt

    %%% Create mask over all physical velocity points, leaving non-zero values 
    %%% Only where adjacent temperature points straddle the current
    %%% temperature (TT(nt)), and signed +1 (-1) if the temperature gradient
    %%% is positive (negative)
    msk_u_theta =  ((theta([1:Nx],:,:)-TT(nt)) .* (theta([Nx 1:Nx-1],:,:)-TT(nt)) < 0) ...
                .* sign(theta([1:Nx],:,:)-theta([Nx 1:Nx-1],:,:)); %%% Flow crosses temp class TT(nt), positive if dT/dx>0
    msk_v_theta =  ((theta(:,[1:Ny],:)-TT(nt)) .* (theta(:,[Ny 1:Ny-1],:)-TT(nt)) < 0) ...
                .* sign(theta(:,[1:Ny],:,:)-theta(:,[Ny 1:Ny-1],:)); %%% Flow crosses temp class TT(nt), positive if dT/dy>0
    msk_w_theta =  ((theta(:,:,[Nr 1:Nr-1])-TT(nt)) .* (theta(:,:,[1:Nr])-TT(nt)) < 0) ...
                .* sign(theta(:,:,[Nr 1:Nr-1])-theta(:,:,[1:Nr])); %%% Flow crosses temp class TT(nt), positive if dT/dz>0
              
    %%% Precompute product of transport with theta mask        
    Tu_msktheta = Tu.*msk_u_theta;    
    Tv_msktheta = Tv.*msk_v_theta;
    Tw_msktheta = Tw.*msk_w_theta;

    %%% Loop over salinity space
    for ns = 1:NS-1  

      %%% Create masks over all physical velocity points that are non-zero
      %%% only where the salinity falls within the current salinity bin
%       msk_u_salt = (salt_u>SS(ns)) & (salt_u<SS(ns+1)); %%% Salinity falls in salt class [SS(ns) SS(ns+1)]     
%       msk_v_salt = (salt_v>SS(ns)) & (salt_v<SS(ns+1));
%       msk_w_salt = (salt_w>SS(ns)) & (salt_w<SS(ns+1));
      msk_u_salt = (salt_u<SS(ns)); %%% Salinity falls in below SS(ns)
      msk_v_salt = (salt_v<SS(ns));
      msk_w_salt = (salt_w<SS(ns));

      %%% Flux in T/S space is equal to sum of physical fluxes over all 
      %%% unmasked points     
%       F_TS(ns,nt,n) = sum(sum(sum(Tu_msktheta.*msk_u_salt + Tv_msktheta.*msk_v_salt + Tw_msktheta.*msk_w_salt)));
      psi_TS(ns,nt) = sum(sum(sum(Tu_msktheta.*msk_u_salt + Tv_msktheta.*msk_v_salt + Tw_msktheta.*msk_w_salt)));      

    end
  end

  
  %%% Streamfunction is equal to sum of fluxes
  % psi_TS(2:NS,:,:) = cumsum(F_TS,1);
  
end












% %%%
% %%% Computes TS streamfunction by integrating w.r.t. temperature at constant
% %%% salinity.
% %%%
% function psi_TS = calcStreamfunction_intT (...
%   uvel,vvel,wvel,salt,theta, ...
%   Nx,Ny,Nr,NS,NT,SS,TT,Ayz,Axz,Axy)
% 
%   %%% Calculate midpoint temperature
%   theta_u = 0.5*(theta([1:Nx],:,:)+theta([Nx 1:Nx-1],:,:));
%   theta_v = 0.5*(theta(:,[1:Ny],:)+theta(:,[Ny 1:Ny-1],:));
%   theta_w = 0.5*(theta(:,:,[1:Nr])+theta(:,:,[Nr 1:Nr-1]));
% 
%   %%% Wrap open boundary transports to close whole-domain transport
%   [Tu,Tv,Tw] = wrapBdryTransport(uvel,vvel,wvel,Nx,Ny,Nr,Ayz,Axz,Axy);
%   
%   %%% Loop over salinity space
%   psi_TS = zeros(NS,NT);
%   for ns = 1:NS 
%     ns
% 
%     %%% Create mask over all physical velocity points, leaving non-zero values 
%     %%% Only where adjacent salinity points straddle the current
%     %%% salinity (SS(ns)), and signed +1 (-1) if the temperature gradient
%     %%% is positive (negative)
%     msk_u_salt =  ((salt([1:Nx],:,:)-SS(ns)) .* (salt([Nx 1:Nx-1],:,:)-SS(ns)) < 0) ...
%                 .* sign(salt([1:Nx],:,:)-salt([Nx 1:Nx-1],:,:)); %%% Flow crosses temp class SS(ns), positive if dT/dx>0
%     msk_v_salt =  ((salt(:,[1:Ny],:)-SS(ns)) .* (salt(:,[Ny 1:Ny-1],:)-SS(ns)) < 0) ...
%                 .* sign(salt(:,[1:Ny],:,:)-salt(:,[Ny 1:Ny-1],:)); %%% Flow crosses temp class SS(ns), positive if dT/dy>0
%     msk_w_salt =  ((salt(:,:,[Nr 1:Nr-1])-SS(ns)) .* (salt(:,:,[1:Nr])-SS(ns)) < 0) ...
%                 .* sign(salt(:,:,[Nr 1:Nr-1])-salt(:,:,[1:Nr])); %%% Flow crosses temp class SS(ns), positive if dT/dz>0
%               
%     %%% Precompute product of transport with salt mask        
%     Tu_msksalt = Tu.*msk_u_salt;    
%     Tv_msksalt = Tv.*msk_v_salt;
%     Tw_msksalt = Tw.*msk_w_salt;
% 
%     %%% Loop over temperature space
%     for nt = 1:NT-1
% 
%       %%% Create masks over all physical velocity points that are non-zero
%       %%% only where the temperature falls within the current temperature bin
%       msk_u_temp = (theta_u<TT(nt)); 
%       msk_v_temp = (theta_v<TT(nt));
%       msk_w_temp = (theta_w<TT(nt));
% 
%       %%% Flux in T/S space is equal to sum of physical fluxes over all 
%       %%% unmasked points    
%       psi_TS(ns,nt) = sum(sum(sum(Tu_msksalt.*msk_u_temp + Tv_msksalt.*msk_v_temp + Tw_msksalt.*msk_w_temp)));      
% 
%     end
%   end
% 
% end













%%%
%%% Computes TS streamfunction by integrating w.r.t. temperature at constant
%%% salinity.
%%%
function psi_TS = calcStreamfunction_intT (...
  uvel,vvel,wvel,salt,theta, ...
  Nx,Ny,Nr,NS,NT,SS,TT,Ayz,Axz,Axy)

  %%% Calculate midpoint temperature
  theta_u = 0.5*(theta([1:Nx],:,:)+theta([Nx 1:Nx-1],:,:));
  theta_v = 0.5*(theta(:,[1:Ny],:)+theta(:,[Ny 1:Ny-1],:));
  theta_w = 0.5*(theta(:,:,[1:Nr])+theta(:,:,[Nr 1:Nr-1]));

  %%% Wrap open boundary transports to close whole-domain transport
  [Tu,Tv,Tw] = wrapBdryTransport(uvel,vvel,wvel,Nx,Ny,Nr,Ayz,Axz,Axy);
  
  %%% Loop over physical space and assign fluxes
  psi_TS = zeros(NS,NT);
  for i=1:Nx
    i
    for j=1:Ny
      for k=1:Nr
        
        %%% u-fluxes
        if ((i > 1) && (Tu(i,j,k)~=0))
 
          %%% Local temperature and adjacent salinities
          Tval = theta_u(i,j,k);
          Sm = salt(i-1,j,k);
          Sp = salt(i,j,k);
          
          %%% Add flux to all temperatures bins smaller than local
          %%% temperature, and to all salinity bins between adjacent
          %%% salinities
          nt = find(Tval<TT);
          ns = find((Sp-SS) .* (Sm-SS) < 0); %%% N.B. can be multiple indices
          psi_TS(ns,nt) = psi_TS(ns,nt) + Tu(i,j,k)*sign(Sp-Sm);
          
        end        
        
        %%% v-fluxes
        if ((j > 1) && (Tv(i,j,k)~=0))
          
          %%% Local temperature and adjacent salinities
          Tval = theta_v(i,j,k);
          Sm = salt(i,j-1,k);
          Sp = salt(i,j,k);
          
          %%% Add flux to all temperatures bins smaller than local
          %%% temperature, and to all salinity bins between adjacent
          %%% salinities
          nt = find(Tval<TT);
          ns = find((Sp-SS) .* (Sm-SS) < 0); %%% N.B. can be multiple indices
          psi_TS(ns,nt) = psi_TS(ns,nt) + Tv(i,j,k)*sign(Sp-Sm);
          
        end
        
        %%% w-fluxes
        if ((k > 1) && (Tw(i,j,k)~=0))
          
          %%% Local temperature and adjacent salinities
          Tval = theta_w(i,j,k);
          Sm = salt(i,j,k);
          Sp = salt(i,j,k-1);
          
          %%% Add flux to all temperatures bins smaller than local
          %%% temperature, and to all salinity bins between adjacent
          %%% salinities
          nt = find(Tval<TT);
          ns = find((Sp-SS) .* (Sm-SS) < 0); %%% N.B. can be multiple indices
          psi_TS(ns,nt) = psi_TS(ns,nt) + Tw(i,j,k)*sign(Sp-Sm);
          
        end
        
      end
    end
  end

end