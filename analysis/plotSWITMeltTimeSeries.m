%%%
%%% plotSWITMeltTimeSeries.m
%%%
%%% Plots a time series of melting of the Stancomb-Willis Ice Tonge.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Read experiment data
loadexp;

%%% Vertical grid spacing matrix
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

%%% Diagnostic indix corresponding to instantaneous velocity
diagnum = length(diag_frequency);

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(diagnum);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round(nIter0 + (1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
nDumps = length(dumpIters);

%%% Indices defining area of integration for SWIT
imid_SWIT = find(XC(:,1)>-24,1);
jmax_SWIT = find(hFacC(imid_SWIT,:,1)>0,1);
jmin_SWIT = spongethickness + 1;
imax_SWIT = find(XC(:,1)>-18,1);
imin_SWIT = find(XC(:,1)>-26,1);
jmax_SWIT_edge =  find(YC(1,:)>-73.8,1);
jmin_SWIT_edge =  find(YC(1,:)>-74.3,1);
imax_SWIT_edge = find(XC(:,1)>-23.7,1);
imin_SWIT_edge = find(XC(:,1)>-26,1);

%%% Indices defining area of integration for 
icorner_RL = find(XC(:,1)>-16.56,1);
jcorner_RL = find(YC(1,:)>-72.83,1);
imax_RL = find(XC(:,1)>-15.25,1);
jmin_RL = find(YC(1,:)>-74,1);
imin_RL = find(XC(:,1)>-21,1);
jmax_RL = find(YC(1,:)>-72.38,1);
msk_RL = zeros(Nx,Ny);
msk_RL(imin_RL:icorner_RL,jmin_RL:jmax_RL) = 1;
msk_RL(icorner_RL:imax_RL,jmin_RL:jcorner_RL) = 1;


%%% Physical parameters
rho_i = 920;
Sref = 35;

%%% To store melt rates
tt = nan(1,nDumps);
RL_melt = nan(1,nDumps);
SWIT_melt = nan(1,nDumps);
SWIT_edge_melt = nan(1,nDumps);
SWIT_KE = nan(1,nDumps);
SWIT_temp = nan(1,nDumps);
SWIT_OW = nan(1,nDumps);
SWIT_speed = nan(1,nDumps);
RL_speed = nan(1,nDumps);
SWIT_uvel = nan(1,nDumps);
RL_uvel = nan(1,nDumps);
SWIT_Ro = nan(1,nDumps);

%%% Grid cell volume, only beneath SWIT
vol_msk = RAC.*DRF.*hFacC.*repmat(hFacC(:,:,1)==0,[1 1 Nr]);
SWIT_vol = sum(sum(sum(vol_msk(imin_SWIT:imax_SWIT,jmin_SWIT:jmax_SWIT,:))));

%%% Coriolis parameter
Omega = 2*pi*366/365/86400;
ff = 2*Omega*sind(YG);  
  
%%% Loop through iterations
for n = 1:288
 
  tt(n) =  dumpIters(n)*deltaT;
  tt(n)/t1day
  
  %%% Attempt to load either instantaneous velocities or their squares
  SIarea = rdmdsWrapper(fullfile(exppath,'/results/SIarea_inst'),dumpIters(n)); 
  if (isempty(SIarea))   
    break;
  end  
  SHIfwFlx = rdmdsWrapper(fullfile(exppath,'/results/SHIfwFlx_inst'),dumpIters(n)); 
  if (isempty(SHIfwFlx))   
    break;
  end
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n));      
  if (isempty(uvel))   
    break;
  end
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n)); 
  if (isempty(vvel))   
    break;
  end
  theta = rdmdsWrapper(fullfile(exppath,'/results/THETA_inst'),dumpIters(n));      
  if (isempty(theta))   
    break;
  end
  salt = rdmdsWrapper(fullfile(exppath,'/results/SALT_inst'),dumpIters(n));      
  if (isempty(salt))   
    break;
  end
  

  %%% Slice of horizontal velocity field
  zlev = 1;
  uvel_lev = uvel(:,:,zlev);
  vvel_lev = vvel(:,:,zlev);
  uvel_lev(hFacW(:,:,zlev)==0) = NaN;
  vvel_lev(hFacS(:,:,zlev)==0) = NaN;

  %%% Compute area-averaged melt of SWIT
  massloss = -SHIfwFlx.*RAC;
  SWIT_melt(n) = sum(sum(massloss(imin_SWIT:imax_SWIT,jmin_SWIT:jmax_SWIT)));
  SWIT_edge_melt(n) = sum(sum(massloss(imin_SWIT_edge:imax_SWIT_edge,jmin_SWIT_edge:jmax_SWIT_edge)));
  RL_melt(n) = sum(sum(massloss.*msk_RL));
  
  
  
  %%% Calculate kinetic energy and maximum horizontal speed
  KE = 0.25*(uvel(1:Nx,1:Ny,:).^2+uvel([2:Nx 1],1:Ny,:).^2) ...
     + 0.25*(vvel(1:Nx,1:Ny,:).^2+vvel(1:Nx,[2:Ny 1],:).^2);     
  KE = KE.*vol_msk;
  SWIT_KE(n) = sum(sum(sum(KE(imin_SWIT:imax_SWIT,jmin_SWIT:jmax_SWIT,:))));
  
  %%% Okubo-Weiss on a z-level upstream of SWIT
  dudx = (uvel_lev([2:Nx 1],:) - uvel_lev(1:Nx,:)) ./ DXG; %%% du/dx on cell centers
  dvdy = (vvel_lev(:,[2:Ny 1]) - vvel_lev(:,1:Ny)) ./ DYG; %%% dv/dy on cell centers
  dvdx = (vvel_lev(1:Nx,:) - vvel_lev([Nx 1:Nx-1],:)) ./ DXC; %%% dv/dx on cell corners
  dudy = (uvel_lev(:,1:Ny) - uvel_lev(:,[Ny 1:Ny-1])) ./ DYC; %%% du/du on cell corners
  Sn = dudx - dvdy;
  Ss = dudy + dvdx;
  omega = dvdx - dudy;
  OW = Ss.^2 - omega.^2 + 0.25*(Sn(1:Nx,1:Ny).^2 + Sn([Nx 1:Nx-1],1:Ny).^2 + Sn(1:Nx,[Ny 1:Ny-1]).^2 + Sn([Nx 1:Nx-1],[Ny 1:Ny-1]).^2);
  msk_OW = (XG > -24) & (XG < -22) & (YG > -74) & (YG < -73.5);
  msk_OW(isnan(OW)) = 0;
  OWdist = OW(msk_OW);
  SWIT_OW(n) = prctile(OWdist,5);
%   SWIT_OW(n) = sum(sum(OW.*RAZ.*msk_OW)) / sum(sum(RAZ.*msk_OW));
  
  
  %%% Compute RMS Rossby number next to SWIT
  Ro = omega ./ ff;
  Ro(isnan(Ro)) = 0;
  SWIT_Ro(n) = sqrt(sum(sum(Ro.^2.*RAZ.*msk_OW)) / sum(sum(RAZ.*msk_OW)));
  
  %%% Mean flow speed next to cavities
  uabs_sq = 0.5 * (uvel_lev(1:Nx,1:Ny).^2 + uvel_lev([2:Nx 1],1:Ny).^2) + 0.5 * (vvel_lev(1:Nx,1:Ny).^2 + vvel_lev(1:Nx,[2:Ny 1]).^2);  
  msk_speed_SWIT = (XC > -26) & (XC < -24) & (YC > -74) & (YC < -73.5);
  msk_speed_RL = (XC > -20) & (XC < -18) & (YC > -72.75) & (YC < -72.25);
  msk_speed_SWIT(isnan(uabs_sq)) = 0;
  msk_speed_RL(isnan(uabs_sq)) = 0;
  uabs_sq(isnan(uabs_sq)) = 0;
  SWIT_speed(n) = sqrt(sum(sum(uabs_sq.*RAC.*msk_speed_SWIT)) / sum(sum(RAC.*msk_speed_SWIT)));
  RL_speed(n) = sqrt(sum(sum(uabs_sq.*RAC.*msk_speed_RL)) / sum(sum(RAC.*msk_speed_RL)));
  uvel_lev(isnan(uvel_lev)) = 0;
  SWIT_uvel(n) = sum(sum(uvel_lev.*RAC.*msk_speed_SWIT)) / sum(sum(RAC.*msk_speed_SWIT));
  RL_uvel(n) = sum(sum(uvel_lev.*RAC.*msk_speed_RL)) / sum(sum(RAC.*msk_speed_RL));
    
  %%% Cavity-averaged temperature
  theta = theta.*vol_msk;
  SWIT_temp(n) = sum(sum(sum(theta(imin_SWIT:imax_SWIT,jmin_SWIT:jmax_SWIT,:)))) / SWIT_vol;
   
end

datevec = datenum('01-Jan-2008 00:00') + tt/t1day;

save('./products/SWIT_timeseries.mat',...
  'datevec','SWIT_edge_melt','SWIT_melt','RL_melt','SWIT_KE','SWIT_temp',...
  'SWIT_Ro','SWIT_speed','RL_speed','SWIT_uvel','RL_uvel');


figure(9);
plot(datevec,SWIT_edge_melt/1e12*t1year);
datetick('x');
grid on;

figure(10);
plot(datevec,SWIT_melt/1e12*t1year);
hold on;
plot(datevec,RL_melt/1e12*t1year);
hold off;
datetick('x');
grid on;

figure(11);
plot(datevec,SWIT_KE);
datetick('x');
grid on;

figure(12);
plot(datevec,SWIT_temp);
datetick('x');
grid on;

figure(13);
plotyy(datevec,SWIT_edge_melt/1e12*t1year,datevec,-SWIT_OW);
datetick('x');
grid on;

figure(14);
scatter(SWIT_KE(:),SWIT_melt/1e12*t1year)

figure(15);
plot(datevec,SWIT_Ro);
datetick('x');
grid on;

figure(16);
plot(datevec,SWIT_speed);
datetick('x');
grid on;

figure(17);
plot(datevec,RL_speed);
datetick('x');
grid on;

figure(18);
plot(datevec,SWIT_uvel);
datetick('x');
grid on;

figure(19);
plot(datevec,RL_uvel);
datetick('x');
grid on;
