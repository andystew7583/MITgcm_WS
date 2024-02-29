%%%
%%% createTRMvelocityFile.m
%%%
%%% Creates a 3D velocity file with the components of the TRM velocity.
%%%

%%% Options
expdir = '../experiments';
expname = 'hires_seq_onethird_RTOPO2';
tmin = 18.05;
tmax = 27.05;
% expname = 'hires_seq_onesixth_RTOPO2';
% tmin = 9.05;
% tmax = 18.05;
% expname = 'hires_seq_onetwelfth_RTOPO2';
% tmin = 1.05;
% tmax = 9.05;
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% tmin = 1.01;
% tmax = 7.01;

%%% Load experiment
loadexp;

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










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TRM VELOCITY CALCULATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load required time-mean model outputs
u_mean  = readIters(exppath,'UVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  
v_mean  = readIters(exppath,'VVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);     
t_mean  = readIters(exppath,'THETA',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  
s_mean  = readIters(exppath,'SALT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  
w_mean  = readIters(exppath,'WVEL',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  
ut_tavg = readIters(exppath,'UVELTH',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  
vt_tavg = readIters(exppath,'VVELTH',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  
wt_tavg = readIters(exppath,'WVELTH',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  
us_tavg = readIters(exppath,'UVELSLT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  
vs_tavg = readIters(exppath,'VVELSLT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  
ws_tavg = readIters(exppath,'WVELSLT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);  

%%% Compute eddy-induced velocities
[u_eddy,v_eddy,w_eddy] = calcTRMvelocity (...
  uvel_mean,vvel_mean,w_mean,t_mean,s_mean, ...
  ut_tavg,vt_tavg,wt_tavg, ...
  us_tavg,vs_tavg,ws_tavg, ...
  hFacC,hFacW,hFacS, ...
  DXG,DYG,RAC,DXC,DYC, ...
  DRF,DRC,RC,RF,...
  rhoConst,gravity);











%%%%%%%%%%%%%%%%%%%%
%%% WRITE OUTPUT %%%
%%%%%%%%%%%%%%%%%%%%



