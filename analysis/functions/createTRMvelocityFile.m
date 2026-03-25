%%%
%%% createTRMvelocityFile.m
%%%
%%% Creates a 3D velocity file with the components of the TRM velocity.
%%%

%%% Options
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% tmin = 18.05;
% tmax = 27.05;
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% tmin = 9.05;
% tmax = 18.05;
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
% tmin = 1.05;
% tmax = 9.05;
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
tmin = 1.01;
tmax = 7.01;

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
  u_mean,v_mean,w_mean,t_mean,s_mean, ...
  ut_tavg,vt_tavg,wt_tavg, ...
  us_tavg,vs_tavg,ws_tavg, ...
  hFacC,hFacW,hFacS, ...
  DXG,DYG,RAC,DXC,DYC, ...
  DRF,DRC,RC,RF,...
  rhoConst,gravity);

%%% Free up memory
clear('ut_tavg','vt_tavg','wt_tavg','us_tavg','vs_tavg','ws_tavg');

%%% Assemble mean+eddy TRM velocities
u_tot = u_mean + u_eddy;
v_tot = v_mean + v_eddy;
w_tot = w_mean + w_eddy;











%%%%%%%%%%%%%%%%%%%%
%%% WRITE OUTPUT %%%
%%%%%%%%%%%%%%%%%%%%

ncfname = fullfile('products',[expname,'_TRM.nc']);
nccreate(ncfname,'XC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'XC',XC);
nccreate(ncfname,'YC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'YC',YC);
nccreate(ncfname,'XG','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'XG',XG);
nccreate(ncfname,'YG','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'YG',YG);
nccreate(ncfname,'DXC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'DXC',DXC);
nccreate(ncfname,'DYC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'DYC',DYC);
nccreate(ncfname,'DXG','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'DXG',DXG);
nccreate(ncfname,'DYG','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'DYG',DYG);
nccreate(ncfname,'DRF','Dimensions',{'z',Nr},'FillValue','disable');
ncwrite(ncfname,'DRF',squeeze(DRF));
nccreate(ncfname,'DRC','Dimensions',{'zf',Nr+1},'FillValue','disable');
ncwrite(ncfname,'DRC',squeeze(DRC));
nccreate(ncfname,'RC','Dimensions',{'z',Nr},'FillValue','disable');
ncwrite(ncfname,'RC',squeeze(RC));
nccreate(ncfname,'RF','Dimensions',{'zf',Nr+1},'FillValue','disable');
ncwrite(ncfname,'RF',squeeze(RF));
nccreate(ncfname,'bathy','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'bathy',bathy);
nccreate(ncfname,'SHELFICEtopo','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'SHELFICEtopo',SHELFICEtopo);
nccreate(ncfname,'TRM_U','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'TRM_U',u_tot);
nccreate(ncfname,'TRM_V','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'TRM_V',v_tot);
nccreate(ncfname,'TRM_W','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'TRM_W',w_tot);
nccreate(ncfname,'TRM_U_MEAN','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'TRM_U_MEAN',u_mean);
nccreate(ncfname,'TRM_V_MEAN','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'TRM_V_MEAN',v_mean);
nccreate(ncfname,'TRM_W_MEAN','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'TRM_W_MEAN',w_mean);
nccreate(ncfname,'TRM_U_EDDY','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'TRM_U_EDDY',u_eddy);
nccreate(ncfname,'TRM_V_EDDY','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'TRM_V_EDDY',v_eddy);
nccreate(ncfname,'TRM_W_EDDY','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'TRM_W_EDDY',w_eddy);
nccreate(ncfname,'THETA','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'THETA',t_mean);
nccreate(ncfname,'SALT','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
ncwrite(ncfname,'SALT',s_mean);


