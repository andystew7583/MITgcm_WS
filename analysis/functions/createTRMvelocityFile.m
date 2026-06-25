%%%
%%% createTRMvelocityFile.m
%%%
%%% Creates a 3D velocity file with the components of the TRM velocity.
%%%

%%% Options
expdir = '../experiments';
% expname = 'hires_seq_onethird_RTOPO2';
% startyr = 2007;
% model_yrs = 20:27;
% cycle_len = 9;
% expname = 'hires_seq_onesixth_RTOPO2';
% model_yrs = 11:18;
% startyr = 2007;
% cycle_len = 9;
expname = 'hires_seq_onetwelfth_RTOPO2';
model_yrs = 2:9;
startyr = 2007;
cycle_len = 9;
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% model_yrs = 2:7;
% cycle_len = 7;

%%% Load experiment
loadexp;
% startyr = num2str(startDate_1); %%% Starting year for the simulation
% startyr = str2num(startyr(1:4));
Nyrs = length(model_yrs);
cal_yrs = startyr + mod(model_yrs-1,cycle_len);

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);



%%% TODO year vectors

%%% Create NC file
ncfname = fullfile('products',[expname,'_TRM_annual.nc']);
ncid = netcdf.create(ncfname, 'NETCDF4');

%%% Define grid dimensions
x_dimid = netcdf.defDim(ncid, 'x', Nx);
y_dimid = netcdf.defDim(ncid, 'y', Ny);
z_dimid = netcdf.defDim(ncid, 'z', Nr);
t_dimid = netcdf.defDim(ncid, 'time', Nyrs);

%%% Define variables
varIDmodelYrs = netcdf.defVar(ncid, 'model_yrs', 'NC_DOUBLE', [t_dimid]);
varIDcalYrs = netcdf.defVar(ncid, 'cal_yrs', 'NC_DOUBLE', [t_dimid]);

varidXC = netcdf.defVar(ncid, 'XC', 'NC_DOUBLE', [x_dimid, y_dimid]);
varidYC = netcdf.defVar(ncid, 'YC', 'NC_DOUBLE', [x_dimid, y_dimid]);
varidXG = netcdf.defVar(ncid, 'XG', 'NC_DOUBLE', [x_dimid, y_dimid]);
varidYG = netcdf.defVar(ncid, 'YG', 'NC_DOUBLE', [x_dimid, y_dimid]);
varidDXC = netcdf.defVar(ncid, 'DXC', 'NC_DOUBLE', [x_dimid, y_dimid]);
varidDYC = netcdf.defVar(ncid, 'DYC', 'NC_DOUBLE', [x_dimid, y_dimid]);
varidDXG = netcdf.defVar(ncid, 'DXG', 'NC_DOUBLE', [x_dimid, y_dimid]);
varidDYG = netcdf.defVar(ncid, 'DYG', 'NC_DOUBLE', [x_dimid, y_dimid]);
varIDbathy = netcdf.defVar(ncid, 'bathy', 'NC_DOUBLE', [x_dimid, y_dimid]);
varIDshelf = netcdf.defVar(ncid, 'SHELFICEtopo', 'NC_DOUBLE', [x_dimid, y_dimid]);

varidDRF = netcdf.defVar(ncid, 'DRF', 'NC_DOUBLE', [z_dimid]);
varidRC = netcdf.defVar(ncid, 'RC', 'NC_DOUBLE', [z_dimid]);
varidDRC = netcdf.defVar(ncid, 'DRC', 'NC_DOUBLE', [z_dimid]);
varidRF = netcdf.defVar(ncid, 'RF', 'NC_DOUBLE', [z_dimid]);

varidTRMU = netcdf.defVar(ncid, 'TRM_U', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidTRMV = netcdf.defVar(ncid, 'TRM_V', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidTRMW = netcdf.defVar(ncid, 'TRM_W', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidTRMUM = netcdf.defVar(ncid, 'TRM_U_MEAN', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidTRMVM = netcdf.defVar(ncid, 'TRM_V_MEAN', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidTRMWM = netcdf.defVar(ncid, 'TRM_W_MEAN', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidTRMUE = netcdf.defVar(ncid, 'TRM_U_EDDY', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidTRMVE = netcdf.defVar(ncid, 'TRM_V_EDDy', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidTRMWE = netcdf.defVar(ncid, 'TRM_W_EDDY', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidTHETA = netcdf.defVar(ncid, 'TRM_THETA', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);
varidSALT = netcdf.defVar(ncid, 'TRM_SALT', 'NC_DOUBLE', [x_dimid, y_dimid z_dimid t_dimid]);

netcdf.endDef(ncid);

%%% Time grids
netcdf.putVar(ncid, varIDmodelYrs, 0, Nyrs, model_yrs);
netcdf.putVar(ncid, varIDcalYrs, 0, Nyrs, cal_yrs);

%%% For 2D grid fields
ncstart = [0 0];
nccount = [Nx Ny];
netcdf.putVar(ncid, varidXC, ncstart, nccount, XC);
netcdf.putVar(ncid, varidYC, ncstart, nccount, YC);
netcdf.putVar(ncid, varidXG, ncstart, nccount, XG);
netcdf.putVar(ncid, varidYG, ncstart, nccount, YG);
netcdf.putVar(ncid, varidDXC, ncstart, nccount, DXC);
netcdf.putVar(ncid, varidDYC, ncstart, nccount, DYC);
netcdf.putVar(ncid, varidDXG, ncstart, nccount, DXG);
netcdf.putVar(ncid, varidDYG, ncstart, nccount, DYG);
netcdf.putVar(ncid, varIDbathy, ncstart, nccount, bathy);
netcdf.putVar(ncid, varIDshelf, ncstart, nccount, SHELFICEtopo);

%%% Vertical grids
netcdf.putVar(ncid, varidDRF, 0, Nr, DRF);
netcdf.putVar(ncid, varidRC, 0, Nr, RC);
netcdf.putVar(ncid, varidDRC, 0, Nr, DRC(1:Nr));
netcdf.putVar(ncid, varidRF, 0, Nr, RF(1:Nr));

% nccreate(ncfname,'XC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
% ncwrite(ncfname,'XC',XC);
% nccreate(ncfname,'YC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
% ncwrite(ncfname,'YC',YC);
% nccreate(ncfname,'XG','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
% ncwrite(ncfname,'XG',XG);
% nccreate(ncfname,'YG','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
% ncwrite(ncfname,'YG',YG);
% nccreate(ncfname,'DXC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
% ncwrite(ncfname,'DXC',DXC);
% nccreate(ncfname,'DYC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
% ncwrite(ncfname,'DYC',DYC);
% nccreate(ncfname,'DXG','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
% ncwrite(ncfname,'DXG',DXG);
% nccreate(ncfname,'DYG','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
% ncwrite(ncfname,'DYG',DYG);
% nccreate(ncfname,'DRF','Dimensions',{'z',Nr},'FillValue','disable');
% ncwrite(ncfname,'DRF',squeeze(DRF));
% nccreate(ncfname,'RC','Dimensions',{'z',Nr},'FillValue','disable');
% ncwrite(ncfname,'RC',squeeze(RC));
% nccreate(ncfname,'DRC','Dimensions',{'zf',Nr},'FillValue','disable');
% ncwrite(ncfname,'DRC',squeeze(DRC(1:Nr)));
% nccreate(ncfname,'RF','Dimensions',{'zf',Nr+1},'FillValue','disable');
% ncwrite(ncfname,'RF',squeeze(RF(1:Nr)));



%%% Loop over desired years
for m = 1:length(model_yrs)

  %%% Start/end of analysis in model years
  tmin = model_yrs(m) - 1 + 0.02;
  tmax = model_yrs(m) + 0.02;




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

  ncstart = [0 0 0 m-1];
  nccount = [Nx Ny Nr 1];

  netcdf.putVar(ncid, varidTRMU, ncstart, nccount, u_tot);
  netcdf.putVar(ncid, varidTRMV, ncstart, nccount, v_tot);
  netcdf.putVar(ncid, varidTRMW, ncstart, nccount, w_tot);
  netcdf.putVar(ncid, varidTRMUM, ncstart, nccount, u_mean);
  netcdf.putVar(ncid, varidTRMVM, ncstart, nccount, v_mean);
  netcdf.putVar(ncid, varidTRMWM, ncstart, nccount, w_mean);
  netcdf.putVar(ncid, varidTRMUE, ncstart, nccount, u_eddy);
  netcdf.putVar(ncid, varidTRMVE, ncstart, nccount, v_eddy);
  netcdf.putVar(ncid, varidTRMWE, ncstart, nccount, w_eddy);
  netcdf.putVar(ncid, varidTHETA, ncstart, nccount, t_mean);
  netcdf.putVar(ncid, varidSALT, ncstart, nccount, s_mean);

  % nccreate(ncfname,'TRM_U','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'TRM_U',u_tot);
  % nccreate(ncfname,'TRM_V','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'TRM_V',v_tot);
  % nccreate(ncfname,'TRM_W','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'TRM_W',w_tot);
  % nccreate(ncfname,'TRM_U_MEAN','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'TRM_U_MEAN',u_mean);
  % nccreate(ncfname,'TRM_V_MEAN','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'TRM_V_MEAN',v_mean);
  % nccreate(ncfname,'TRM_W_MEAN','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'TRM_W_MEAN',w_mean);
  % nccreate(ncfname,'TRM_U_EDDY','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'TRM_U_EDDY',u_eddy);
  % nccreate(ncfname,'TRM_V_EDDY','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'TRM_V_EDDY',v_eddy);
  % nccreate(ncfname,'TRM_W_EDDY','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'TRM_W_EDDY',w_eddy);
  % nccreate(ncfname,'THETA','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'THETA',t_mean);
  % nccreate(ncfname,'SALT','Dimensions',{'x',Nx,'y',Ny,'z',Nr},'FillValue','disable');
  % ncwrite(ncfname,'SALT',s_mean);

end

%%% Close up the netcdf file
netcdf.close(ncid);