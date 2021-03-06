%%%
%%% regridOBCS.m
%%%
%%% Regrids output from a lower-resolution run to serve as boundary
%%% conditions for a higher-resolution run.
%%%

%%% Matlab utilities 
addpath ../newexp_utils
addpath ../utils/matlab
addpath ../analysis

%%% Load core simulation parameters
defineGrid;

%%% Defines date range over which to extract data - N.B. need to ensure
%%% overlaps so that entire hi-res simulation period has boundary data
startdate_obcs = datenum([num2str(start_year-1),'-12-15']);
enddate_obcs = datenum([num2str(endyr+1),'-02-15']);

%%% Load low-resolution grid
expdir_lo = '../experiments/';
expname_lo = 'hires_seq_onetwelfth_RTOPO2';
resultspath_lo = fullfile(expdir_lo,expname_lo,'results');
XC_lo = rdmds(fullfile(resultspath_lo,'XC'));
YC_lo = rdmds(fullfile(resultspath_lo,'YC'));
RC_lo = rdmds(fullfile(resultspath_lo,'RC'));
hFacC_lo = rdmds(fullfile(resultspath_lo,'hFacC'));

%%% Temporal parameters for low-res run
startdate_lo = datenum('2007-01-01');
diag_frequency_lo = 2.62960000e+06;
endTime_lo = 283996920;
deltaT_lo = 1.20000000e+02;
nIter0_lo = 1;
dumpFreq_lo = abs(diag_frequency_lo);
nDumps_lo = round(endTime_lo/dumpFreq_lo);
dumpIters_lo = round((1:nDumps_lo)*dumpFreq_lo/deltaT_lo);
dumpIters_lo = dumpIters_lo(dumpIters_lo > nIter0_lo);
dumpTimes_lo = startdate_lo + dumpIters_lo*deltaT_lo/t1day;
idx_obcs = find((dumpTimes_lo>startdate_obcs) & (dumpTimes_lo<enddate_obcs));
Nt_obcs = length(idx_obcs);

%%% Find slices from which to extract data from low-res run
yidx_N = find(YC_lo(1,:)>ymax,1,'first');
xidx_E = find(XC_lo(:,1)>xmax,1,'first');
landidx_N = find(hFacC_lo(:,yidx_N,1)==0);
landidx_E = find(hFacC_lo(xidx_E,:,1)==0);

%%% Interpolation grids
[XX_xz,ZZ_xz] = meshgrid(xmc,mrc);
[YY_yz,ZZ_yz] = meshgrid(ymc,mrc);
[XX_xz_lo,ZZ_xz_lo] = meshgrid(XC_lo(:,1),RC_lo);
[YY_yz_lo,ZZ_yz_lo] = meshgrid(YC_lo(1,:),RC_lo);

%%% Lists of 2D and 1D OB variables and file names on each boundary
vars_2D = {'THETA','SALT','UVEL','VVEL'};
vars_1D = {'ETAN','SIuice','SIvice','SIarea','SIheff','SIhsnow','SIhsalt'};
obcsFilesN_2D = {OBNtFile,OBNsFile,OBNuFile,OBNvFile};
obcsFilesE_2D = {OBEtFile,OBEsFile,OBEuFile,OBEvFile};
obcsFilesN_1D = {OBNetaFile,OBNuiceFile,OBNviceFile,OBNaFile,OBNhFile,OBNsnFile,OBNslFile};
obcsFilesE_1D = {OBEetaFile,OBEuiceFile,OBEviceFile,OBEaFile,OBEhFile,OBEsnFile,OBEslFile};

%%% Loop over 2D OB variables
for m = 1:length(vars_2D)
  
  disp(vars_2D{m});
  
  %%% Initialize 3D OB arrays
  OBNarray = zeros(Nx,Nr,Nt_obcs);
  OBEarray = zeros(Ny,Nr,Nt_obcs);
  
  %%% Loop through snapshots and interpolate to hi-res model grid
  for n = 1:Nt_obcs
    
    disp(n);
    
    tmp = rdmdsWrapper(fullfile(resultspath_lo,vars_2D{m}),dumpIters_lo(idx_obcs(n)));
    tmpN = squeeze(tmp(:,yidx_N,:));
    tmpE = squeeze(tmp(xidx_E,:,:));
    OBNarray(:,:,n) = interpBdyData2(tmpN,XX_xz_lo,ZZ_xz_lo,XX_xz,ZZ_xz);    
    OBEarray(:,:,n) = interpBdyData2(tmpE,YY_yz_lo,ZZ_yz_lo,YY_yz,ZZ_yz);    
    
  end
  
  %%% Write to OB files
  writeDataset(OBNarray,fullfile(inputconfigdir,obcsFilesN_2D{m}),ieee,prec);
  writeDataset(OBEarray,fullfile(inputconfigdir,obcsFilesE_2D{m}),ieee,prec);
  
end

%%% Loop over 1D OB variables 
for m = 1:length(vars_1D)
  
  disp(vars_1D{m});
  
  %%% Initialize 2D OB arrays
  OBNarray = zeros(Nx,Nt_obcs);
  OBEarray = zeros(Ny,Nt_obcs);
  
  %%% Loop through snapshots and interpolate to hi-res model grid
  for n = 1:Nt_obcs
    
    disp(n);
    
    tmp = rdmdsWrapper(fullfile(resultspath_lo,vars_1D{m}),dumpIters_lo(idx_obcs(n)));
    tmpN = squeeze(tmp(:,yidx_N));
    tmpE = squeeze(tmp(xidx_E,:));
    OBNarray(:,n) = interpBdyData1(tmpN,XC_lo(:,1),xmc,landidx_N);    
    OBEarray(:,n) = interpBdyData1(tmpE,YC_lo(1,:),ymc,landidx_E);    
    
  end
  
  %%% Write to OB files
  writeDataset(OBNarray,fullfile(inputconfigdir,obcsFilesN_1D{m}),ieee,prec);
  writeDataset(OBEarray,fullfile(inputconfigdir,obcsFilesE_1D{m}),ieee,prec);
  
end




%%%
%%% interpBdyData2
%%%
%%% Convenience function to interpolate boundary data onto the model grid;
%%%
function data_interp = interpBdyData2 (data,XD,ZD,XI,ZI)
 
  
  %%% Remove land points
  data(data==0) = NaN;
  
  %%% To store interpolated data
  data_interp = zeros(size(XI,2),size(XI,1));
  
  
  data_interp(:,:) = interp2(XD,ZD,data(:,:)',XI,ZI,'linear')';
  data_interp(:,:) = inpaint_nans(data_interp(:,:),4); %%% Crude extrapolation

end


%%%
%%% interpBdyData1
%%%
%%% Convenience function to interpolate boundary data onto the model grid.
%%%
function data_interp = interpBdyData1 (data,xd,xi,landidx)
 
  %%% Remove land points
  xd(landidx) = [];
  
  %%% To store interpolated data
  data_interp = 0*xi;
  
  %%% Do interpolation
  data_tmp = data(:);
  data_tmp(landidx) = [];
  data_interp(:) = interp1(xd,data_tmp,xi,'nearest','extrap')';
  
end