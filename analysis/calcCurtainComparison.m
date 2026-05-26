%%%
%%% calcCurtainComparison.m
%%%
%%% Computes diagnostics to compare runs with and without a "curtain"
%%% across the Filchner Trough.
%%%

%%% Pointer to experiment directory
expdir = '../experiments';

%%% Pointer to storage directory for output .mat files
proddir = './products_WCbatch';

%%% Set true to overwrite existing files
overwrite = false;

%%% Seconds in 1 year, averaged over the 9-year simulation period
t1year = 31555200;

%%% Read Google Drive spreadsheet with experiment data
url = "https://docs.google.com/spreadsheets/d/19mLumSIpCtc6AsPfuzLBkuPxVXoyXuaD3yWbCmxe8FY/export?format=csv&gid=0";
T = readtable(url);

%%% List of experiments
expnames = {...
  'WC_onethird_ref',...
  'WC_onethird_strat1e-4',...
  'WC_onethird_C450',...
  'WC_onethird_strat1e-4_C450'};

%%% Descriptive names/labels for experiments
titles = {...
  'Reference case', ...
  'Upstream freshwater pert.', ...
  'Filchner Trough curtain', ...
  'Freshwater + curtain'};

%%% To store output structures
diags_curtain = cell(1,4);

%%% Narrow down to production experiments
T = T(T.production_=="Y",:);

%%% Narrow down to experiments that are complete and downloaded
T = T(T.completed_=="Y" & T.downloaded_=="Y",:);

%%% Pre-allocate storage
strat = zeros(1,size(T,1));
batchnum = zeros(1,size(T,1));
dpyc = zeros(1,size(T,1));
FRISmelt_batch = zeros(1,size(T,1));
initFRISmelt_batch = zeros(1,size(T,1));
SIprod_batch = zeros(1,size(T,1));
has_WC = zeros(1,size(T,1));

%%% Loop through experiments and compute MOC
Nexps = size(T,1);
for m = 1:length(diags_curtain)

  %%% Locate experiment parameters
  expname = expnames{m};  
  n = -1;
  for i = 1:length(T.Name)
    if (strcmp(T.Name{i},expname))
      n = i;
    end
  end
  if (n == -1)
    error(['Could not find experiment: ',expname]);
  end

  %%% Extract time range for analysis
  tmax = T.EndTime_yr_(n) + 0.05;
  tmin = tmax - 9; %%% Last 3 years

  %%% Load pre-computed melt rates
  datafname = [expname,'_FRISMeltRate.mat'];
  load(fullfile(proddir,datafname));
  tt = tt / t1year;

  %%% Store melt rate time series
  diags_curtain{m}.FRISmelt = -SHImelt*t1year/1e12;
  diags_curtain{m}.tt = tt;

  %%% Load simulation parameters
  loadexp;

  %%% Store grids
  diags_curtain{m}.bathy = bathy;
  diags_curtain{m}.SHELFICEtopo = SHELFICEtopo;
  diags_curtain{m}.XC = XC;
  diags_curtain{m}.YC = YC;
  diags_curtain{m}.RC = RC;
  diags_curtain{m}.RF = RF;
  diags_curtain{m}.DRF = DRF;
  diags_curtain{m}.Nx = Nx;
  diags_curtain{m}.Ny = Ny;
  diags_curtain{m}.Nr = Nr;
  diags_curtain{m}.hFacC = hFacC;

  %%% Compute iterations at which model output is written
  dumpFreq = abs(diag_frequency(4));
  nDumps = round(endTime/dumpFreq);
  dumpIters = round((1:nDumps)*dumpFreq/deltaT);
  dumpIters = dumpIters(dumpIters > nIter0);

  %%% Read time-mean model hydrography
  diags_curtain{m}.theta_tavg = readIters(exppath,'THETA',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);
  diags_curtain{m}.salt_tavg = readIters(exppath,'SALT',dumpIters,deltaT,tmin*t1year,tmax*t1year,Nx,Ny,Nr);


end

%%% Save to .mat file
save(fullfile(proddir,'curtain_diags.mat'),'expnames','titles','diags_curtain');
