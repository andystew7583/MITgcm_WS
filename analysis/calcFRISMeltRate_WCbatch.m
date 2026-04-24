%%%
%%% calcFRISMeltRate_WCbatch.m
%%%
%%% Computes FRIS melt rate and continental shelf buoyancy loss for all experiments in our Weddell
%%% Catastrophe batch.
%%%

%%% Pointer to experiment directory
expdir = '../experiments';

%%% Pointer to storage directory for output .mat files
outdir = './products_WCbatch';

%%% Set true to overwrite existing files
overwrite = false;

%%% Set true to compute effective melt rate (heat flux into the cavity)
compute_eff_melt = false;

%%% Set true to compute upstream freshwater flux
compute_upstream = false;



%%% Read Google Drive spreadsheet with experiment data
url = "https://docs.google.com/spreadsheets/d/19mLumSIpCtc6AsPfuzLBkuPxVXoyXuaD3yWbCmxe8FY/export?format=csv&gid=0";
T = readtable(url);

%%% Narrow down to production experiments
T = T(T.production_=="Y",:);

%%% Narrow down to experiments that are complete and downloaded
T = T(T.completed_=="Y" & T.downloaded_=="Y",:);

%%% Loop through experiments and compute MOC
Nexps = size(T,1);
for n = 1:Nexps
  expname = T.Name(n);
  expname = expname{1};
  tmax = T.EndTime_yr_(n) + 0.05;
  tmin = tmax - 3; %%% Last 3 years

  %%% Keep track of progress
  disp(expname);

  %%% Output data .mat file name
  outfname = [expname,'_FRISMeltRate.mat'];

  %%% If we can't overwrite output files then we need to check whether they
  %%% already exist
  if (~overwrite)

    %%% Construct output file name to check whether it exists    
    if (exist(fullfile(outdir,outfname)))
      disp('Skipping ...');
      continue;
    end

  end

  %%% Compute melt rates  
  [tt,SHImelt,SIprod,SHImelt_mean, ...
    SIprod_mean,XC,YC,bathy,SHELFICEtopo, ...
    SHImelt_eff,SHImelt_tend,SHImelt_diff, ...
    SHImelt_tflux,SHImelt_conv,FWupstream] ...
    = calcFRISMeltTimeSeries (expdir,expname,compute_eff_melt,compute_upstream,tmax,tmax-tmin);

  %%% Save to .matfile
  save(fullfile(outdir,outfname),...
    'tt','SHImelt','SIprod','SHImelt_mean', ...
    'SIprod_mean','XC','YC','bathy','SHELFICEtopo', ...
    'SHImelt_eff','SHImelt_tend','SHImelt_diff', ...
    'SHImelt_tflux','SHImelt_conv','FWupstream');

end

