%%%
%%% calcOverturning_WCbatch.m
%%%
%%% Computes overturning circulation for all experiments in our Weddell
%%% Catastrophe batch.
%%%

%%% Pointer to experiment directory
expdir = '../experiments';

%%% Pointer to storage directory for output .mat files
outdir = './products_WCbatch';

%%% MOC calculation options: 
calc_psi_eddy = false; %%% just compute mean MOC
deform_cavity = false; %%% Don't use coordinate system centered on Ronne polynya
use_PsiBT = false; %%% Don't use barotropic streamfunction coordinats
use_layers = false; %%% Don't use LAYERS diagnostics 
gl_coord = true; %%% Do use coordinate system aligned with ice front

%%% Set true to overwrite existing files
overwrite = false;


%%% Read Google Drive spreadsheet with experiment data
url = "https://docs.google.com/spreadsheets/d/19mLumSIpCtc6AsPfuzLBkuPxVXoyXuaD3yWbCmxe8FY/export?format=csv&gid=0";
T = readtable(url);

%%% Narrow down to production experiments
T = T(T.production_=="Y",:);

%%% Narrow down to experiments that are complete and downloaded
T = T(T.completed_=="Y" & T.downloaded_=="Y",:);

%%% Loop through experiments and compute MOC
Nexps = size(T,1);
for n = 20:Nexps
  expname = T.Name(n);
  expname = expname{1};
  tmax = T.EndTime_yr_(n) + 0.05;
  tmin = tmax - 3; %%% Last 3 years
  calcOverturning (expdir,expname,outdir,tmin,tmax,calc_psi_eddy,deform_cavity,use_PsiBT,use_layers,gl_coord)
end

