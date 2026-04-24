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
for n = 1:Nexps
  expname = T.Name(n);
  expname = expname{1};
  tmax = T.EndTime_yr_(n) + 0.05;
  tmin = tmax - 3; %%% Last 3 years

  %%% If we can't overwrite output files then we need to check whether they
  %%% already exist
  if (~overwrite)

    %%% Construct output file name to check whether it exists
    outfname = [expname,'_MOC_PD0'];
    if (calc_psi_eddy)
      if (use_layers)
        estr = '_layers';
      else
        estr = '_TRM';
      end
    else
      estr = '_noeddy';
    end
    outfname = [outfname,estr];
    if (use_PsiBT)
      outfname = [outfname,'_PsiBT'];
    else
      if (deform_cavity)
        outfname = [outfname,'_deform'];
      elseif (gl_coord)
        outfname = [outfname,'_GLcoord'];
      end
    end
    outfname = [outfname,'.mat'];

    if (exist(fullfile(outdir,outfname)))
      continue;
    end

  end

  calcOverturning (expdir,expname,outdir,tmin,tmax,calc_psi_eddy,deform_cavity,use_PsiBT,use_layers,gl_coord)
end

