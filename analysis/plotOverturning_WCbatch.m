%%%
%%% plotOverturnin_WCbatch.m
%%%
%%% Plots FRIS melt rate and continental shelf buoyancy loss for all experiments in our Weddell
%%% Catastrophe batch.
%%%

%%% Pointer to experiment directory
expdir = '../experiments';

%%% Pointer to storage directory for output .mat files
proddir = './products_WCbatch';

%%% Set true to overwrite existing files
overwrite = false;


%%% Read Google Drive spreadsheet with experiment data
url = "https://docs.google.com/spreadsheets/d/19mLumSIpCtc6AsPfuzLBkuPxVXoyXuaD3yWbCmxe8FY/export?format=csv&gid=0";
T = readtable(url);

%%% Narrow down to production experiments
T = T(T.production_=="Y",:);

%%% Narrow down to experiments that are complete and downloaded
T = T(T.completed_=="Y" & T.downloaded_=="Y",:);

%%% Restrict to a specific batch of experiments
T = T(T.batch >= 0 & T.batch<=12,:);

%%% Constant physical parameters
beta = 8e-4;
g = 9.81;
f0 = 2*(2*pi*366/365/86400)*sind(65);
Sref = 34.67;
rhofresh = 1000;

%%% Pre-allocate storage
strat = zeros(1,size(T,1));
batchnum = zeros(1,size(T,1));
dpyc = zeros(1,size(T,1));
F_in = zeros(1,size(T,1));
Psimax_neg = zeros(1,size(T,1));
Psimax_pos = zeros(1,size(T,1));
Psimax_cavity = zeros(1,size(T,1));
has_WC = zeros(1,size(T,1));
has_pos_MOC = zeros(1,size(T,1));

%%% Loop through experiments and compute MOC
Nexps = size(T,1);
for n = 1:Nexps
  expname = T.Name(n);
  expname = expname{1};
  tmax = T.EndTime_yr_(n) + 0.05;
  tmin = tmax - 9; %%% Last 3 years

  %%% Record upstream stratification and pycnocline depth
  strat(n) = T.strat(n);
  dpyc(n) = 400+T.Dpyc(n);
  batchnum(n) = T.batch(n);    
  if (T.WeddellCatastrophe_(n)=="Y")
    has_WC(n) = 1;
  elseif (T.WeddellCatastrophe_(n)=="N")
    has_WC(n) = 0;
  else
    has_WC(n) = -1;
  end

  %%% Load pre-computed melt rates
  datafname = [expname,'_MOC_PD0_layers_GLcoord.mat'];
  load(fullfile(proddir,datafname));
  tt = tt / t1year;
  psi_tot_tavg = mean(psi_tot,3)/1e6;

  jidx_cavity = find(eta==0);
  jidx_shelfbreak = find(eta==3.5);
  Psimax_neg(n) = -min(psi_tot_tavg(jidx_shelfbreak,:,:));
  Psimax_pos(n) = max(psi_tot_tavg(jidx_shelfbreak,:,:));
  kidx_neg = find(psi_tot_tavg(jidx_shelfbreak,:)==-Psimax_neg(n));
  kidx_pos = find(psi_tot_tavg(jidx_shelfbreak,:)==Psimax_pos(n));
  Psimax_cavity(n) = max(psi_tot_tavg(jidx_cavity,:,:));
  has_pos_MOC(n) = Psimax_pos(n) > 0.5*Psimax_cavity(n);




end


Sstrat = strat./(g*beta);
FWupstream = (1/32)*(g.*beta.*dpyc.^4.*Sstrat.^2)/(12*f0*Sref) * rhofresh * t1year / 1e12;