%%%
%%% calcSubmesStatistics.m
%%%
%%% Computes statistics relevant to submesoscale motions.
%%%

%%% Need scripts from analysis base directory
addpath ..;

%%% Read experiment data
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
expdir = '../experiments/';
loadexp;

%%% Vertical grid spacing matrix
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

%%% Diagnostic index corresponding to instantaneous velocity
% diagnum = length(diag_frequency);
diagnum = 1;
% diagnum = 69;

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(diagnum);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round(nIter0 + (1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
nDumps = length(dumpIters);

%%% High-frequency 1/24 degree run
dumpStart = 1578240;
dumpStep = 43200/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;

%%% Quasi-meridional grid
deform_cavity = false;
gl_coord = false;
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity,gl_coord);

%%% Mask region of interest - south of shelf break, outside cavities
msk = (hFacC(:,:,1)>0) & (ETA < 3);

%%% Physical parameters
Omega = 2*pi*366/365/86400;
ff = 2*Omega*sind(YG);
g = 9.81;

%%% Grid for storing vorticity/divergence statistics
binedges = [-3:0.01:3];
binmid = 0.5*(binedges(1:end-1)+binedges(2:end));
Nbins = length(binmid);
vorthist = zeros(nDumps,Nbins);
divhist = zeros(nDumps,Nbins);

%%% Loop through iterations
for n=1:length(dumpIters)
 
  tt(n) =  dumpIters(n)*deltaT/86400;
  tt(n)
  
  %%% Attempt to load either instantaneous velocities or their squares
  zlev = 1;
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),dumpIters(n)) ;      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),dumpIters(n));
  uvel = uvel(:,:,zlev);
  vvel = vvel(:,:,zlev);
  uvel(hFacW(:,:,zlev)==0) = NaN;
  vvel(hFacS(:,:,zlev)==0) = NaN;

  %%% Vorticity on a z-level
  vort = zeros(Nx,Ny);
  vort(:,2:Ny) = - (uvel(:,2:Ny)-uvel(:,1:Ny-1))./DYC(:,2:Ny);
  vort(2:Nx,:) = vort(2:Nx,:) + (vvel(2:Nx,:)-vvel(1:Nx-1,:))./DXC(2:Nx,:);
  vort = vort ./ abs(ff);

  %%% Divergence on a z-level
  div = zeros(Nx,Ny);
  zlev = 1;
  div(:,1:Ny-1) = (vvel(:,2:Ny,zlev)-vvel(:,1:Ny-1,zlev))./DYG(:,1:Ny-1);
  div(1:Nx-1,:) = div(2:Nx,:) + (uvel(2:Nx,:,zlev)-uvel(1:Nx-1,:,zlev))./DXG(1:Nx-1,:);
  div = div ./ abs(ff);

  %%% Compute histograms for this time slice
  vorthist(n,:) = histcounts(vort(msk),binedges);
  divhist(n,:) = histcounts(div(msk),binedges);

end

%%% Write to output file
outfname = [expname,'_SubmesStats.mat'];
save(fullfile('products',outfname), ...
    'tt','binedges','binmid',...
    'vorthist','divhist',...    
    '-v7.3');