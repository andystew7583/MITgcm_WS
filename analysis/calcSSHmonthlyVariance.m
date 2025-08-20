%%%
%%% calcSSHmonthlyVariance.m
%%%
%%% Calculates variance of SSH using instantaneous model output, averaged
%%% over months.
%%%


%%% This needs to be set to ensure we are using the correct output
%%% frequency
loadexp;
diagfreq = diag_frequency(end);

%%% High-frequency output
dumpFreq = abs(diag_frequency(69));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Monthly output
dumpFreq_mon = abs(diag_frequency(1));
nDumps_mon = round(endTime/dumpFreq_mon);
dumpIters_mon = round((1:nDumps_mon)*dumpFreq_mon/deltaT);
dumpIters_mon = dumpIters_mon(dumpIters_mon > nIter0);
nDumps_mon = length(dumpIters_mon)

%%% To store the result
ssh_mon_var = zeros(Nx,Ny,nDumps_mon);
ssh_mon_mean = zeros(Nx,Ny,nDumps_mon);

%%% Times corresponding to the middle of each month
times_brk = [0 dumpIters_mon*deltaT];
Navg_mon = zeros(1,nDumps_mon);
times = 0.5*(times_brk(1:end-1)+times_brk(2:end));

%%% Pre-compute number of snapshots in each averaging window
for n=1:nDumps

  tt(n) =  dumpIters(n)*deltaT;
  monidx = find(tt(n)<=times_brk,1,'first')-1;
  Navg_mon(monidx) = Navg_mon(monidx) + 1;
 
end


%%% Loop over outputs
for n=1:nDumps

  tt(n) =  dumpIters(n)*deltaT;
  tt(n)

  %%% Load SSH data
  eta = rdmdsWrapper(fullfile(exppath,'/results/ETAN_inst'),dumpIters(n));  
  if (isempty(eta))
    error(['No data at n = ',num2str(n)]);
  end

  monidx = find(tt(n)<=times_brk,1,'first')-1;
  ssh_mon_mean(:,:,monidx) = ssh_mon_mean(:,:,monidx) + eta/Navg_mon(monidx);
  ssh_mon_var(:,:,monidx) = ssh_mon_var(:,:,monidx) + eta.^2/Navg_mon(monidx);
 
end

%%% Subtract squared means
ssh_mon_var = ssh_mon_var - ssh_mon_mean.^2;

%%% Write to data file
ncfname = fullfile('products',[expname,'_SSHmonVar.nc']);
nccreate(ncfname,'XC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'XC',XC);
nccreate(ncfname,'YC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'YC',YC);
nccreate(ncfname,'bathy','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'bathy',bathy);
nccreate(ncfname,'SHELFICEtopo','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'SHELFICEtopo',SHELFICEtopo);
nccreate(ncfname,'ssh_mon_mean','Dimensions',{'x',Nx,'y',Ny,'t',nDumps_mon},'FillValue','disable');
ncwrite(ncfname,'ssh_mon_mean',ssh_mon_mean);
nccreate(ncfname,'ssh_mon_var','Dimensions',{'x',Nx,'y',Ny,'t',nDumps_mon},'FillValue','disable');
ncwrite(ncfname,'ssh_mon_var',ssh_mon_var);
