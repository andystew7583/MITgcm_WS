%%%
%%% calcSSHvariance.m
%%%
%%% Calculates variance of SSH using instantaneous model output.
%%%


%%% This needs to be set to ensure we are using the correct output
%%% frequency
loadexp;
diagfreq = diag_frequency(end);

% %%% For daily/12-hourly outputs
% dumpStart = 1578240;
% dumpStep = 86400/2/60;
% nDumps = 731;
% dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;

%%% For daily/12-hourly outputs, SSH-specific re-run
dumpStart = 720;
dumpStep = 86400/2/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;

%%% To store the result
ssh = zeros(Nx,Ny);
ssh_var = zeros(Nx,Ny);
ish = zeros(Nx,Ny);
ish_var = zeros(Nx,Ny);

L_runmean = 10;
ssh_runmean = zeros(Nx,Ny,nDumps);
ssh_runpert = zeros(Nx,Ny,nDumps);
ssh_var_runmean = zeros(Nx,Ny,nDumps);
ssh_runpert_var = zeros(Nx,Ny);
ish_runmean = zeros(Nx,Ny,nDumps);
ish_runpert = zeros(Nx,Ny,nDumps);
ish_var_runmean = zeros(Nx,Ny,nDumps);
ish_runpert_var = zeros(Nx,Ny);
N_runmean = zeros(1,nDumps);

%%% Pre-compute lengths of running means
for n=1:nDumps
  runidx = max(n-L_runmean,1):1:min(n+L_runmean,nDumps);
  N_runmean(n) = length(runidx);
end

%%% Loop over outputs
for n=1:nDumps

  tt(n) =  dumpIters(n)*deltaT;
  tt(n)

  %%% Attempt to load melt ave per month
  SIheff = rdmdsWrapper(fullfile(exppath,'/results/SIheff_12hourly'),dumpIters(n));  
  SIhsnow = rdmdsWrapper(fullfile(exppath,'/results/SIhsnow_12hourly'),dumpIters(n));  
  SIarea = rdmdsWrapper(fullfile(exppath,'/results/SIarea_12hourly'),dumpIters(n));  
  eta = rdmdsWrapper(fullfile(exppath,'/results/ETAN_12hourly'),dumpIters(n));  
  if (isempty(eta))
    error(['No data at n = ',num2str(n)]);
  end

n

  %%% Ice surface height
  eta_i = eta;
  idx = find(SIarea > 0);
  eta_i(idx) = eta_i(idx) + (SIhsnow(idx) + SIheff(idx))./SIarea(idx);
  
  %%% Add to averages
  ish = ish + eta_i/nDumps;
  ish_var = ish_var + eta_i.^2/nDumps;
  ssh = ssh + eta/nDumps;
  ssh_var = ssh_var + eta.^2/nDumps;
  
  %%% Add to running means and variances
  runidx = max(n-L_runmean,1):1:min(n+L_runmean,nDumps);
  N_runmean(n) = length(runidx);
  for m = runidx
    ssh_runmean(:,:,m) = ssh_runmean(:,:,m) + eta/N_runmean(m);
    ssh_var_runmean(:,:,m) = ssh_var_runmean(:,:,m) + eta.^2/N_runmean(m);
    ish_runmean(:,:,m) = ish_runmean(:,:,m) + eta_i/N_runmean(m);    
    ish_var_runmean(:,:,m) = ish_var_runmean(:,:,m) + eta_i.^2/N_runmean(m);
  end
  ssh_runpert(:,:,m) = eta;
  ish_runpert(:,:,m) = eta_i;

end

%%% Compute perturbations relative to running means
ssh_runpert = ssh_runpert - ssh_runmean;
ish_runpert = ish_runpert - ish_runmean;

%%% Subtract squared means
ssh_var = ssh_var - ssh.^2;
ish_var = ish_var - ish.^2;
ssh_var_runmean = ssh_var_runmean - ssh_runmean.^2;
ish_var_runmean = ish_var_runmean - ish_runmean.^2;

%%% Write to data file
ncfname = fullfile('products',[expname,'_SSHvar.nc']);
nccreate(ncfname,'XC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'XC',XC);
nccreate(ncfname,'YC','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'YC',YC);
nccreate(ncfname,'bathy','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'bathy',bathy);
nccreate(ncfname,'SHELFICEtopo','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'SHELFICEtopo',SHELFICEtopo);
nccreate(ncfname,'ssh','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'ssh',ssh);
nccreate(ncfname,'ssh_var','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'ssh_var',ssh_var);
nccreate(ncfname,'ish','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'ish',ish);
nccreate(ncfname,'ish_var','Dimensions',{'x',Nx,'y',Ny},'FillValue','disable');
ncwrite(ncfname,'ish_var',ish_var);
nccreate(ncfname,'ssh_runmean','Dimensions',{'x',Nx,'y',Ny,'t',nDumps},'FillValue','disable');
ncwrite(ncfname,'ssh_runmean',ssh_runmean);
nccreate(ncfname,'ssh_runpert','Dimensions',{'x',Nx,'y',Ny,'t',nDumps},'FillValue','disable');
ncwrite(ncfname,'ssh_runpert',ssh_runpert);
nccreate(ncfname,'ssh_var_runmean','Dimensions',{'x',Nx,'y',Ny,'t',nDumps},'FillValue','disable');
ncwrite(ncfname,'ssh_var_runmean',ssh_var_runmean);
nccreate(ncfname,'ish_runmean','Dimensions',{'x',Nx,'y',Ny,'t',nDumps},'FillValue','disable');
ncwrite(ncfname,'ish_runmean',ish_runmean);
nccreate(ncfname,'ish_runpert','Dimensions',{'x',Nx,'y',Ny,'t',nDumps},'FillValue','disable');
ncwrite(ncfname,'ish_runpert',ish_runpert);
nccreate(ncfname,'ish_var_runmean','Dimensions',{'x',Nx,'y',Ny,'t',nDumps},'FillValue','disable');
ncwrite(ncfname,'ish_var_runmean',ish_var_runmean);

%%% Make plots

bathy_plot = bathy;
bathy_plot(sum(hFacC,3)==0) = NaN;
figure(1);
pcolor(XC,YC,sqrt(abs(ssh_var)));
hold on;
[C,h] = contour(XC,YC,-bathy_plot,[0 250 500 1000 2000 3000 4000],'EdgeColor',[.3 .3 .3]); 
clabel(C,h);
shading interp;
colormap(cmocean('amp',16));
colorbar;
caxis([0 0.08]);
xlabel('Longitude');
ylabel('Latitude');
title('SSH standard deviation (m)');
set(gca,'FontSize',14);


figure(2);
pcolor(XC,YC,sqrt(abs(ish_var)));
shading interp;
colormap(cmocean('amp',50));
colorbar;
caxis([0 1]);

figure(4);
pcolor(XC,YC,eta_i);
shading interp;
colormap(cmocean('amp',50));
colorbar;
caxis([0 5]);

xidx = find(XC(:,1)>-40 & XC(:,1)<-30);
yidx = find(YC(1,:)>-75 & YC(1,:)<-74);
ssh_var_runmean_avg = squeeze(sum(sum(ssh_var_runmean(xidx,yidx,:).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);
ish_var_runmean_avg = squeeze(sum(sum(ish_var_runmean(xidx,yidx,:).*DXG(xidx,yidx).*DYG(xidx,yidx),1),2)) / sum(sum(DXG(xidx,yidx).*DYG(xidx,yidx),1),2);

figure(5);
plot(tt/t1year,sqrt(ssh_var_runmean_avg));

figure(6);
plot(tt/t1year,sqrt(ish_var_runmean_avg));
