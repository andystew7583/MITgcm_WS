%%%
%%% syntheticMooring.m
%%%
%%% Creates synthetic mooring data.
%%%

%%% Mooring location and output file name
% load('SealComparison/mooringA253_337.mat');
% outfname = 'modelA253.mat';
load('SealComparison/mooringA254_337.mat');
outfname = 'modelA254.mat';

%%% Read experiment data
setExpname; %%% To set directory paths
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2'; %%% Enforce high-resolution experiment
loadexp;
depth = -squeeze(zz);

%%% Frequency of diagnostic output
dumpStart = 1578240;
dumpStep = 43200/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;

%%% Find index of closest location to this mooring
im = find(XC(:,1)<mooring_lat_lon(2),1,'last');
jm = find(YC(1,:)<mooring_lat_lon(1),1,'last');
if ((XC(im+1,1)-mooring_lat_lon(2)) < (mooring_lat_lon(2)-XC(im,1)))
  im = im+1;
end
if ((YC(1,jm+1)-mooring_lat_lon(1)) < (mooring_lat_lon(1)-YC(1,jm)))
  jm = jm+1;
end
  
%%% Loop through iterations
uvel = zeros(length(dumpIters),Nr);
vvel = zeros(length(dumpIters),Nr);
salt = zeros(length(dumpIters),Nr);
theta = zeros(length(dumpIters),Nr);
time = zeros(length(dumpIters),1);
for n=1:length(dumpIters)
  
  %%% Load instantaneous model state
  u_model = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),dumpIters(n));      
  v_model = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),dumpIters(n));
  t_model = rdmdsWrapper(fullfile(exppath,'/results/THETA_12hourly'),dumpIters(n));      
  s_model = rdmdsWrapper(fullfile(exppath,'/results/SALT_12hourly'),dumpIters(n));
  if (isempty(u_model) || isempty(v_model) || isempty(s_model) || isempty(t_model))   
    break;
  end
  
  %%% Extract data at mooring location
  uvel(n,:) = squeeze(u_model(im,jm,:));
  vvel(n,:) = squeeze(v_model(im,jm,:));
  theta(n,:) = squeeze(t_model(im,jm,:));
  salt(n,:) = squeeze(s_model(im,jm,:));
  time(n) = datenum('2008-01-01')+dumpIters(n)*deltaT/86400;
    
end

%%% Write product to output file
save(fullfile('SealComparison',outfname), ...
  'uvel','vvel','salt','theta','time','depth' ...
  );