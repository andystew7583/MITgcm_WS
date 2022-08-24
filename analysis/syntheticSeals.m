%%%
%%% syntheticSeals.m
%%%
%%% Creates synthetic seal measurements data.
%%%

%%% Seal tags to process
seal_numbers = [21126,40020];
Nseals = length(seal_numbers);
sealdir = 'SealComparison';

%%% Extract seal measurement positions/times
seal_lons = cell(1,Nseals);
seal_lats = cell(1,Nseals);
seal_times = cell(1,Nseals);
seal_yeardays = cell(1,Nseals);
for m = 1:Nseals
  load(fullfile(sealdir,['seal',num2str(seal_numbers(m)),'.mat']));
  seal_times{m} = tag_datenum';
  yeardays = 0*tag_datenum;
  for j=1:length(yeardays)
    yeardays(j) = tag_datenum(j) - datenum([datestr(tag_datenum(j),'yyyy'),'-01-01']);
  end
  seal_yeardays{m} = yeardays; 
  seal_lons{m} = tag_lat_lon(:,2);
  seal_lats{m} = tag_lat_lon(:,1);
end

%%% Read experiment data
setExpname; %%% To set directory paths
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2'; %%% Enforce high-resolution experiment
loadexp;
model_depth = -squeeze(zz);
model_lons = XC(:,1); %%% Regular lat/lon grid so we can do this
model_lats = YC(1,:);

%%% Frequency of diagnostic output
dumpStart = 1578240;
dumpStep = 43200/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;
 
%%% To store model state at seal measurement locations
seal_uvel = cell(1,Nseals);
seal_vvel = cell(1,Nseals);
seal_salt = cell(1,Nseals);
seal_theta = cell(1,Nseals);
for m = 1:Nseals
  seal_uvel{m} = zeros(length(seal_times{m}),Nr);
  seal_vvel{m} = zeros(length(seal_times{m}),Nr);
  seal_salt{m} = zeros(length(seal_times{m}),Nr);
  seal_theta{m} = zeros(length(seal_times{m}),Nr);
end

%%% Loop through time outputs from the model and interpolate to seal
%%% lat/lon/times 
next_n = 0; 
for n = 1:length(dumpIters)-1
  
  %%% Keep track of where we are in the computation
  disp(n)
  
  %%% Days since start of model calendar year
  this_model_yearday = datenum('2008-01-01') + dumpIters(n)*deltaT/86400 - datenum('2011-01-01');
  next_model_yearday = datenum('2008-01-01') + dumpIters(n+1)*deltaT/86400 - datenum('2011-01-01'); 
  
  %%% Determine whether any seal measurements were taken between this model
  %%% yearday and the next one.
  times_found = false;
  seal_idx = cell(1,Nseals);
  for m = 1:Nseals
    seal_yearday = seal_yeardays{m};
    idx = find((seal_yearday>=this_model_yearday) & (seal_yearday<next_model_yearday));    
    if (~isempty(idx))
      times_found = true;      
    end
    seal_idx{m} = idx;
  end
  %%% If not then we don't need to load the model output
  if (~times_found)
    continue;
  end
    
  %%% Load instantaneous model state at the current model time 
  this_n = n;
  if (this_n == next_n) %%% If 'next' variables store data for 'this' time step then just copy it over
    this_u = next_u;
    this_v = next_v;
    this_t = next_t;
    this_s = next_s;
  else %%% Otherwise we need to read in the data for 'this' time step
    this_u = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),dumpIters(n));      
    this_v = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),dumpIters(n));
    this_t = rdmdsWrapper(fullfile(exppath,'/results/THETA_12hourly'),dumpIters(n));      
    this_s = rdmdsWrapper(fullfile(exppath,'/results/SALT_12hourly'),dumpIters(n));
    if (isempty(this_u) || isempty(this_v) || isempty(this_s) || isempty(this_t))   
      error(['No data found at iteration number ',num2str(dumpIters(n))]);
    end
  end
  
  %%% Load instantaneous model state at the next model time
  next_n = n + 1;
  next_u = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),dumpIters(n+1));      
  next_v = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),dumpIters(n+1));
  next_t = rdmdsWrapper(fullfile(exppath,'/results/THETA_12hourly'),dumpIters(n+1));      
  next_s = rdmdsWrapper(fullfile(exppath,'/results/SALT_12hourly'),dumpIters(n+1));
  if (isempty(next_u) || isempty(next_v) || isempty(next_s) || isempty(next_t))   
    error(['No data found at iteration number ',num2str(dumpIters(n+1))]);
  end
  
  %%% Interpolate trilinearly to obtain model state at seal location
  for m = 1:Nseals
    
    %%% Look over indices (in seal time) at which the measurement time lies
    %%% between 'this' model time and the 'next' model time
    for p = seal_idx{m}
      
      %%% Extract seal spatio-temporal location
      seal_lon = seal_lons{m}(p);
      seal_lat = seal_lats{m}(p);
      seal_yearday = seal_yeardays{m}(p);
      
      %%% Nearest model gridpoints
      i = find(model_lons<seal_lon,1,'last');
      j = find(model_lats<seal_lat,1,'last');
      
      %%% Fractional location of seal relative to model spatial and
      %%% temporal points
      r_lon = (seal_lon-model_lons(i))/(model_lons(i+1)-model_lons(i)); 
      r_lat = (seal_lat-model_lats(j))/(model_lats(j+1)-model_lats(j)); 
      r_time = (seal_yearday-this_model_yearday)/(next_model_yearday-this_model_yearday);
      
      %%% Interpolation weights
      wmmm = (1-r_lon).*(1-r_lat).*(1-r_time);
      wpmm = (r_lon).*(1-r_lat).*(1-r_time);
      wmpm = (1-r_lon).*(r_lat).*(1-r_time);
      wppm = (r_lon).*(r_lat).*(1-r_time);
      wmmp = (1-r_lon).*(1-r_lat).*(r_time);
      wpmp = (r_lon).*(1-r_lat).*(r_time);
      wmpp = (1-r_lon).*(r_lat).*(r_time);
      wppp = (r_lon).*(r_lat).*(r_time);
      
      %%% Perform trilinear interpolation
      seal_uvel{m}(p,:) = ...
          wmmm.*squeeze(this_u(i,j,:)) ...
        + wpmm.*squeeze(this_u(i+1,j,:)) ...
        + wmpm.*squeeze(this_u(i,j+1,:)) ...
        + wppm.*squeeze(this_u(i+1,j+1,:)) ...
        + wmmp.*squeeze(next_u(i,j,:)) ...
        + wpmp.*squeeze(next_u(i+1,j,:)) ...
        + wmpp.*squeeze(next_u(i,j+1,:)) ...
        + wppp.*squeeze(next_u(i+1,j+1,:));
      seal_vvel{m}(p,:) = ...
          wmmm.*squeeze(this_v(i,j,:)) ...
        + wpmm.*squeeze(this_v(i+1,j,:)) ...
        + wmpm.*squeeze(this_v(i,j+1,:)) ...
        + wppm.*squeeze(this_v(i+1,j+1,:)) ...
        + wmmp.*squeeze(next_v(i,j,:)) ...
        + wpmp.*squeeze(next_v(i+1,j,:)) ...
        + wmpp.*squeeze(next_v(i,j+1,:)) ...
        + wppp.*squeeze(next_v(i+1,j+1,:));
      seal_salt{m}(p,:) = ...
          wmmm.*squeeze(this_s(i,j,:)) ...
        + wpmm.*squeeze(this_s(i+1,j,:)) ...
        + wmpm.*squeeze(this_s(i,j+1,:)) ...
        + wppm.*squeeze(this_s(i+1,j+1,:)) ...
        + wmmp.*squeeze(next_s(i,j,:)) ...
        + wpmp.*squeeze(next_s(i+1,j,:)) ...
        + wmpp.*squeeze(next_s(i,j+1,:)) ...
        + wppp.*squeeze(next_s(i+1,j+1,:));
      seal_theta{m}(p,:) = ...
          wmmm.*squeeze(this_t(i,j,:)) ...
        + wpmm.*squeeze(this_t(i+1,j,:)) ...
        + wmpm.*squeeze(this_t(i,j+1,:)) ...
        + wppm.*squeeze(this_t(i+1,j+1,:)) ...
        + wmmp.*squeeze(next_t(i,j,:)) ...
        + wpmp.*squeeze(next_t(i+1,j,:)) ...
        + wmpp.*squeeze(next_t(i,j+1,:)) ...
        + wppp.*squeeze(next_t(i+1,j+1,:));
      
    end
  end
 
end

%%% Write products to output files
for m = 1:Nseals
  outfname = ['model',num2str(seal_numbers(m)),'.mat'];
  lon = seal_lons{m};
  lat = seal_lats{m};
  time = seal_times{m};
  yearday = seal_yeardays{m};
  uvel = seal_uvel{m};
  vvel = seal_vvel{m};
  salt = seal_salt{m};
  theta = seal_theta{m};
  save(fullfile(sealdir,outfname),'lon','lat','time','yearday','uvel','vvel','salt','theta','depth');
end