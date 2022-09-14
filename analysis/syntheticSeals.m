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
seal_lons = cell(1,2*Nseals);
seal_lats = cell(1,2*Nseals);
seal_times = cell(1,2*Nseals);
seal_yeardays = cell(1,2*Nseals);
Ce = 40000; %%% Approximate Earth circumference in km
for m = 1:2:2*Nseals-1
  
  %%% First we just store the seal lats, lons, datenums and yeardays
  load(fullfile(sealdir,['seal',num2str(seal_numbers((m+1)/2)),'.mat']));
  seal_times{m} = tag_datenum';
  yeardays = 0*tag_datenum;
  for j=1:length(yeardays)
    yeardays(j) = tag_datenum(j) - datenum([datestr(tag_datenum(j),'yyyy'),'-01-01']);
  end
  seal_yeardays{m} = yeardays; 
  seal_lons{m} = tag_lat_lon(:,2);
  seal_lats{m} = tag_lat_lon(:,1);
  
  %%% Now we create a synthetic seal at higher resolution
  seal_times{m+1} = [tag_datenum(1)];  
  seal_lons{m+1} = [tag_lat_lon(1,2)];
  seal_lats{m+1} = [tag_lat_lon(1,1)];
  for j = 1:length(yeardays)-1
    seal_dist = distance(tag_lat_lon(j,1),tag_lat_lon(j,2),tag_lat_lon(j+1,1),tag_lat_lon(j+1,2))*Ce/360; %%% Distance between seal measurements in km
    Nsplit = ceil(seal_dist/1); %%% Split into subdivisions no longer than 1 km
    dlat = (tag_lat_lon(j+1,1)-tag_lat_lon(j,1)) / Nsplit; %%% Subdivision lengths in lat/lon/time
    dlon = (tag_lat_lon(j+1,2)-tag_lat_lon(j,2)) / Nsplit; 
    dtime = (tag_datenum(j+1)-tag_datenum(j))/Nsplit;
    for i = 1:Nsplit-1
      seal_times{m+1} = [seal_times{m+1} tag_datenum(j)+i*dtime];
      seal_lons{m+1} = [seal_lons{m+1} tag_lat_lon(j,2)+i*dlon];
      seal_lats{m+1} = [seal_lats{m+1} tag_lat_lon(j,1)+i*dlat];  
    end
    seal_times{m+1} = [seal_times{m+1} tag_datenum(j+1)];
    seal_lons{m+1} = [seal_lons{m+1} tag_lat_lon(j+1,2)];
    seal_lats{m+1} = [seal_lats{m+1} tag_lat_lon(j+1,1)];
  end
  yeardays = 0*seal_times{m+1};
  for j=1:length(yeardays)
    yeardays(j) = seal_times{m+1}(j) - datenum([datestr(seal_times{m+1}(j),'yyyy'),'-01-01']);
  end
  seal_yeardays{m+1} = yeardays; 
  
end

%%% Read experiment data
setExpname; %%% To set directory paths
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2'; %%% Enforce high-resolution experiment
loadexp;
model_depth = -squeeze(zz);
model_lons_C = XC(:,1); %%% Regular lat/lon grid so we can do this
model_lats_C = YC(1,:);
model_lons_G = XG(:,1); 
model_lats_G = YG(1,:);
DXC3D = repmat(DXC,[1 1 Nr]);
DYC3D = repmat(DYC,[1 1 Nr]);
RAZ3D = repmat(RAZ,[1 1 Nr]);

%%% Model iteration numbers at which time-averaged output is saved
dumpIters_avg = [1534200 1578034 1621869 1665703 1709537 1753371 1797206 ...
               1841040 1884874 1928709 1972543 2016377 2060211 2104046 2147880];
nDumps_avg = length(dumpIters_avg);
yeardays_avg = round(datenum('2008-01-01') + ((dumpIters_avg(2:end)+dumpIters_avg(1:end-1))/2)*deltaT/86400 ...
                  - datenum('2011-01-01')); %%% Actual times to which time-averaged output correspond, i.e. mid-month
dumpIters_avg(1) = []; %%% No longer need the first output file, which corresponds to November 2010
this_avg_day = 0;
next_avg_day = yeardays_avg(1); %%% Must be smaller than the first instantaneous model yearday
this_avg_n = 0;
next_avg_n = 1;

%%% Frequency of diagnostic output
dumpStart_inst = 1578240;
dumpStep_inst = 43200/60;
nDumps_inst = 731;
dumpIters_inst = dumpStart_inst:dumpStep_inst:dumpStart_inst+(nDumps_inst-1)*dumpStep_inst;
 
%%% To store model state at seal measurement locations
seal_uvel = cell(1,2*Nseals);
seal_vvel = cell(1,2*Nseals);
seal_salt = cell(1,2*Nseals);
seal_theta = cell(1,2*Nseals);
seal_zeta = cell(1,2*Nseals);
seal_wteddy = cell(1,2*Nseals);
seal_wseddy = cell(1,2*Nseals);
for m = 1:2*Nseals
  seal_uvel{m} = zeros(length(seal_times{m}),Nr);
  seal_vvel{m} = zeros(length(seal_times{m}),Nr);
  seal_salt{m} = zeros(length(seal_times{m}),Nr);
  seal_theta{m} = zeros(length(seal_times{m}),Nr);
  seal_zeta{m} = zeros(length(seal_times{m}),Nr);
  seal_wteddy{m} = zeros(length(seal_times{m}),Nr);
  seal_wseddy{m} = zeros(length(seal_times{m}),Nr);
end

%%% Loop through time outputs from the model and interpolate to seal
%%% lat/lon/times 
next_n = 0; 

for n = 1:length(dumpIters_inst)-1
  
  %%% Keep track of where we are in the computation
  disp(n)
  
  %%% Days since start of model calendar year
  this_model_yearday = datenum('2008-01-01') + dumpIters_inst(n)*deltaT/86400 - datenum('2011-01-01');
  next_model_yearday = datenum('2008-01-01') + dumpIters_inst(n+1)*deltaT/86400 - datenum('2011-01-01'); 
  
  %%% If we have passed a model average output time step, load the next month's vertical eddy fluxes
  if (this_model_yearday >= next_avg_yearday)
    
    if (n > 1) %%% A bit sloppy, but this condition ensures that we load the first month's vertical eddy fluxes. 
               %%%It assumes that the first couple of instantaneous outputs precede any of the seal measurements
               
      %%% 'next' month becomes 'this' month
      this_avg_yearday = next_avg_yearday;
      this_avg_wteddy = next_avg_wteddy;
      this_avg_wseddy = next_avg_wseddy;      
      this_avg_n = next_avg_n;
      
    end
    
    next_avg_w = rdmdsWrapper(fullfile(exppath,'/results/WVEL'),dumpIters_avg(next_avg_n));   
    next_avg_t = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters_avg(next_avg_n));   
    next_avg_wt = rdmdsWrapper(fullfile(exppath,'/results/WVELTH'),dumpIters_avg(next_avg_n));          
    next_avg_wteddy = next_avg_wt - next_avg_w.*(0.5*(next_avg_t(:,:,1:Nr)+next_avg_t(:,:,[Nr 1:Nr-1]))); 
    clear('next_avg_t','next_avg_wt');
    next_avg_s = rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters_avg(next_avg_n));   
    next_avg_ws = rdmdsWrapper(fullfile(exppath,'/results/WVELSLT'),dumpIters_avg(next_avg_n));          
    next_avg_wseddy = next_avg_ws - next_avg_w.*(0.5*(next_avg_s(:,:,1:Nr)+next_avg_s(:,:,[Nr 1:Nr-1]))); 
    clear('next_avg_s','next_avg_ws','next_avg_w');    
    
    %%% Increment index for average output
    if (n > 1)
      next_avg_n = next_avg_n + 1;
      next_avg_yearday = yeardays_avg(next_avg_n);
    end
    
  end
  
  %%% Determine whether any seal measurements were taken between this model
  %%% yearday and the next one.
  times_found = false;
  seal_idx = cell(1,2*Nseals);
  for m = 1:2*Nseals
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
    this_z = next_z;
  else %%% Otherwise we need to read in the data for 'this' time step
    this_u = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),dumpIters_inst(n));      
    this_v = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),dumpIters_inst(n));
    this_t = rdmdsWrapper(fullfile(exppath,'/results/THETA_12hourly'),dumpIters_inst(n));      
    this_s = rdmdsWrapper(fullfile(exppath,'/results/SALT_12hourly'),dumpIters_inst(n));
    if (isempty(this_u) || isempty(this_v) || isempty(this_s) || isempty(this_t))   
      error(['No data found at iteration number ',num2str(dumpIters_inst(n))]);
    end
    this_u(hFacW==0) = NaN;
    this_v(hFacS==0) = NaN;
    this_t(hFacC==0) = NaN;
    this_s(hFacC==0) = NaN;    
    this_z = ( this_u(:,[Ny 1:Ny-1],:).*DXC3D(:,[Ny 1:Ny-1],:) ...
             + this_v([1:Nx],:,:).*DYC3D([1:Nx],:,:) ...
             - this_u(:,[1:Ny],:).*DXC3D(:,[1:Ny],:) ...
             - this_v([Nx 1:Nx-1],:,:).*DYC3D([Nx 1:Nx-1],:,:) ) ./ RAZ3D;
  end
  
  %%% Load instantaneous model state at the next model time
  next_n = n + 1;
  next_u = rdmdsWrapper(fullfile(exppath,'/results/UVEL_12hourly'),dumpIters_inst(n+1));      
  next_v = rdmdsWrapper(fullfile(exppath,'/results/VVEL_12hourly'),dumpIters_inst(n+1));
  next_t = rdmdsWrapper(fullfile(exppath,'/results/THETA_12hourly'),dumpIters_inst(n+1));      
  next_s = rdmdsWrapper(fullfile(exppath,'/results/SALT_12hourly'),dumpIters_inst(n+1));
  if (isempty(next_u) || isempty(next_v) || isempty(next_s) || isempty(next_t))   
    error(['No data found at iteration number ',num2str(dumpIters_inst(n+1))]);
  end
  next_u(hFacW==0) = NaN;
  next_v(hFacS==0) = NaN;
  next_t(hFacC==0) = NaN;
  next_s(hFacC==0) = NaN;    
  next_z = ( next_u(:,[Ny 1:Ny-1],:).*DXC3D(:,[Ny 1:Ny-1],:) ...
           + next_v([1:Nx],:,:).*DYC3D([1:Nx],:,:) ...
           - next_u(:,[1:Ny],:).*DXC3D(:,[1:Ny],:) ...
           - next_v([Nx 1:Nx-1],:,:).*DYC3D([Nx 1:Nx-1],:,:) ) ./ RAZ3D;
  
  %%% Interpolate trilinearly to obtain model state at seal location
  for m = 1:2*Nseals
    
    %%% Look over indices (in seal time) at which the measurement time lies
    %%% between 'this' model time and the 'next' model time
    for p = seal_idx{m}
      
      %%% Extract seal spatio-temporal location
      seal_lon = seal_lons{m}(p);
      seal_lat = seal_lats{m}(p);
      seal_yearday = seal_yeardays{m}(p);
      
      %%% Nearest model gridpoints
      iC = find(model_lons_C<seal_lon,1,'last');
      jC = find(model_lats_C<seal_lat,1,'last');
      iG = find(model_lons_G<seal_lon,1,'last');
      jG = find(model_lats_G<seal_lat,1,'last');
      
      %%% Fractional distance of seal across grid cells (in space and time)
      r_lonC = (seal_lon-model_lons_C(iC))/(model_lons_C(iC+1)-model_lons_C(iC)); 
      r_lonG = (seal_lon-model_lons_G(iG))/(model_lons_G(iG+1)-model_lons_G(iG)); 
      r_latC = (seal_lat-model_lats_C(jC))/(model_lats_C(jC+1)-model_lats_C(jC));
      r_latG = (seal_lat-model_lats_G(jG))/(model_lats_G(jG+1)-model_lats_G(jG)); 
      r_time_inst = (seal_yearday-this_model_yearday)/(next_model_yearday-this_model_yearday);   
      r_time_avg = (seal_yearday-this_avg_yearday)/(next_avg_yearday-this_avg_yearday);   
      
      %%% Interpolation weights (cell centers) for instantaneous output
      WmmmC_inst = (1-r_lonC).*(1-r_latC).*(1-r_time_inst);
      WpmmC_inst = (r_lonC).*(1-r_latC).*(1-r_time_inst);
      WmpmC_inst = (1-r_lonC).*(r_latC).*(1-r_time_inst);
      WppmC_inst = (r_lonC).*(r_latC).*(1-r_time_inst);
      WmmpC_inst = (1-r_lonC).*(1-r_latC).*(r_time_inst);
      WpmpC_inst = (r_lonC).*(1-r_latC).*(r_time_inst);
      WmppC_inst = (1-r_lonC).*(r_latC).*(r_time_inst);
      WpppC_inst = (r_lonC).*(r_latC).*(r_time_inst);
      
      %%% Interpolation weights (cell centers)      
      WmmmC_avg = (1-r_lonC).*(1-r_latC).*(1-r_time_avg);
      WpmmC_avg = (r_lonC).*(1-r_latC).*(1-r_time_avg);
      WmpmC_avg = (1-r_lonC).*(r_latC).*(1-r_time_avg);
      WppmC_avg = (r_lonC).*(r_latC).*(1-r_time_avg);
      WmmpC_avg = (1-r_lonC).*(1-r_latC).*(r_time_avg);
      WpmpC_avg = (r_lonC).*(1-r_latC).*(r_time_avg);
      WmppC_avg = (1-r_lonC).*(r_latC).*(r_time_avg);
      WpppC_avg = (r_lonC).*(r_latC).*(r_time_avg);
      
      %%% Interpolation weights (cell west faces)
      WmmmW = (1-r_lonG).*(1-r_latC).*(1-r_time_inst);
      WpmmW = (r_lonG).*(1-r_latC).*(1-r_time_inst);
      WmpmW = (1-r_lonG).*(r_latC).*(1-r_time_inst);
      WppmW = (r_lonG).*(r_latC).*(1-r_time_inst);
      WmmpW = (1-r_lonG).*(1-r_latC).*(r_time_inst);
      WpmpW = (r_lonG).*(1-r_latC).*(r_time_inst);
      WmppW = (1-r_lonG).*(r_latC).*(r_time_inst);
      WpppW = (r_lonG).*(r_latC).*(r_time_inst);
      
      %%% Interpolation weights (cell north faces)      
      WmmmS = (1-r_lonC).*(1-r_latG).*(1-r_time_inst);
      WpmmS = (r_lonC).*(1-r_latG).*(1-r_time_inst);
      WmpmS = (1-r_lonC).*(r_latG).*(1-r_time_inst);
      WppmS = (r_lonC).*(r_latG).*(1-r_time_inst);
      WmmpS = (1-r_lonC).*(1-r_latG).*(r_time_inst);
      WpmpS = (r_lonC).*(1-r_latG).*(r_time_inst);
      WmppS = (1-r_lonC).*(r_latG).*(r_time_inst);
      WpppS = (r_lonC).*(r_latG).*(r_time_inst);
      
      %%% Interpolation weights (cell corners)             
      WmmmZ = (1-r_lonG).*(1-r_latG).*(1-r_time_inst);
      WpmmZ = (r_lonG).*(1-r_latG).*(1-r_time_inst);
      WmpmZ = (1-r_lonG).*(r_latG).*(1-r_time_inst);
      WppmZ = (r_lonG).*(r_latG).*(1-r_time_inst);
      WmmpZ = (1-r_lonG).*(1-r_latG).*(r_time_inst);
      WpmpZ = (r_lonG).*(1-r_latG).*(r_time_inst);
      WmppZ = (1-r_lonG).*(r_latG).*(r_time_inst);
      WpppZ = (r_lonG).*(r_latG).*(r_time_inst);
      
      %%% Perform trilinear interpolation
      seal_uvel{m}(p,:) = ...
          WmmmC_inst.*squeeze(this_u(iG,jC,:)) ...
        + WpmmC_inst.*squeeze(this_u(iG+1,jC,:)) ...
        + WmpmC_inst.*squeeze(this_u(iG,jC+1,:)) ...
        + WppmC_inst.*squeeze(this_u(iG+1,jC+1,:)) ...
        + WmmpC_inst.*squeeze(next_u(iG,jC,:)) ...
        + WpmpC_inst.*squeeze(next_u(iG+1,jC,:)) ...
        + WmppC_inst.*squeeze(next_u(iG,jC+1,:)) ...
        + WpppC_inst.*squeeze(next_u(iG+1,jC+1,:));
      seal_vvel{m}(p,:) = ...
          WmmmC_inst.*squeeze(this_v(iC,jG,:)) ...
        + WpmmC_inst.*squeeze(this_v(iC+1,jG,:)) ...
        + WmpmC_inst.*squeeze(this_v(iC,jG+1,:)) ...
        + WppmC_inst.*squeeze(this_v(iC+1,jG+1,:)) ...
        + WmmpC_inst.*squeeze(next_v(iC,jG,:)) ...
        + WpmpC_inst.*squeeze(next_v(iC+1,jG,:)) ...
        + WmppC_inst.*squeeze(next_v(iC,jG+1,:)) ...
        + WpppC_inst.*squeeze(next_v(iC+1,jG+1,:));
      seal_salt{m}(p,:) = ...
          WmmmC_inst.*squeeze(this_s(iC,jC,:)) ...
        + WpmmC_inst.*squeeze(this_s(iC+1,jC,:)) ...
        + WmpmC_inst.*squeeze(this_s(iC,jC+1,:)) ...
        + WppmC_inst.*squeeze(this_s(iC+1,jC+1,:)) ...
        + WmmpC_inst.*squeeze(next_s(iC,jC,:)) ...
        + WpmpC_inst.*squeeze(next_s(iC+1,jC,:)) ...
        + WmppC_inst.*squeeze(next_s(iC,jC+1,:)) ...
        + WpppC_inst.*squeeze(next_s(iC+1,jC+1,:));
      seal_theta{m}(p,:) = ...
          WmmmC_inst.*squeeze(this_t(iC,jC,:)) ...
        + WpmmC_inst.*squeeze(this_t(iC+1,jC,:)) ...
        + WmpmC_inst.*squeeze(this_t(iC,jC+1,:)) ...
        + WppmC_inst.*squeeze(this_t(iC+1,jC+1,:)) ...
        + WmmpC_inst.*squeeze(next_t(iC,jC,:)) ...
        + WpmpC_inst.*squeeze(next_t(iC+1,jC,:)) ...
        + WmppC_inst.*squeeze(next_t(iC,jC+1,:)) ...
        + WpppC_inst.*squeeze(next_t(iC+1,jC+1,:));
      seal_zeta{m}(p,:) = ...
          WmmmZ.*squeeze(this_z(iG,jG,:)) ...
        + WpmmZ.*squeeze(this_z(iG+1,jG,:)) ...
        + WmpmZ.*squeeze(this_z(iG,jG+1,:)) ...
        + WppmZ.*squeeze(this_z(iG+1,jG+1,:)) ...
        + WmmpZ.*squeeze(next_z(iG,jG,:)) ...
        + WpmpZ.*squeeze(next_z(iG+1,jG,:)) ...
        + WmppZ.*squeeze(next_z(iG,jG+1,:)) ...
        + WpppZ.*squeeze(next_z(iG+1,jG+1,:));
      seal_wseddy{m}(p,:) = ...
          WmmmC_avg.*squeeze(this_avg_wseddy(iC,jC,:)) ...
        + WpmmC_avg.*squeeze(this_avg_wseddy(iC+1,jC,:)) ...
        + WmpmC_avg.*squeeze(this_avg_wseddy(iC,jC+1,:)) ...
        + WppmC_avg.*squeeze(this_avg_wseddy(iC+1,jC+1,:)) ...
        + WmmpC_avg.*squeeze(next_avg_wseddy(iC,jC,:)) ...
        + WpmpC_avg.*squeeze(next_avg_wseddy(iC+1,jC,:)) ...
        + WmppC_avg.*squeeze(next_avg_wseddy(iC,jC+1,:)) ...
        + WpppC_avg.*squeeze(next_avg_wseddy(iC+1,jC+1,:));
      seal_wteddy{m}(p,:) = ...
          WmmmC_avg.*squeeze(this_avg_wteddy(iC,jC,:)) ...
        + WpmmC_avg.*squeeze(this_avg_wteddy(iC+1,jC,:)) ...
        + WmpmC_avg.*squeeze(this_avg_wteddy(iC,jC+1,:)) ...
        + WppmC_avg.*squeeze(this_avg_wteddy(iC+1,jC+1,:)) ...
        + WmmpC_avg.*squeeze(next_avg_wteddy(iC,jC,:)) ...
        + WpmpC_avg.*squeeze(next_avg_wteddy(iC+1,jC,:)) ...
        + WmppC_avg.*squeeze(next_avg_wteddy(iC,jC+1,:)) ...
        + WpppC_avg.*squeeze(next_avg_wteddy(iC+1,jC+1,:));
      
    end
  end  
 
end

%%% Write products to output files
for m = 1:2*Nseals
  if (mod(m,2)==1)
    outfname = ['model',num2str(seal_numbers((m+1)/2)),'.mat'];
  else
    outfname = ['interp',num2str(seal_numbers(m/2)),'.mat'];
  end
  lon = seal_lons{m};
  lat = seal_lats{m};
  time = seal_times{m};
  yearday = seal_yeardays{m};
  uvel = seal_uvel{m};
  vvel = seal_vvel{m};
  salt = seal_salt{m};
  theta = seal_theta{m};
  zeta = seal_zeta{m};
  wteddy = seal_wteddy{m};
  wseddy = seal_wseddy{m};
  save(fullfile(sealdir,outfname),'lon','lat','time','yearday','uvel','vvel','salt','theta','zeta','wteddy','wseddy','depth');
end