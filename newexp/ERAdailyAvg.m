%%%
%%% ERAdailyAvg.m
%%%
%%% Helper function that averages ERA-Interim data over a single day and interpolates it onto the model grid.
%%%
function avgdata = ERAdailyAvg (years,months,days,datadir,filename,varname,username,XMC,YMC,time_weight) 
  
  %%% Data file path
  if (years<=2010)

    %%% 2007-2010 format (yearly files)
    datafname = fullfile(datadir,[filename,'.',num2str(years),'.',username,'.nc']);
    
  else
    
    %%% 2010-onward format (monthly files)
    datafname = fullfile(datadir,[filename,'.',num2str(years),num2str(months,'%.2d'),'.',username,'.nc']);
    
  end

  %%% Grid for current time period
  lats=ncread(datafname,'g4_lat_1');
  lons=ncread(datafname,'g4_lon_2');
  time=ncread(datafname,'initial_time0_hours');

  %%% Switch longitude convention
  idx_west = find(lons>=180);
  idx_east = find(lons<180);
  lons = [lons(idx_west)-360 lons(idx_east)];

  %%% Meshgrid for interpolation
  [LO,LA] = meshgrid(lons,lats);

  %%% Convert current day to time in hours since 1800-01-01
  thedate = (datenum([num2str(years),'-',num2str(months),'-',num2str(days)]) - datenum('1800-01-01')) * 24;

  %%% Find index of this day in the dataset
  dayidx = find(time==thedate,1,'first');

  %%% Calculate daily average 
  %%% Need to read two outputs per day because output is 12-hourly
  avgdata = ncread(datafname,varname,[1 1 dayidx],[Inf Inf 2]);
  avgdata = mean(avgdata,3);
  
  %%% If the data is already weighted by time then divide by the weight
  %%% (0.5 days)
  if (time_weight)
    avgdata = avgdata / (86400/2);
  end

  %%% Switch longitude convention
  avgdata = [avgdata(idx_west,:); avgdata(idx_east,:)];

  %%% Interpolate onto our grid
  avgdata = interp2(LO,LA,avgdata',XMC,YMC,'linear')';
  
end
  