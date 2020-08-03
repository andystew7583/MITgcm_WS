%%%
%%% genOBCSagain.m
%%%
%%% Extracts slices of SOSE data along the boundaries of our model domain.
%%%

%%% Matlab utilities 

addpath ../newexp_utils
addpath ../utils/matlab


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Data format parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ieee='b';
prec='real*8';
realdigits = 8;
realfmt=['%.',num2str(realdigits),'e'];


%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%% SET-UP %%%%%%%%%%%%%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%  



%%% Defines our model domain dimensions
defineGrid

% name of SOSE data file
OBCS_data_dir = '../data/SOSEdata/08-10';
OBCS_storage_dir = './OBCS';
OBCS_data_dir_snow = '../data/SOSEdata/13-17';







%%% Default grid for SOSE variables
load sosegrid.mat

%%% Switch longitude convention and generate required grid indices
[idx_west,idx_east,idx_OBN,idx_OBE,XC,YC] = switchLons (XC,YC,xmin,xmax,ymin,ymax);

%%% SOSE sub-grid dimensions
SX = size(XC,1);
SY = size(XC,2);
SR = length(DRF);




%%% 5-day output, 2008-2010
Nrec = 219; %%% 219 records
data_freq = 5; %%% Output every 5 days
data_off = 3; %%% Starts on 3rd day of first month
days_per_month = calcDaysPerMonth(2008,2010);
month_off = 1; %%% First month is Jan 2008

disp('THETA');
[theta_obcs_north,theta_obcs_east,theta_obcs_west] = loadSOSErecs(OBCS_data_dir,'THETA',SX,SY,SR,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
theta_obcs_north = binByMonth(theta_obcs_north,data_freq,data_off,days_per_month,month_off);
theta_obcs_east = binByMonth(theta_obcs_east,data_freq,data_off,days_per_month,month_off);
theta_obcs_west = binByMonth(theta_obcs_west,data_freq,data_off,days_per_month,month_off);

disp('SALT');
[salt_obcs_north,salt_obcs_east,salt_obcs_west] = loadSOSErecs(OBCS_data_dir,'SALT',SX,SY,SR,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
salt_obcs_north = binByMonth(salt_obcs_north,data_freq,data_off,days_per_month,month_off);
salt_obcs_east = binByMonth(salt_obcs_east,data_freq,data_off,days_per_month,month_off);
salt_obcs_west = binByMonth(salt_obcs_west,data_freq,data_off,days_per_month,month_off);

disp('UVEL');
[uvel_obcs_north,uvel_obcs_east,uvel_obcs_west] = loadSOSErecs(OBCS_data_dir,'UVEL',SX,SY,SR,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
uvel_obcs_north = binByMonth(uvel_obcs_north,data_freq,data_off,days_per_month,month_off);
uvel_obcs_east = binByMonth(uvel_obcs_east,data_freq,data_off,days_per_month,month_off);
uvel_obcs_west = binByMonth(uvel_obcs_west,data_freq,data_off,days_per_month,month_off);

disp('VVEL');
[vvel_obcs_north,vvel_obcs_east,vvel_obcs_west] = loadSOSErecs(OBCS_data_dir,'VVEL',SX,SY,SR,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
vvel_obcs_north = binByMonth(vvel_obcs_north,data_freq,data_off,days_per_month,month_off);
vvel_obcs_east = binByMonth(vvel_obcs_east,data_freq,data_off,days_per_month,month_off);
vvel_obcs_west = binByMonth(vvel_obcs_west,data_freq,data_off,days_per_month,month_off);

disp('PHIHYD');
[PHIHYD_obcs_north,PHIHYD_obcs_east,PHIHYD_obcs_west] = loadSOSErecs(OBCS_data_dir,'PHIHYD',SX,SY,SR,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
PHIHYD_obcs_north = binByMonth(PHIHYD_obcs_north,data_freq,data_off,days_per_month,month_off);
PHIHYD_obcs_east = binByMonth(PHIHYD_obcs_east,data_freq,data_off,days_per_month,month_off);
PHIHYD_obcs_west = binByMonth(PHIHYD_obcs_west,data_freq,data_off,days_per_month,month_off);

disp('SIheff');
[SIThick_obcs_north,SIThick_obcs_east,SIThick_obcs_west] = loadSOSErecs(OBCS_data_dir,'SIheff',SX,SY,1,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
SIThick_obcs_north = binByMonth(SIThick_obcs_north,data_freq,data_off,days_per_month,month_off);
SIThick_obcs_east = binByMonth(SIThick_obcs_east,data_freq,data_off,days_per_month,month_off);
SIThick_obcs_west = binByMonth(SIThick_obcs_west,data_freq,data_off,days_per_month,month_off);

%%% daily output, 2008-2010
Nrec = 1096;
data_freq = 1; %%% Output every 1 day
data_off = 1; %%% Starts on 1st day of first month
days_per_month = calcDaysPerMonth(2008,2010);
month_off = 1; %%% First month is Jan 2008

disp('IceConc');
[SIArea_obcs_north,SIArea_obcs_east,SIArea_obcs_west] = loadSOSErecs(OBCS_data_dir,'IceConc',SX,SY,1,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
SIArea_obcs_north = binByMonth(SIArea_obcs_north,data_freq,data_off,days_per_month,month_off);
SIArea_obcs_east = binByMonth(SIArea_obcs_east,data_freq,data_off,days_per_month,month_off);
SIArea_obcs_west = binByMonth(SIArea_obcs_west,data_freq,data_off,days_per_month,month_off);

%%% 5-day output, 2005-2010
Nrec = 438;
data_freq = 5; %%% Output every 5 days
data_off = 3; %%% Starts on 3rd day of first month
days_per_month = calcDaysPerMonth(2005,2010);
month_off = 1; %%% First month is Jan 2005

disp('SIuice');
[SIUvel_obcs_north,SIUvel_obcs_east,SIUvel_obcs_west] = loadSOSErecs(OBCS_data_dir,'SIuice',SX,SY,1,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
SIUvel_obcs_north = binByMonth(SIUvel_obcs_north,data_freq,data_off,days_per_month,month_off);
SIUvel_obcs_east = binByMonth(SIUvel_obcs_east,data_freq,data_off,days_per_month,month_off);
SIUvel_obcs_west = binByMonth(SIUvel_obcs_west,data_freq,data_off,days_per_month,month_off);

disp('SIvice');
[SIVvel_obcs_north,SIVvel_obcs_east,SIVvel_obcs_west] = loadSOSErecs(OBCS_data_dir,'SIvice',SX,SY,1,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
SIVvel_obcs_north = binByMonth(SIVvel_obcs_north,data_freq,data_off,days_per_month,month_off);
SIVvel_obcs_east = binByMonth(SIVvel_obcs_east,data_freq,data_off,days_per_month,month_off);
SIVvel_obcs_west = binByMonth(SIVvel_obcs_west,data_freq,data_off,days_per_month,month_off);

%%% daily output, 2013-2017
Nrec = 1826;
data_freq = 1; %%% Output every 1 day
data_off = 1; %%% Starts on 1st day of first month
days_per_month = calcDaysPerMonth(2013,2017);
month_off = 1; %%% First month is Jan 2013
XC = ncread('../data/SOSEdata/13-17/SIhsnow.nc','XC'); %%% 2013-2017 solution is on a different grid
YC = ncread('../data/SOSEdata/13-17/SIhsnow.nc','YC');
XC = repmat(XC,[1 length(YC)]);
YC = repmat(YC',[size(XC,1) 1]);
[idx_west,idx_east,idx_OBN,idx_OBE,XC,YC] = switchLons (XC,YC,xmin,xmax,ymin,ymax);
SX = size(XC,1);
SY = size(XC,2);

disp('SIhsnow');
[SIHsnow_obcs_north,SIHsnow_obcs_east,SIHsnow_obcs_west] = loadSOSErecsNC(OBCS_data_dir_snow,'SIhsnow.nc','SIhsnow',SX,SY,1,Nrec,idx_west,idx_east,idx_OBN,idx_OBE);
SIHsnow_obcs_north = binByMonth(SIHsnow_obcs_north,data_freq,data_off,days_per_month,month_off);
SIHsnow_obcs_east = binByMonth(SIHsnow_obcs_east,data_freq,data_off,days_per_month,month_off);
SIHsnow_obcs_west = binByMonth(SIHsnow_obcs_west,data_freq,data_off,days_per_month,month_off);


%%% Write to .mat file
save(fullfile(OBCS_storage_dir,'OBCS.mat'), ...
  'theta_obcs_north','theta_obcs_east','theta_obcs_west', ...
  'salt_obcs_north','salt_obcs_east','salt_obcs_west', ...
  'uvel_obcs_north','uvel_obcs_east','uvel_obcs_west', ...
  'vvel_obcs_north','vvel_obcs_east','vvel_obcs_west', ...
  'PHIHYD_obcs_north','PHIHYD_obcs_east','PHIHYD_obcs_west', ...
  'SIThick_obcs_north','SIThick_obcs_east','SIThick_obcs_west', ...
  'SIArea_obcs_north','SIArea_obcs_east','SIArea_obcs_west', ...
  'SIUvel_obcs_north','SIUvel_obcs_east','SIUvel_obcs_west', ...
  'SIVvel_obcs_north','SIVvel_obcs_east','SIVvel_obcs_west', ...
  'SIHsnow_obcs_north','SIHsnow_obcs_east','SIHsnow_obcs_west');


 


%%% 
%%% loadSOSErecsNC
%%%
%%% Convenience function to load records from SOSE output in NetCDF format, switch
%%% longitudinal coordinates, and then retain only a subsection of the grid.
%%%
function [data_north,data_east,data_west] = loadSOSErecsNC(datadir,filename,varname,Nx,Ny,Nr,Nrec,lon_idx_west,lon_idx_east,region_lon_idx,region_lat_idx)

  %%% To store boundary data  
  data_north = zeros(Nx,1,Nr,Nrec);
  data_east = zeros(1,Ny,Nr,Nrec);
  data_west = zeros(1,Ny,Nr,Nrec);

  %%% Loop through records
  for i = 1:Nrec
    
    %%% Load next record
    data = ncread(fullfile(datadir,filename),varname,[1 1 i],[Inf Inf 1]);
    
    %%% Switch to correct longitude convention
    data = [data(lon_idx_west,:,:,:) ; data(lon_idx_east,:,:,:)];

    %%% Crop data to region of interest
    data = data(region_lon_idx,region_lat_idx,:,:); 
    
    %%% Extract boundary data
    data_north(:,:,:,i) = data(:,end,:,:);
    data_east(:,:,:,i) = data(end,:,:,:);
    data_west(:,:,:,i) = data(1,:,:,:);
    
  end

end


%%% 
%%% loadSOSErecs
%%%
%%% Convenience function to load records from SOSE output, switch
%%% longitudinal coordinates, and then retain only a subsection of the grid.
%%%
function [data_north,data_east,data_west] = loadSOSErecs(datadir,varname,Nx,Ny,Nr,Nrec,lon_idx_west,lon_idx_east,region_lon_idx,region_lat_idx)

  %%% To store boundary data  
  data_north = zeros(Nx,1,Nr,Nrec);
  data_east = zeros(1,Ny,Nr,Nrec);
  data_west = zeros(1,Ny,Nr,Nrec);

  %%% Loop through records
  for i = 1:Nrec
    
    %%% Load next record
    data = rdmds(fullfile(datadir,varname),'rec',i);
    
    %%% Switch to correct longitude convention
    data = [data(lon_idx_west,:,:,:) ; data(lon_idx_east,:,:,:)];

    %%% Crop data to region of interest
    data = data(region_lon_idx,region_lat_idx,:,:); 
    
    %%% Extract boundary data
    data_north(:,:,:,i) = data(:,end,:,:);
    data_east(:,:,:,i) = data(end,:,:,:);
    data_west(:,:,:,i) = data(1,:,:,:);
    
  end

end


%%%
%%% binByMonth
%%%
%%% Crudely averages data into monthly bins.
%%%
%%% data - Nx x Ny x Nr x Nrec array of data to be binned
%%% data_freq - Number of days between data outputs
%%% data_off - Number of days into the first month of the first record
%%% days_per_month - 1 x Nrec vector of numbers of days in each month over 
%%%                  the period spanned by the data records
%%% month_off - Index (1-12) of the month in which the data records start
%%%
function data_binned = binByMonth (data,data_freq,data_off,days_per_month,month_off)
 
  %%% Defined for clarity
  months_per_year = 12;  
  Nrec = size(data,4);
  
  %%% Total number of days up until the end of the current month
  days_cumul = cumsum(days_per_month);
  
  %%% To store binned data
  data_binned = zeros(size(data,1),size(data,2),size(data,3),months_per_year);
  
  %%% Counts the number of records added for each month of the year
  data_cntr = zeros(1,months_per_year);

  %%% Loop through the records and add each to the appropriate monthly bin
  day_cntr = data_off;
  month_cntr = 1;
  month_idx = mod(month_off-1,12)+1;
  for n = 1:Nrec
    
    %%% Add this record's data to the current monthly bin, and increment
    %%% the counter for that bin
    data_binned(:,:,:,month_idx) = data_binned(:,:,:,month_idx) + data(:,:,:,n);
    data_cntr(month_idx) = data_cntr(month_idx) + 1;
    
    %%% Increment day counter
    day_cntr = day_cntr + data_freq;
    
    %%% If we enter the next month, increment the month counter and month
    %%% index
    if (day_cntr > days_cumul(month_cntr))
      month_cntr = month_cntr + 1;
      month_idx = mod(month_idx,12) + 1;
    end
    
  end
  
  %%% Divide by number of records added to each bin to obtain the average
  for m = 1:months_per_year
    data_binned(:,:,:,m) = data_binned(:,:,:,m) / data_cntr(m);
  end  
 
end


%%% 
%%% switchLons
%%%
%%% Convenience function to switch grids from SOSE's [0,360] convention to
%%% our [-180,180] convention, generate indices for switching other
%%% matrices, and generate indices that restrict SOSE data to our model
%%% domain.
%%%
function [idx_west,idx_east,idx_OBN,idx_OBE,XC,YC] = switchLons (XC,YC,xmin,xmax,ymin,ymax)

  %%% Indices of "western" and "eastern" halves of SOSE domain
  idx_west = find(XC(:,1)>=180);
  idx_east = find(XC(:,1)<180);

  %%% Switch longitude and latitude grids around
  XC = [XC(idx_west,:)-360 ; XC(idx_east,:)];
  YC = [YC(idx_west,:) ; YC(idx_east,:)];

  %%% Indices defining the subset of SOSE that contains our model grid
  idx_OBN = find(XC(:,1)>xmin & XC(:,1)<xmax);
  idx_OBE = find(YC(1,:)>ymin & YC(1,:)<ymax);
  
  %%% Restrict grid matrices to our model domain
  XC = XC(idx_OBN,idx_OBE);
  YC = YC(idx_OBN,idx_OBE);

end


%%%
%%% calcDaysPerMonth
%%%
%%% Calculates days in each month over a specified range of months
%%%
function days_per_month = calcDaysPerMonth (startYear,endYear)

 days_per_month = zeros(1,(endYear-startYear+1)*12);
 
 year = startYear;
 for n = 1:length(days_per_month)
   switch (mod(n-1,12)+1)
     case 1
       days_per_month(n) = 31;
     case 2
       if (mod(year,4)==0)         
         days_per_month(n) = 29;
       else
         days_per_month(n) = 28;
       end
     case 3
       days_per_month(n) = 31;
     case 4
       days_per_month(n) = 30;
     case 5
       days_per_month(n) = 31;
     case 6
       days_per_month(n) = 30;
     case 7
       days_per_month(n) = 31;
     case 8
       days_per_month(n) = 31;
     case 9
       days_per_month(n) = 30;
     case 10
       days_per_month(n) = 31;
     case 11
       days_per_month(n) = 30;
     case 12
       days_per_month(n) = 31;
       year = year + 1;
   end
 end

end



