%%%
%%% RiddingZerosDailyForcing.m
%%%
%%% This is basially a bunch of hacks to fix missing or spurious data in
%%% AMPS
%%%

%%% Load grid
defineGrid

%%% Select which fields to fix
fix_pressure = false;
fix_temp = false;
fix_rh = false;
fix_sw = false;
fix_lw = false;
fix_precip = false;
fix_uwind = true;
fix_vwind = true;

%%%%%% write dataset path
addpath ../newexp_utils

% days = 26292;
days = 3287;
gendir = '/data3';

inputfr = '/data3/MITgcm_WS/newexp/DEFAULTS/input/';









%%% Temperature
if (fix_temp)

  forcingvar1 = zeros(EXF_Nx,EXF_Ny,days);
  fid = fopen(fullfile(inputfr,aTemp),'r','b');
  for k=1:days
    forcingvar1(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
  end
  fclose(fid);

  forcingvar1_mean = squeeze(mean(mean(forcingvar1,2),1));
  k = 1;
  while (k <= Ndays)

    %%% Check whether this day is missing data
    if (forcingvar1_mean(k) < -60)
      k
      kstart = k;

      %%% If so, find the last day following this one that is also missing
      %%% data
      kend = kstart;
      while (forcingvar1_mean(kend) < -60)
        kend = kend + 1;
      end
      kend = kend -1;

      %%% Interpolate linearly over the group of days that are missing data
      for n = kstart:kend
        forcingvar1(:,:,n) = (forcingvar1(:,:,kstart-1) * ((kend+1) - n) + forcingvar1(:,:,kend+1) * (n - (kstart-1)) ) / ((kend+1) - (kstart-1));
      end

      %%% Continue the search after the end of this group of days
      k = kend;

    end

    %%% Increment day index
    k = k + 1;

  end

  writeDataset(forcingvar1,fullfile(inputfr,aTemp),ieee,prec);

  clear forcingvar1

end




%%% SW
if (fix_sw)
  
  forcingvar1 = zeros(EXF_Nx,EXF_Ny,days);
  fid = fopen(fullfile(inputfr,aSW),'r','b');
  for k=1:days
    forcingvar1(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
  end
  fclose(fid);

  forcingvar1_mean = squeeze(mean(mean(forcingvar1,2),1));
  k = 1;
  while (k <= Ndays)

    %%% Check whether this day is missing data
    if ((forcingvar1_mean(k)==0) || ((k>1) && (abs(forcingvar1_mean(k)-forcingvar1_mean(k-1))>50)))
      k
      kstart = k;

      %%% If so, find the last day following this one that is also missing
      %%% data
      kend = kstart+1;
      while (forcingvar1_mean(kend)==0)
        kend = kend + 1;
      end
      kend = kend -1;

      %%% Interpolate linearly over the group of days that are missing data
      if (kend-kstart < 10) %%% Exclude long periods to avoid interpolating across austral winter
        for n = kstart:kend
          forcingvar1(:,:,n) = (forcingvar1(:,:,kstart-1) * ((kend+1) - n) + forcingvar1(:,:,kend+1) * (n - (kstart-1)) ) / ((kend+1) - (kstart-1));
        end
      end

      %%% Continue the search after the end of this group of days
      k = kend;

    end

    %%% Increment day index
    k = k + 1;

  end
  
  %%% SW values are way too high in later parts of time series, so replace
  %%% with first year
  

  writeDataset(forcingvar1,fullfile(inputfr,aSW),ieee,prec);

  clear forcingvar1

end



%%% Relative humidity
if (fix_rh)

  forcingvar1 = zeros(EXF_Nx,EXF_Ny,days);
  fid = fopen(fullfile(inputfr,anewAQ),'r','b');
  for k=1:days
    forcingvar1(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
  end
  fclose(fid);

  forcingvar1_mean = squeeze(mean(mean(forcingvar1,2),1));
  k = 1;
  while (k <= Ndays)

    %%% Check whether this day is missing data
    if (forcingvar1_mean(k) == 0)
      k
      kstart = k;

      %%% If so, find the last day following this one that is also missing
      %%% data
      kend = kstart;
      while (forcingvar1_mean(kend) == 0)
        kend = kend + 1;
      end
      kend = kend -1;

      %%% Interpolate linearly over the group of days that are missing data
      for n = kstart:kend
        forcingvar1(:,:,n) = (forcingvar1(:,:,kstart-1) * ((kend+1) - n) + forcingvar1(:,:,kend+1) * (n - (kstart-1)) ) / ((kend+1) - (kstart-1));
      end

      %%% Continue the search after the end of this group of days
      k = kend;

    end

    %%% Increment day index
    k = k + 1;

  end

  writeDataset(forcingvar1,fullfile(inputfr,anewAQ),ieee,prec);

  clear forcingvar1

end








%%% Zonal wind
if (fix_uwind)

  forcingvar1 = zeros(EXF_Nx,EXF_Ny,days);
  fid = fopen(fullfile(inputfr,zwind),'r','b');
  for k=1:days
    forcingvar1(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
  end
  fclose(fid);

  forcingvar1_mean = squeeze(mean(mean(forcingvar1,2),1));
  k = 1;
  while (k <= Ndays)

    %%% Check whether this day is missing data
    if (forcingvar1_mean(k) == 0)
      k
      kstart = k;

      %%% If so, find the last day following this one that is also missing
      %%% data
      kend = kstart;
      while (forcingvar1_mean(kend) == 0)
        kend = kend + 1;
      end
      kend = kend -1;

      %%% Interpolate linearly over the group of days that are missing data
      for n = kstart:kend
        forcingvar1(:,:,n) = (forcingvar1(:,:,kstart-1) * ((kend+1) - n) + forcingvar1(:,:,kend+1) * (n - (kstart-1)) ) / ((kend+1) - (kstart-1));
      end

      %%% Continue the search after the end of this group of days
      k = kend;

    end

    %%% Increment day index
    k = k + 1;

  end

  writeDataset(forcingvar1,fullfile(inputfr,zwind),ieee,prec);

  clear forcingvar1

end






%%% Meridional wind
if (fix_vwind)

  forcingvar1 = zeros(EXF_Nx,EXF_Ny,days);
  fid = fopen(fullfile(inputfr,mwind),'r','b');
  for k=1:days
    forcingvar1(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
  end
  fclose(fid);

  forcingvar1_mean = squeeze(mean(mean(forcingvar1,2),1));
  k = 1;
  while (k <= Ndays)

    %%% Check whether this day is missing data
    if (forcingvar1_mean(k) == 0)
      k
      kstart = k;

      %%% If so, find the last day following this one that is also missing
      %%% data
      kend = kstart;
      while (forcingvar1_mean(kend) == 0)
        kend = kend + 1;
      end
      kend = kend -1;

      %%% Interpolate linearly over the group of days that are missing data
      for n = kstart:kend
        forcingvar1(:,:,n) = (forcingvar1(:,:,kstart-1) * ((kend+1) - n) + forcingvar1(:,:,kend+1) * (n - (kstart-1)) ) / ((kend+1) - (kstart-1));
      end

      %%% Continue the search after the end of this group of days
      k = kend;

    end

    %%% Increment day index
    k = k + 1;

  end

  writeDataset(forcingvar1,fullfile(inputfr,mwind),ieee,prec);

  clear forcingvar1

end







%%% Longwave
if (fix_lw)

  forcingvar1 = zeros(EXF_Nx,EXF_Ny,days);
  fid = fopen(fullfile(inputfr,aLW),'r','b');
  for k=1:days
    forcingvar1(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
  end
  fclose(fid);

  forcingvar1_mean = squeeze(mean(mean(forcingvar1,2),1));
  k = 1;
  while (k <= Ndays)

    %%% Check whether this day is missing data
    if (forcingvar1_mean(k) == 0)
      k
      kstart = k;

      %%% If so, find the last day following this one that is also missing
      %%% data
      kend = kstart;
      while (forcingvar1_mean(kend) == 0)
        kend = kend + 1;
      end
      kend = kend -1;

      %%% Interpolate linearly over the group of days that are missing data
      for n = kstart:kend
        forcingvar1(:,:,n) = (forcingvar1(:,:,kstart-1) * ((kend+1) - n) + forcingvar1(:,:,kend+1) * (n - (kstart-1)) ) / ((kend+1) - (kstart-1));
      end

      %%% Continue the search after the end of this group of days
      k = kend;

    end

    %%% Increment day index
    k = k + 1;

  end

  writeDataset(forcingvar1,fullfile(inputfr,aLW),ieee,prec);

  clear forcingvar1

end









%%% Precip
if (fix_precip)

  forcingvar1 = zeros(EXF_Nx,EXF_Ny,days);
  fid = fopen(fullfile(inputfr,aPrecip),'r','b');
  for k=1:days
    forcingvar1(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
  end
  fclose(fid);

  %%% First couple of years (20km AMPS) has negative precip, so replace with
  %%% years 3-4
  forcingvar1(:,:,1:731) = forcingvar1(:,:,732:1462);

  writeDataset(forcingvar1,fullfile(inputfr,aPrecip),ieee,prec);

  clear forcingvar1

end











%%% Pressure
if (fix_pressure)
  
  forcingvar1 = zeros(EXF_Nx,EXF_Ny,days);
  fid = fopen(fullfile(inputfr,pressure),'r','b');
  for k=1:days
    tmpvar = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    forcingvar1(:,:,k) = tmpvar;
  end
  fclose(fid);

  forcingvar1_mean = squeeze(mean(mean(forcingvar1,2),1));
  k = 1;
  while (k <= Ndays)

    %%% Check whether this day is missing data
    if (forcingvar1_mean(k) < 8.6e4)
      k
      kstart = k;

      %%% If so, find the last day following this one that is also missing
      %%% data
      kend = kstart;
      while (forcingvar1_mean(kend) < 8.6e4)
        kend = kend + 1;
      end
      kend = kend -1;

      %%% Interpolate linearly over the group of days that are missing data
      for n = kstart:kend
        forcingvar1(:,:,n) = (forcingvar1(:,:,kstart-1) * ((kend+1) - n) + forcingvar1(:,:,kend+1) * (n - (kstart-1)) ) / ((kend+1) - (kstart-1));
      end

      %%% Continue the search after the end of this group of days
      k = kend;

    end

    %%% Increment day index
    k = k + 1;

  end

  writeDataset(forcingvar1,fullfile(inputfr,pressure),ieee,prec);

  clear forcingvar1

end









