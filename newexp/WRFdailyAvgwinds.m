function [avgdatawinds,F,iszero] = WRFdailyAvgwinds (years,months,days,griddir,varname,ncname,F,err_val)
  %%% get our grid to interpolate onto
 
  load('XMC');
  load('YMC');
  
  Weddell_Lon = XMC';
  Weddell_Lat = YMC';
  
  

  
  %%% 2006-10/2008 format
  if (years<2008) || (years==2008 && months <=10)        
               
    %%% Grid for current time period
    LA=ncread('/data1/MITgcm_WS/data/PolarWrf/2006_LAT.nc','LAT');
    LO=ncread('/data1/MITgcm_WS/data/PolarWrf/2006_LON.nc','LON');        
         
    %%% Calculate daily average
    avgdatawinds = 0*LA;    
    navg = 0;
    for n=0:3:21
      datafname = fullfile(griddir,[num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d',days)) num2str(sprintf('%02d',n)),varname,'.nc']);
      if (exist(datafname))
        data_temp = ncread(datafname,ncname);
        if ( isempty(err_val) ==0)
          if ( (~isequal(data_temp,err_val*ones(size(data_temp)))) )

                avgdatawinds = avgdatawinds + ncread(datafname,ncname);
                 navg = navg + 1;

          end
        end
      
      end
    end
    
    if (navg > 0)
      avgdatawinds = avgdatawinds / navg;
    else
      iszero = true;
    end
    
                    
  end
                    
    
  
  if (years==2008 && months >=11) || (years >=2009 && years <=2012)
              
    %%% Grid for current time period
    LA=ncread('/data1/MITgcm_WS/data/PolarWrf/2008_LAT.nc','LAT');
    LO=ncread('/data1/MITgcm_WS/data/PolarWrf/2008_LON.nc','LON');
    
    %%% Calculate daily average
    avgdatawinds = 0*LA;    
    navg = 0;
    for n=6:3:27
      datafname = fullfile(griddir,[num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d',days)) num2str(sprintf('%02d',0)),varname,'_f',num2str(sprintf('%03d',n)),'.nc']);
      if (exist(datafname))
        data_temp = ncread(datafname,ncname);
        if ( isempty(err_val) ==0)
          if ( (~isequal(data_temp,err_val*ones(size(data_temp)))) )

                avgdatawinds = avgdatawinds + ncread(datafname,ncname);
                 navg = navg + 1;

          end
        end
      
      end
    end
    
    if (navg > 0)
      avgdatawinds = avgdatawinds / navg;
    else
      iszero = true;
    end

  end
                
                                                       
 
                           
  if (years >=2013)
    %%% Grid for current time period
    LA=ncread('/data1/MITgcm_WS/data/PolarWrf/2013_LAT.nc','LAT');
    LO=ncread('/data1/MITgcm_WS/data/PolarWrf/2013_LON.nc','LON');
    
    %%% Calculate daily average
    avgdatawinds = 0*LA;    
    navg = 0;
    for n=6:3:27
      datafname = fullfile(griddir,[num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d',days)) num2str(sprintf('%02d',0)),varname,'_f',num2str(sprintf('%03d',n)),'.nc']);
      if (exist(datafname))
        data_temp = ncread(datafname,ncname);
        if ( isempty(err_val) ==0)
          if ( (~isequal(data_temp,err_val*ones(size(data_temp)))) )

                avgdatawinds = avgdatawinds + ncread(datafname,ncname);
                 navg = navg + 1;

          end
        end
      
      end
    end
    
    if (navg > 0)
      avgdatawinds = avgdatawinds / navg;
    else
      iszero = true;
    end     
    
  end  
  
  
  
  
  
  
  %%% Reset scattered interpolant object
  if (~isempty(F) && (length(F.Points) ~= size(LA,1)*size(LA,2)))
    F = [];   
  end
 
  %%% Do interpolation
  [avgdatawinds,F] = force(avgdatawinds,LO,LA,Weddell_Lon,Weddell_Lat,F);  
    
end