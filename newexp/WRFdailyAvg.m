%%%
%%% WRFdailyAvg.m
%%%
%%% Helper function that averages AMPS data over a single day and interpolates it onto the model grid.
%%%
function [avgdata,F,iszero] = WRFdailyAvg (years,months,days,griddir,varname,ncname,F,err_val,XMC,YMC)
  
  Weddell_Lon = XMC';
  Weddell_Lat = YMC';
  iszero = false;
  
  %%% 2006-10/2008 format
 if (years<2008)||(years==2008 && months <=10)
    %%% Grid for current time period
    LA=ncread('/data3/MITgcm_WS/data/PolarWrf/2006_LAT.nc','LAT');
    LO=ncread('/data3/MITgcm_WS/data/PolarWrf/2006_LON.nc','LON');
         
    %%% Calculate daily average
    avgdata = 0*LA;    
    navg = 0;

    
    for n=0:3:21
      datafname = fullfile(griddir,[num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d',days)) num2str(sprintf('%02d',n)),'_WRF_d2_',varname,'.nc']);
      if (exist(datafname))==2
        data_temp = ncread(datafname,ncname);
       
        if ( ~isempty(err_val))
          
          
          if ( ~isequal(data_temp,err_val*ones(size(data_temp))))
            
          


             avgdata = avgdata + data_temp;  
      
             navg = navg + 1;

         end
       end
      end
    end
%     
    if navg>0
      avgdata = avgdata/navg;
    else
      iszero = true;
    end
    

 
 end

                     
    
  
 if (years==2008 && months >=11) || (years >=2009 && years <=2012) 
    %%% Grid for current time period
    LA=ncread('/data3/MITgcm_WS/data/PolarWrf/2008_LAT.nc','LAT');
    LO=ncread('/data3/MITgcm_WS/data/PolarWrf/2008_LON.nc','LON');
    
    %%% Calculate daily average
    avgdata = 0*LA;    
    navg = 0;

    
    for n=6:3:27
      datafname = fullfile(griddir,[num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d',days)) num2str(sprintf('%02d',0)),'_WRF_d2_',varname,'_f',num2str(sprintf('%03d',n)),'.nc']);
      if (exist(datafname))==2
        data_temp = ncread(datafname,ncname);
       
        if ( ~isempty(err_val))
          
          
          if ( ~isequal(data_temp,err_val*ones(size(data_temp))))
            
          


             avgdata = avgdata + data_temp;  
      
             navg = navg + 1;

         end
       end
      end
    end
    if navg>0
      avgdata = avgdata/navg;
    else
      iszero = true;
    end
    


 end
 
                
                                                     
 
                           
 if (years >=2013 && years < 2016)
    %%% Grid for current time period
    LA=ncread('/data3/MITgcm_WS/data/PolarWrf/2013_LAT.nc','LAT');
    LO=ncread('/data3/MITgcm_WS/data/PolarWrf/2013_LON.nc','LON');
    
    %%% Calculate daily average
    avgdata = 0*LA;    
    navg = 0;

    
    
    for n=6:3:27
      datafname = fullfile(griddir,[num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d',days)) num2str(sprintf('%02d',0)),'_WRF_d2_',varname,'_f',num2str(sprintf('%03d',n)),'.nc']);
      if (exist(datafname))==2
        data_temp = ncread(datafname,ncname);
       
        if ( ~isempty(err_val))
          
          
          if ( ~isequal(data_temp,err_val*ones(size(data_temp))))
            
          


             avgdata = avgdata + data_temp;  
      
             navg = navg + 1;
  
         end
       end
      end
    end
    if navg>0
      avgdata = avgdata/navg;
    else
      iszero = true;
    end
      

 
 end

  %%% Reset scattered interpolant object
 if (~isempty(F) && (length(F.Points) ~= size(LA,1)*size(LA,2)));
   F = [];
 end
 
 
  %%% Do interpolation
 [avgdata,F] = force(avgdata,LO,LA,Weddell_Lon,Weddell_Lat,F);
 
   
end
  
%                    
%                    
%    
% 
%         
