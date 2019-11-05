function [avgdata,F] = WRF3hour (years,months,days,time,griddir,varname,ncname,F)
  %%% get our grid to interpolate onto
  load('XMC6');
  load('YMC6');
  
  Weddell_Lon = XMC';
  Weddell_Lat = YMC';
  
  %%% 2006-10/2008 format
 if (years<2008)||(years==2008 && months <=10)
    %%% Grid for current time period
    LA=ncread('/data3/MITgcm_WS/data/PolarWrf/2006_LAT.nc','LAT');
    LO=ncread('/data3/MITgcm_WS/data/PolarWrf/2006_LON.nc','LON');
         
    %%% Calculate daily average
    avgdata = 0*LA;    
    
    
    datafname = fullfile(griddir,[num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d',days)) num2str(sprintf('%02d',time)),'_WRF_d2_',varname,'.nc']);
      if (exist(datafname))
        avgdata = ncread(datafname,ncname);       
      


      end
    


    
 
 end

                     
    
  
 if (years==2008 && months >=11) || (years >=2009 && years <=2012) 
    %%% Grid for current time period
    LA=ncread('/data3/MITgcm_WS/data/PolarWrf/2008_LAT.nc','LAT');
    LO=ncread('/data3/MITgcm_WS/data/PolarWrf/2008_LON.nc','LON');
    
    %%% Calculate daily average
    avgdata = 0*LA;    

    
    
    datafname = fullfile(griddir,[num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d',days)) num2str(sprintf('%02d',0)),'_WRF_d2_',varname,'_f',num2str(sprintf('%03d',time)),'.nc']);
      if (exist(datafname))
        avgdata =  ncread(datafname,ncname);
        

      end 
    
 end
    




 
                
                                                     
 
                           
 if (years >=2013 && years < 2016)
    %%% Grid for current time period
    LA=ncread('/data3/MITgcm_WS/data/PolarWrf/2013_LAT.nc','LAT');
    LO=ncread('/data3/MITgcm_WS/data/PolarWrf/2013_LON.nc','LON');
    
    %%% Calculate daily average
    avgdata = 0*LA;    
   
    
      datafname = fullfile(griddir,[num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d',days)) num2str(sprintf('%02d',0)),'_WRF_d2_',varname,'_f',num2str(sprintf('%03d',time)),'.nc']);
      if (exist(datafname))
        avgdata =  ncread(datafname,ncname);

      end 
    
 end
    
    
    
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