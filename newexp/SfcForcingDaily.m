%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Daily Surface Forcing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defineGrid


%%% Options
load_wind = true;
load_pressure = false;
load_temp = false;
load_rh = false;
load_sw = false;
load_lw = false;
load_precip = false;

%%% Need pressure and temp to convert humidity
if (load_rh)
  load_pressure = true;
  load_temp = true;
end

%%%%%%%%%%% Adding path for humidity and temperature


addpath(fullfile(gendir,'/MITgcm_WS/newexp/convert_humidity/convert_humidity'));

addpath(fullfile(gendir,'MITgcm_WS/newexp/Temperature/Temperature'));


%%% Path to utilities for writing binary files

addpath ../newexp_utils


%%%%%%%%%%%%% load data

datadir = fullfile(gendir,'/MITgcm_WS/data/PolarWrf');






%%% Surface air temperature perturbation

deltaT = 0;



%%% Max componentwise wind speed
wind_max = 30;


%%%%%% Calculating dpm for each month, using leap year condiitonal
dpm = zeros(Nyears,Nmon);
for i=1:Nyears
  dpm(i,1) = 31;
  if (is_leap_year(i))
    dpm(i,2) = 29;
  else
    
    dpm(i,2) = 28;
  end
  dpm(i,3) = 31;
  dpm(i,4) = 30;
  dpm(i,5) = 31;
  dpm(i,6) = 30;
  dpm(i,7) = 31;
  dpm(i,8) = 31;
  dpm(i,9) = 30;
  dpm(i,10) = 31;
  dpm(i,11) = 30;
  dpm(i,12) = 31;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Creating Uwind/Vwind/Pressure/Precip/RAD/Intital Files %%%%%

datadiru = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/RotatedUwind/');
datadirv = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/RotatedVWind/');
datadirRH = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/RHdaily/');
datadirPCP = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/PCP/');
% datadirSWDN = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/SWDN/');
datadirSWDN = fullfile(gendir,'/MITgcm_WS/data/ERAInterim/SfcForcing/SW/');
datadirLWDN = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/LWDN/');
% datadirLWDN = fullfile(gendir,'/MITgcm_WS/data/ERAInterim/SfcForcing/LW/');
datadirTEMP = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/Temp/');
% datadirTEMP = fullfile(gendir,'/MITgcm_WS/data/ERAInterim/SfcForcing/2mTemp/');
datadirPRES = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/PRES/');

F_U = [];
F_V = [];
F_RH = [];
F_PCP = [];
F_SWDN = [];
F_LWDN = [];
F_TEMP = [];
F_PRES = [];

ERR_U = 0;
ERR_V = 0;
ERR_PRES = 0;
ERR_TEMP = -273.15;
ERR_HUMID = 0;
ERR_LW = 0;
ERR_PRECIP = 0;
ERR_SW = 0;



years = start_year:1:endyr;
months = start_month:1:end_month;

days_cntr = 1;

 for i=1:Nyears
    
    for j= 1:Nmon
        
        for days = 1:dpm(i,j)
          
          disp([num2str(years(i)),'-',num2str(months(j)),'-',num2str(days)]);
          
           %%% Determines whether to append to an existing file
            append = days_cntr > 1;             
            
            
            
            
            
            
            
            
            
            
            
            if (load_wind)
            
            %%% UWIND
              
             %%% NC file names depend on calendar date
              varname_u = 'UGRD_10m';       
              ncname_u = 'RotatedU';
  

              %%% Calculate next zonal wind daily average 
              [Uwindinterpolated,F_U] = WRFdailyAvgwinds(years(i),months(j),days,datadiru,varname_u,ncname_u,F_U,ERR_U,EXF_XMC,EXF_YMC); 
              Uwindinterpolated(isnan(Uwindinterpolated))=0;
              
              %%% VWIND
             
              %%% NC file names depend on calendar date
              varname_v = 'VGRD_10m'; 
              ncname_v = 'RotatedV';


              %%% Calculate next zonal wind daily average                       
              [Vwindinterpolated,F_V] = WRFdailyAvgwinds(years(i),months(j),days,datadirv,varname_v,ncname_v,F_V,ERR_V,EXF_XMC,EXF_YMC);
              Vwindinterpolated(isnan(Vwindinterpolated)) = 0;
              
              %%% Limit wind magnitudes
              Uwindinterpolated(Uwindinterpolated>wind_max) = wind_max;
              Uwindinterpolated(Uwindinterpolated<-wind_max) = -wind_max;
              Vwindinterpolated(Vwindinterpolated>wind_max) = wind_max;
              Vwindinterpolated(Vwindinterpolated<-wind_max) = -wind_max;
              
              %%% Append to file
              writeDatasetA(Uwindinterpolated,fullfile(inputfolder,zwind),ieee,prec,append);            
              writeDatasetA(Vwindinterpolated,fullfile(inputfolder,mwind),ieee,prec,append);
              
            end

            
            
            
            
            

            
            %%% PRESSURE
            if (load_pressure)

              [Pressure9years,F_PRES] = WRFdailyAvg(years(i),months(j),days,datadirPRES,'PRES','PRES',F_PRES,ERR_PRES,EXF_XMC,EXF_YMC);
              Pressure9years(isnan(Pressure9years)) = 0;
              
             writeDatasetA(Pressure9years,fullfile(inputfolder,pressure),ieee,prec,append);
              
            end

            
            
            %%% TEMPERATURE
            if (load_temp)

             % %% NC file names depend on calendar date
              varname_t = 'TMP_2m';
              if (years(i)<2008) || (years(i)==2008 && months(j) <=10)
                ncname_t = 'TMP_2m';
              elseif (years(i) ==2008 && months(j) >11) || years(i) <2013
                ncname_t = 'TMP';
              elseif (years(i)==2013 && months(j)==1)
                ncname_t = 'TMP_2m';
              elseif (years(i) ==2013 && months(j) >1) || years(i)>2013
                ncname_t = 'TMP';
              end

            %%% Calculate next daily average
              [Temp,F_TEMP] = WRFdailyAvg(years(i),months(j),days,datadirTEMP,varname_t,ncname_t,F_TEMP,ERR_TEMP,EXF_XMC,EXF_YMC);
              Temp(isnan(Temp)) = 0;
              
%               Temp = ERAdailyAvg (years(i),months(j),days,datadirTEMP,'ei.oper.an.sfc.regn128sc.2d_168','2D_GDS4_SFC','hazel323879',EXF_XMC,EXF_YMC,false);

              %%%%%%%%%%%%%%%%%%%%%%
              %Converting Temp    %%
              %%%%%%%%%%%%%%%%%%%%%%

              Temp9years = kel2cel(Temp);
              Temp9years = Temp9years+deltaT;              
              
              writeDatasetA(Temp9years,fullfile(inputfolder,aTemp),ieee,prec,append);          
              
            end
              
            
            
            %%% RELATIVE HUMIDITY
            if (load_rh)

             % %% NC file names depend on calendar date
              varname_rh = 'RH_2m'; 
              if (years(i)<2008) || (years(i)==2008 && months(j) <=10)             
                ncname_rh = 'RH_2m';
              elseif (years(i) ==2008 &&months(j) >11) || years(i) <2013
                ncname_rh = 'RH';
              elseif (years(i)==2013 && months(j)==1)
                ncname_rh = 'RH_2m';
              elseif (years(i) ==2013 &&months(j) >1) || years(i) >2013
                ncname_rh = 'RH';
              end

              %%% Calculate next zonal wind daily average 
              [Humid,F_RH] = WRFdailyAvg(years(i),months(j),days,datadirRH,varname_rh,ncname_rh,F_RH,ERR_HUMID,EXF_XMC,EXF_YMC);           
              Humid9years = convert_humidity(Pressure9years,Temp,Humid,'relative humidity','specific humidity');          
              Humid9years(isnan(Humid9years)) = 0;
              
              writeDatasetA(Humid9years,fullfile(inputfolder,anewAQ),ieee,prec,append);
            
            end
            
            
            
            
            
            
            
            
            %%% PRECIPITATION
            if (load_precip)

              [Precip,F_PCP] = WRFdailyAvg(years(i),months(j),days,datadirPCP,'PCP','PCP',F_PCP,0,EXF_XMC,EXF_YMC);
              Precip(isnan(Precip)) = 0;


              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %%%Converting Precip  (Currently in kg/m^2 (mm) need m/s  %%%
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

              %%% First get it into mm/hr (it is 3hr accumulated precip) 
              Precip = Precip/3;
              %%% Conversion to m/s 
              Precip9years = Precip*.000000278;         
            
              writeDatasetA(Precip9years,fullfile(inputfolder,aPrecip),ieee,prec,append);
             
            end            
             
            
            
            
            
            %%% LONGWAVE
            if (load_lw)
             
              [LW9years,F_LWDN] = WRFdailyAvg(years(i),months(j),days,datadirLWDN,'LWDN','LWDN',F_LWDN,ERR_LW,EXF_XMC,EXF_YMC);
              LW9years(isnan(LW9years)) = 0;
             
%               LW9years = ERAdailyAvg (years(i),months(j),days,datadirLWDN,'ei.oper.fc.sfc.regn128sc.strd_175','STRD_GDS4_SFC','hazel323880',EXF_XMC,EXF_YMC,true);
              
              writeDatasetA(LW9years,fullfile(inputfolder,aLW),ieee,prec,append);

            end
            
            
            
            
            
            %%% SHORTWAVE
            if (load_sw)

%               [SW9years,F_SWDN] = WRFdailyAvg(years(i),months(j),days,datadirSWDN,'SWDN','SWDN',F_SWDN,ERR_SW,EXF_XMC,EXF_YMC);
%               SW9years(isnan(SW9years)) = 0; 
              
              SW9years = ERAdailyAvg (years(i),months(j),days,datadirSWDN,'ei.oper.fc.sfc.regn128sc.ssrd_169','SSRD_GDS4_SFC','hazel323880',EXF_XMC,EXF_YMC,true);
              
              writeDatasetA(SW9years,fullfile(inputfolder,aSW),ieee,prec,append);

            end
             
             
  
             %%% Increment day counter
             days_cntr = days_cntr + 1     
        
        end
       
    end
 end
 


