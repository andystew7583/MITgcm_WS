%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Daily Surface Forcing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defineGrid
%%%%%%%%%%% Adding path for humidity and temperature


addpath(fullfile(gendir,'/MITgcm_WS/newexp/convert_humidity/convert_humidity'));

addpath(fullfile(gendir,'MITgcm_WS/newexp/Temperature/Temperature'));


%%%%%%%%%%%%% load data

datadir = fullfile(gendir,'/MITgcm_WS/data/PolarWrf');


%%%%%%%%% Highest resolution of AMPS

AMPS_x = 666;
AMPS_y = 627;
AMPS_dataspan = 3287;




%%%%%%%%%% alpha
%%% choose amount by which to perturb winds.

alpha = 1;

%%% scale will be (x(i)-min(x)/max(x)-minx) * alpha
%%%


%%% Surface air temperature perturbation

deltaT = 0;




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
  dpm(i,9) = 31;
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
datadirSWDN = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/SWDN/');
datadirLWDN = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/LWDN/');
datadirTEMP = fullfile(gendir,'/MITgcm_WS/data/PolarWrf/dailydata/Temp/');
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
          
           %%% Determines whether to append to an existing file
            append = days_cntr > 1;             
            
%            %%% NC file names depend on calendar date
            varname_u = 'UGRD_10m';       
            ncname_u = 'RotatedU';
% 
            
%             %%% Calculate next zonal wind daily average 
%             [Uwindinterpolated,F_U] = WRFdailyAvgwinds(years(i),months(j),days,datadiru,varname_u,ncname_u,F_U,ERR_U); 
%             Uwindinterpolated(isnan(Uwindinterpolated))=0;
% %             
% %             %%% Adding Wind Perturbation
%             for k=1:Nx
%                 for p = 1:Ny
% 
%                         Uwindinterpolated(k,p) = (Uwindinterpolated(k,p)*alpha);
%                         if Uwindinterpolated(k,p)<-20
%                             Uwindinterpolated(k,p)=-20;
%                         end
%                         if Uwindinterpolated(k,p)>20
%                             Uwindinterpolated(k,p)=20;
%                         end                   
%                 end
%             end
% 
%             
%             
%              
%             %%% NC file names depend on calendar date
%             varname_v = 'VGRD_10m'; 
%             ncname_v = 'RotatedV';
% 
%             
%             %%% Calculate next zonal wind daily average                       
%             [Vwindinterpolated,F_V] = WRFdailyAvgwinds(years(i),months(j),days,datadirv,varname_v,ncname_v,F_V,ERR_V);
%             Vwindinterpolated(isnan(Vwindinterpolated)) = 0;
%             
%                         %%% Adding Wind Perturbation
%             for k=1:Nx
%                 for p = 1:Ny
% 
%                         Vwindinterpolated(k,p) =  (Vwindinterpolated(k,p)*alpha);
%                         if Vwindinterpolated(k,p)<-20
%                             Vwindinterpolated(k,p)=-20;
%                         end
%                         if Vwindinterpolated(k,p)>20
%                             Vwindinterpolated(k,p)=20;
%                         end                        
%                         
%                 end
%             end
            
% % % %
% %%             [Pressure9years,F_PRES] = WRFdailyAvg(years(i),months(j),days,datadirPRES,'PRES','PRES',F_PRES,ERR_PRES);
% %%             Pressure9years(isnan(Pressure9years)) = 0;
% %%             
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

                %%% Calculate next zonal wind daily average
                  [Temp,F_TEMP] = WRFdailyAvg(years(i),months(j),days,datadirTEMP,varname_t,ncname_t,F_TEMP,ERR_TEMP);
                  Temp(isnan(Temp)) = 0;
              
%             
% %             
%            % %% NC file names depend on calendar date
%             varname_rh = 'RH_2m'; 
%             if (years(i)<2008) || (years(i)==2008 && months(j) <=10)             
%               ncname_rh = 'RH_2m';
%             elseif (years(i) ==2008 &&months(j) >11) || years(i) <2013
%               ncname_rh = 'RH';
%             elseif (years(i)==2013 && months(j)==1)
%               ncname_rh = 'RH_2m';
%             elseif (years(i) ==2013 &&months(j) >1) || years(i) >2013
%               ncname_rh = 'RH';
%             end
%             
%             %%% Calculate next zonal wind daily average 
%             [Humid,F_RH] = WRFdailyAvg(years(i),months(j),days,datadirRH,varname_rh,ncname_rh,F_RH,ERR_HUMID);           
%             Humid9years = convert_humidity(Pressure9years,Temp,Humid,'relative humidity','specific humidity');          
%             
%             Humid9years(isnan(Humid9years)) = 0;
%             
            
            
            %%%%%%%%%%%%%%%%%%%%%%
            %Converting Temp    %%
            %%%%%%%%%%%%%%%%%%%%%%
 
 
            Temp9years = kel2cel(Temp);
            Temp9years = Temp9years+deltaT;
%             
            
%              
%             [Precip,F_PCP] = WRFdailyAvg(years(i),months(j),days,datadirPCP,'PCP','PCP',F_PCP,0);
%             Precip(isnan(Precip)) = 0;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%Converting Precip  (Currently in kg/m^2 (mm) need m/s  %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% First get it into mm/hr (it is 3hr accumulated precip) 
%             Precip = Precip/3;
%             %%% Conversion to m/s 
%             Precip9years = Precip*.000000278;         
%             
%             
%              
%              
%             [LW9years,F_LWDN] = WRFdailyAvg(years(i),months(j),days,datadirLWDN,'LWDN','LWDN',F_LWDN,ERR_LW);
%             LW9years(isnan(LW9years)) = 0;
%               
%             [SW9years,F_SWDN] = WRFdailyAvg(years(i),months(j),days,datadirSWDN,'SWDN','SWDN',F_SWDN,ERR_SW);
%             SW9years(isnan(SW9years)) = 0; 

             
             

             

            
             

            %%% Add to output file
             addpath ../newexp_utils
%              writeDatasetA(Uwindinterpolated,fullfile(inputfolder,zwind),ieee,prec,append);
%              
%              writeDatasetA(Vwindinterpolated,fullfile(inputfolder,mwind),ieee,prec,append);
%              
             writeDatasetA(Temp9years,fullfile(inputfolder,aTemp),ieee,prec,append);
%              
%              writeDatasetA(Precip9years,fullfile(inputfolder,aPrecip),ieee,prec,append);
%              
%              writeDatasetA(LW9years,fullfile(inputfolder,aLW),ieee,prec,append);
%                  
%              writeDatasetA(SW9years,fullfile(inputfolder,aSW),ieee,prec,append);
             
%              writeDatasetA(Pressure9years,fullfile(inputfolder,pressure),ieee,prec,append);
% 
%              writeDatasetA(Humid9years,fullfile(inputfolder,anewAQ),ieee,prec,append);

             
  
             %%% Increment day counter
             days_cntr = days_cntr + 1     
        
        end
       
    end
 end
 


