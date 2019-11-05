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

%%% Beginning Data from 2007  ( first full year we have)

base_year = 2006;
start_year = 2007;
endyr = 2015;

start_month = 1;
end_month = 12;

%%%%  total years, months

Nmon = 12;
Nyears = 9;
Nmonths = 108;


%%%% Finding if current year is a leap year
dpm = zeros(Nyears,Nmon);
is_leap_year = zeros(Nyears,1);
for i=1:Nyears
  if (mod(i+2,4)==0)
    is_leap_year(i) = 1;
  end
end



%%%%%%%%%% alpha
%%% choose amount by which to perturb winds.

% alpha = 1.5;

%%% scale will be (x(i)-min(x)/max(x)-minx) * alpha
%%%


% Alpha = zeros(Nx,Ny);
% max_lon = 20;
% min_lon = 0;
% interval = .0025;
% % interval = .0085;
% alpha = 1;
% for i = 1:Nx
%     if xmc(i) < 0
%        Alpha(i,:) = 0;
%     else
%        Alpha(i,:) = alpha;
%        alpha = alpha +interval;
%     end
%     
% end




%%%%%% Calculating dpm for each month, using leap year condiitonal





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


years = start_year:1:endyr;
months = start_month:1:end_month;  



days_cntr = 1;
time_cntr = 1;

 for i=8:Nyears

        
    for j= 1:Nmon
        
        if (years(i)<2008) || (years(i)==2008 && months(j) <=10)
            starttime = 0;
            endtime = 21;
        end
        if (years(i)==2008 && months(j) >=11) || (years(i) >=2009)

            starttime = 6;
            endtime = 27;
        end
            
        
        
        for days = 1:dpm(i,j)
            
            for time = starttime:03:endtime
          
                %%% Determines whether to append to an existing file
                append = time_cntr > 1;             
            
%                 %%% NC file names depend on calendar date
                varname_u = 'UGRD_10m';       
                ncname_u = 'RotatedU';

%             
%                 %%% Calculate next zonal wind daily average 
                [Uwindinterpolated,F_U] = WRF3hrwinds(years(i),months(j),days,time,datadiru,varname_u,ncname_u,F_U); 
                 Uwindinterpolated(isnan(Uwindinterpolated))=0;
            
                %%%% Adding Wind Perturbation
%             for k=1:Nx
%                 for p = 1:Ny
%                     if Alpha(k,p)==0
%                            %disp('hello')
%                         Uwindinterpolated(k,p) = Uwindinterpolated(k,p);
%                     elseif Alpha(k,p) > 0
%                           %disp('hellno')
%                         Uwindinterpolated(k,p) = Uwindinterpolated(k,p)*Alpha(k,p);
%                     end
%                 end
%             end

            
            
             
                %%% NC file names depend on calendar date
                varname_v = 'VGRD_10m'; 
                ncname_v = 'RotatedV';

            
                %%% Calculate next zonal wind daily average                       
                [Vwindinterpolated,F_V] = WRF3hrwinds(years(i),months(j),days,time,datadirv,varname_v,ncname_v,F_V);
                 Vwindinterpolated(isnan(Vwindinterpolated)) = 0;

                [Pressure9years,F_PRES] = WRF3hour(years(i),months(j),days,time,datadirPRES,'PRES','PRES',F_PRES);
                 Pressure9years(isnan(Pressure9years)) = 0;
            
            %%% NC file names depend on calendar date
                 varname_t = 'TMP_2m'; 
                 if (years(i)<2008) || (years(i)==2008 && months(j) <=10)             
                  ncname_t = 'TMP_2m';
                 elseif (years(i) ==2008 && months(j) >11) || years(i) <2013
                  ncname_t = 'TMP';
                 elseif (years(i)==2013 && months(j)==1);
                  ncname_t = 'TMP_2m';
                 elseif (years(i) ==2013 && months(j) >1) ;
                  ncname_t = 'TMP';
                 end
            
            %%% Calculate next zonal wind daily average 
                [Temp,F_TEMP] = WRF3hour(years(i),months(j),days,time,datadirTEMP,varname_t,ncname_t,F_TEMP); 
                 Temp(isnan(Temp)) = 0;
                     
            
            %%% NC file names depend on calendar date
                varname_rh = 'RH_2m'; 
                if (years(i)<2008) || (years(i)==2008 && months(j) <=10)             
                 ncname_rh = 'RH_2m';
                elseif (years(i) ==2008 &&months(j) >11) || (years(i) <2013)
                 ncname_rh = 'RH';
                elseif (years(i)==2013 && months(j)==1)
                 ncname_rh = 'RH_2m';
                elseif (years(i) ==2013 &&months(j) >1) 
                 ncname_rh = 'RH';
                end
            
            %%% Calculate next zonal wind daily average 
                [Humid,F_RH] = WRF3hour(years(i),months(j),days,time,datadirRH,varname_rh,ncname_rh,F_RH);           
                Humid9years = convert_humidity(Pressure9years,Temp,Humid,'relative humidity','specific humidity');          
            
                Humid9years(isnan(Humid9years)) = 0;
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%Converting Temp    %%
            %%%%%%%%%%%%%%%%%%%%%%%%%
 
 
                Temp9years = kel2cel(Temp);
%                 %Temp9years
%             
%              
                [Precip,F_PCP] = WRF3hour(years(i),months(j),days,time,datadirPCP,'PCP','PCP',F_PCP);
                Precip(isnan(Precip)) = 0;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%Converting Precip  (Currently in kg/m^2 (mm) need m/s  %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% First get it into mm/hr (it is 3hr accumulated precip) 
                Precip = Precip/3;
            %%% Conversion to m/s 
                Precip9years = Precip*.000000278;         
            
            
             
             
                [LW9years,F_LWDN] = WRF3hour(years(i),months(j),days,time,datadirLWDN,'LWDN','LWDN',F_LWDN);
                 LW9years(isnan(LW9years)) = 0;
             
                [SW9years,F_SWDN] = WRF3hour(years(i),months(j),days,time,datadirSWDN,'SWDN','SWDN',F_SWDN);
                 SW9years(isnan(SW9years)) = 0; 
%              

             
             

             

            
             

             %%% Add to output file
                addpath ../newexp_utils
%                 writeDatasetA(Uwindinterpolated,fullfile(inputfolder,zwind),ieee,prec,append);
%              
%                 writeDatasetA(Vwindinterpolated,fullfile(inputfolder,mwind),ieee,prec,append);
%              
%                 writeDatasetA(Temp9years,fullfile(inputfolder,aTemp),ieee,prec,append);
%              
%                 writeDatasetA(Precip9years,fullfile(inputfolder,aPrecip),ieee,prec,append);
%              
%                 writeDatasetA(LW9years,fullfile(inputfolder,aLW),ieee,prec,append);
%                  
%                 writeDatasetA(SW9years,fullfile(inputfolder,aSW),ieee,prec,append);
%              
%                 writeDatasetA(Pressure9years,fullfile(inputfolder,'apressurefile.bin'),ieee,prec,append);
% 
%                 writeDatasetA(Humid9years,fullfile(inputfolder,anewAQ),ieee,prec,append);
             time_cntr = time_cntr+1
            end
  
             %%% Increment day counter
             days_cntr = days_cntr + 1     
        
        end
       
    end
 end
 

