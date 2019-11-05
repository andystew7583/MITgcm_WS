%%%% Checking surface forcing for ERAInterim (LW,SW,Temp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Daily Surface Forcing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../newexp_utils

defineGrid
%%%%%%%%%%% Adding path for humidity and temperature
  load('XMC');
  load('YMC');
  
Weddell_Lon = XMC';
Weddell_Lat = YMC';
lat_ERA=ncread('/data3/MITgcm_WS/data/ERAInterim/SfcForcing/2mTemp/ei.oper.an.sfc.regn128sc.2d_168.201807.hazel323879.nc','g4_lat_1');
lon_ERA=ncread('/data3/MITgcm_WS/data/ERAInterim/SfcForcing/2mTemp/ei.oper.an.sfc.regn128sc.2d_168.201807.hazel323879.nc','g4_lon_2'); 

addpath(fullfile(gendir,'MITgcm_WS/newexp/Temperature/Temperature'));

%%%%%%%%%%%%% load data

datadir = fullfile(gendir,'/MITgcm_WS/data/ERAInterim');


%%%%%%%%% Highest resolution of AMPS

ERA_x = 512;
ERA_y = 256;

lon_ERA=lon_ERA';
lat_ERA = lat_ERA';
lat_ERA=flip(lat_ERA);
idx1 = find(lon_ERA>180);
idx2 = find(lon_ERA<=180);
new_lon = [lon_ERA(idx1)-360 lon_ERA(idx2)];
[LO,LA]=meshgrid(new_lon,lat_ERA);

%%% Beginning Data from 2007  ( first full year we have)

base_year = 2006;
start_year = 2007;
endyr = 2016;

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

% alpha = 1.2;

%%% scale will be (x(i)-min(x)/max(x)-minx) * alpha
%%%

ncname_lw = 'STRD_GDS4_SFC';
ncname_sw = 'SSRD_GDS4_SFC';
ncname_temp = '2D_GDS4_SFC';


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

datadirSWDN = fullfile(gendir,'/MITgcm_WS/data/ERAInterim/SfcForcing/SW/');
datadirLWDN = fullfile(gendir,'/MITgcm_WS/data/ERAInterim/SfcForcing/LW/');
datadirTEMP = fullfile(gendir,'/MITgcm_WS/data/ERAInterim/SfcForcing/2mTemp/');


sw_daily = zeros(Nx,Ny,366,Nyears);
lw_daily = zeros(Nx,Ny,366,Nyears);
temp_daily=zeros(Nx,Ny,366,Nyears);

sw_daily_2 = zeros(Nx,Ny,31,Nmon,Nyears);
lw_daily_2 = zeros(Nx,Ny,31,Nmon,Nyears);
temp_daily_2 =zeros(Nx,Ny,31,Nmon,Nyears);

years = start_year:1:endyr;
months = start_month:1:end_month;

days_cntr = 1;

 for i=1:Nyears
    
    for j= 1:Nmon
        
        
           %%% Determines whether to append to an existing file
            append = days_cntr > 1;             
           if years(i)<2011 
            datafname_temp = fullfile(datadirTEMP,['ei.oper.an.sfc.regn128sc.2d_168.', num2str(years(i)),'.hazel323879.nc']);
            datafname_lw = fullfile(datadirLWDN,['ei.oper.fc.sfc.regn128sc.strd_175.', num2str(years(i)),'.hazel323880.nc']);
            datafname_sw = fullfile(datadirSWDN,['ei.oper.fc.sfc.regn128sc.ssrd_169.', num2str(years(i)),'.hazel323880.nc']);

            if (exist(datafname_temp))==2
                    data_temp = ncread(datafname_temp,ncname_temp);
                    dpy = size(data_temp,3);
                    data_temp = reshape(data_temp,ERA_x,ERA_y,4,dpy/4);
                    data_tempnew = squeeze(nanmean(data_temp,3));
                    data_tempnew = data_tempnew([idx1 idx2],:,:);
                    data_tempnew = flip(data_tempnew,2);
                    for day = 1:size(data_tempnew,3)
                        avgdata_t = kel2cel(data_tempnew(:,:,day)');
                        temp_daily_one = interp2(LO,LA,data_tempnew(:,:,day)',XMC',YMC','spline');
                        temp_daily(:,:,day,i)=temp_daily_one;
                    end

            end
            
            if (exist(datafname_lw))==2
                y = 1;
               data_lw = ncread(datafname_lw,ncname_lw);
               data_lw = data_lw([idx1 idx2],:,:);
               data_lw = flip(data_lw,2);

               dpy = size(data_lw,3);
                    for n = 1:2:size(data_lw,3)-1
                        avgdata_l = (data_lw(:,:,(n+1))-data_lw(:,:,(n)))/86400;
                        lw_daily_one = interp2(LO,LA,avgdata_l',XMC,YMC,'spline');
                        lw_daily_one = lw_daily_one';
                        lw_daily(:,:,y,i) = lw_daily_one;
                        y = y+1;

                    end

            end
            if (exist(datafname_sw))==2
                y = 1;
               data_sw = ncread(datafname_sw,ncname_sw);
               data_sw = data_sw([idx1 idx2],:,:);
               data_sw = flip(data_sw,2);
            
               dpy = size(data_sw,3);
                    for n = 1:2:size(data_sw,3)-1
                        avgdata_sw = (data_sw(:,:,(n+1))-data_sw(:,:,(n)))/86400;
                        sw_daily_one = interp2(LO,LA,avgdata_sw',XMC,YMC,'spline');
                        sw_daily(:,:,y,i)=sw_daily_one';
                        y = y+1;

                    end

            end  
            clear y
           elseif years(i)>=2011
            datafname_temp = fullfile(datadirTEMP,['ei.oper.an.sfc.regn128sc.2d_168.', num2str(years(i)) num2str(sprintf('%02d',months(j))),'.hazel323879.nc']);
            datafname_lw = fullfile(datadirLWDN,['ei.oper.fc.sfc.regn128sc.strd_175.', num2str(years(i)) num2str(sprintf('%02d',months(j))),'.hazel323880.nc']);
            datafname_sw = fullfile(datadirSWDN,['ei.oper.fc.sfc.regn128sc.ssrd_169.', num2str(years(i)) num2str(sprintf('%02d',months(j))),'.hazel323880.nc']);

            if (exist(datafname_temp))==2
                    data_temp = ncread(datafname_temp,ncname_temp);
                    dpy = size(data_temp,3);
                    data_temp = reshape(data_temp,ERA_x,ERA_y,4,dpy/4);
                    data_tempnew = squeeze(nanmean(data_temp,3));
                    data_tempnew = data_tempnew([idx1 idx2],:,:);
                    data_tempnew = flip(data_tempnew,2);
             
                    for day = 1:size(data_tempnew,3)
                        avgdata_t = kel2cel(data_tempnew(:,:,day)');
                        temp_daily_one = interp2(LO,LA,avgdata_t,XMC',YMC','spline');
                        temp_daily_2(:,:,day,j,i)=temp_daily_one;
                    end

            end
            
            if (exist(datafname_lw))==2
                y = 1
               data_lw = ncread(datafname_lw,ncname_lw);
               data_lw = data_lw([idx1 idx2],:,:);
               data_lw = flip(data_lw,2);
             
               dpy = size(data_lw,3);
                    for n = 1:2:size(data_lw,3)-1
                        avgdata_l = (data_lw(:,:,(n+1))-data_lw(:,:,(n)))/86400
                        lw_daily_one = interp2(LO,LA,avgdata_l',XMC,YMC,'spline')
                        lw_daily_2(:,:,y,j,i) = lw_daily_one';
                        y = y+1;
                      
                    end

            end
            if (exist(datafname_sw))==2
               data_sw = ncread(datafname_sw,ncname_sw);
               data_sw = data_sw([idx1 idx2],:,:);
               data_sw = flip(data_sw,2);
             
               dpy = size(data_sw,3);
               y = 1;
                    for n = 1:2:size(data_sw,3)-1
                        avgdata_sw = (data_sw(:,:,(n+1))-data_sw(:,:,(n)))/86400;
                        sw_daily_one = interp2(LO,LA,avgdata_sw',XMC,YMC,'nearest');
                        sw_daily_2(:,:,y,j,i)=sw_daily_one';
                        y = y+1;
                        
                    end

            end                  
           end
           
      
    end
 end
 

%   temp_daily = squeeze(nanmean(temp_daily,4));
 temp_daily = reshape(temp_daily,Nx,Ny,366*9);
 s=squeeze(temp_daily(1,1,:));
 index = find(s==0);
 temp_daily(:,:,index)=[];
 clear index s
 
 
 temp_daily_2=reshape(temp_daily_2,Nx,Ny,31*9*12);
 s=squeeze(temp_daily(1,1,:));
 index = find(s==0);
 temp_daily_2(:,:,index)=[];

  
 finaltemp=cat(3,temp_daily,temp_daily_2);
  
  
 %   temp_daily = squeeze(nanmean(temp_daily,4));
 
 lw_daily = reshape(lw_daily,Nx,Ny,366*9);
 p=squeeze(lw_daily(1,1,:));
 index = find(p==0);
 lw_daily(:,:,index)=[];
 
 
 lw_daily_2=reshape(lw_daily_2,Nx,Ny,31*9*12);
 s=squeeze(lw_daily(1,1,:));
 index2 = find(s==0);
 lw_daily_2(:,:,index2)=[];
  finallw=cat(3,lw_daily,lw_daily_2);

 
 
 sw_daily = reshape(sw_daily,Nx,Ny,366*9);
 sw_daily(:,:,index)=[];
 
 
 sw_daily_2=reshape(sw_daily_2,Nx,Ny,31*9*12);
 sw_daily_2(:,:,index2)=[];

  
 finalsw=cat(3,sw_daily,sw_daily_2);
 
 


  
 
%   sw_daily = squeeze(nanmean(sw_daily,4));
%   sw_daily = reshape(sw_daily,Nx,Ny,366*9);    
%   
%   lw_daily = squeeze(nanmean(lw_daily,4));
%   lw_daily = reshape(lw_daily,Nx,Ny,366*9);                
%              
  
  writeDataset(finaltemp,fullfile(inputfolder,aTemp),ieee,prec);

  writeDataset(finalsw,fullfile(inputfolder,aSW),ieee,prec);
  
  writeDataset(finallw,fullfile(inputfolder,aLW),ieee,prec);


  
