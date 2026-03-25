
%%% Script calculates JJA mean salt flux within defined ronne polynya
%%% region for each offshore wind experiment
%%% (ronnregion.m) 


%%% Frequency of diagnostic output

%%% Read experiment data
setExpname
loadexp;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
deltaT = 440;
% nIter0 = 0;
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

deltaT_4 = 200;
nIter0_4 = 1;
nDumps_4 = round(nTimeSteps*deltaT_4/dumpFreq);
dumpIters_4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters_4 = dumpIters_4(dumpIters_4 >= nIter0_4);

deltaT2 = 400;
nDumps2 = round(nTimeSteps*deltaT2/dumpFreq);
dumpIters_2 = round((1:nDumps2)*dumpFreq/deltaT2); 
dumpIters_2 = dumpIters_2(dumpIters_2 >= nIter0);
nDumps2 = length(dumpIters_2);


deltaT_3 = 300;
nIter0_3 = 1;
nDumps_3 = round(nTimeSteps*10*deltaT_3/dumpFreq);
dumpIters_3 = round((1:nDumps_3)*dumpFreq/deltaT_3); 
dumpIters_3 = dumpIters_3(dumpIters_3 >= nIter0_3);

deltaT_5 = 480;
nIter0_5 = 1;
nDumps_5 = round(nTimeSteps*10*deltaT_5/dumpFreq);
dumpIters_5 = round((1:nDumps_5)*dumpFreq/deltaT_5); 
dumpIters_5 = dumpIters_5(dumpIters_5 >= nIter0_5);



exppath1 = '/data3/MITgcm_WS/experiments/s_hist38_c_minus50';
exppath2 = '/data3/MITgcm_WS/experiments/s_hist31_c_meridplus10';
exppath3 = '/data3/MITgcm_WS/experiments/s_hist1_meridplus20';
exppath4 = '/data3/MITgcm_WS/experiments/s_hist32_c_meridplus5';
exppath5 = '/data3/MITgcm_WS/experiments/n_34452';
exppath6 = '/data3/MITgcm_WS/experiments/s_hist37_c_minus15_2';
exppath7 = '/data3/MITgcm_WS/experiments/s_hist34_cminus5_2';
exppath8 = '/data3/MITgcm_WS/experiments/s_hist19_cminus2_2';
exppath9 = '/data3/MITgcm_WS/experiments/s_hist36_c_minus20_2';
exppath10 = '/data3/MITgcm_WS/experiments/s_hist35_cminus10_2';
exppath11 = '/data3/MITgcm_WS/experiments/minus30_conth';
exppath12 = '/data3/MITgcm_WS/experiments/s_hist41_c_minus40';
exppath13 = '/data3/MITgcm_WS/experiments/s_hist45_c_meridplus15';



exppath1w = '/data3/MITgcm_WS/experiments/n_342';
exppath2w = '/data3/MITgcm_WS/experiments/s_hist4_warm_meridplus20';
exppath3w = '/data3/MITgcm_WS/experiments/s_hist18_warm_windsplus2';
exppath4w = '/data3/MITgcm_WS/experiments/w_hist42_wminus50';
exppath5w = '/data3/MITgcm_WS/experiments/s_hist40_w_plus10';
exppath6w = '/data3/MITgcm_WS/experiments/s_hist24_w_mplus15_2';
exppath7w = '/data3/MITgcm_WS/experiments/s_hist32_w_meridminus10';
exppath8w = '/data3/MITgcm_WS/experiments/s_hist32_w_meridminus20';
exppath9w = '/data3/MITgcm_WS/experiments/s_hist32_w_meridminus5';
exppath10w = '/data3/MITgcm_WS/experiments/s_hist24_w_mplus15_2';
exppath11w = '/data3/MITgcm_WS/experiments/s_hist43_w_meridplus5';
exppath12w = '/data3/MITgcm_WS/experiments/s_hist44_w_meridminus30';
exppath13w =  '/data3/MITgcm_WS/experiments/s_hist46_w_minus40';


ronneregion

%%%% finding bottom
kmax = ones(Nx,Ny);
kmin = ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    idx = find(squeeze(hFacC(i,j,:))>0);
    if (~isempty(idx))
      kmin(i,j) = min(idx);
      kmax(i,j) = max(idx);
    end
  end
end


%%%%%%%%%%%%%% Warm Experiments First

% $###########Set time range
% tmin = 15*86400*360;
% tmax = 20*86400*360;
Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Nyears = Nyear_f-Nyear_s;

Cavity_salt_control_warm = NaN(Nx,Ny,nDump_Finish-nDump_start);
%%% n = 1:nDumps (but experiment hasn't finished yet %%%
for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath1w,'/results/SFLUX'),dumpIters_4(n+nDump_start));
  Cavity_salt_control_warm(:,:,n) = (SALT); %%%convert to kg/m2/s
  tt(n) =  (dumpIters(n)*deltaT)/86400;
  %tt(n) is the day of the run
  
end

SALT_season = zeros(size(Cavity_salt_control_warm,1),size(Cavity_salt_control_warm,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_control_warm(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_control_warm(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_control_warm(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_control_warm(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp1w = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));


salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_1= 0
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_1(1) = SALT_winter_tot_1 +(salt_winter(i,j)*RAC(i,j));
                            
                end
        end
            
end

%%%%%%%%%%%%%
% $###########Set time range

Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Nyears = Nyear_f-Nyear_s;


% Cavity_temp_wplus5 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_salt_wplus5 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath11w,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_wplus5(:,:,n) = (SALT);
  
  
end

SALT_season = zeros(size(Cavity_salt_wplus5,1),size(Cavity_salt_wplus5,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = nanmean(Cavity_salt_wplus5(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = nanmean(Cavity_salt_wplus5(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = nanmean(Cavity_salt_wplus5(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = nanmean(Cavity_salt_wplus5(:,:,(n-1)*12+(9:11)),3);
  
end
M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp11w = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));


salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_wplus5=0;
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_wplus5 = SALT_winter_tot_wplus5+ +(salt_winter(i,j)*RAC(i,j));%convert to kg /yr
                            
                end
        end
            
end


% SALT_winter_tot_wplus5=nanmean(SALT_winter_tot_wplus5,2);
% SALT_winter_tot_wplus5=nanmean(SALT_winter_tot_wplus5);
% 



%%%%%%%%%%%
%%%%%NEXT
Nyear_s =9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;
Nyears = Nyear_f-Nyear_s;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_salt_wplus10 = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath5w,'/results/SFLUX'),dumpIters_4(n+nDump_start));
  Cavity_salt_wplus10(:,:,n) = (SALT);
  
end

SALT_season = zeros(size(Cavity_salt_wplus10,1),size(Cavity_salt_wplus10,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_wplus10(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_wplus10(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_wplus10(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_wplus10(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp5w = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));


salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_wplus10= 0;
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_wplus10 = SALT_winter_tot_wplus10+ (salt_winter(i,j)*RAC(i,j));
                          
                end
        end
            
end


%%%%%%%%%%%%%%

%%%%%NEXT
Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Cavity_salt_wplus20 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Nyears = Nyear_f-Nyear_s;

for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath2w,'/results/SFLUX'),dumpIters_2(n+nDump_start));
  Cavity_salt_wplus20(:,:,n) = (SALT);
  
end


SALT_season = zeros(size(Cavity_salt_wplus20,1),size(Cavity_salt_wplus20,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_wplus20(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_wplus20(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_wplus20(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_wplus20(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp2w = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));


salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_wplus20=0;
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1
                                    SALT_winter_tot_wplus20 = SALT_winter_tot_wplus20+ (salt_winter(i,j)*RAC(i,j));
                             
                end
        end
            
end





%%%%%NEXT
Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;
Cavity_salt_wminus10 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Nyears = Nyear_f-Nyear_s;

for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath7w,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_wminus10(:,:,n) = (SALT);
  
end


SALT_season = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_wminus10(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_wminus10(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_wminus10(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_wminus10(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp7w = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));

salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_wminus10=0;
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_wminus10 = SALT_winter_tot_wminus10+ (salt_winter(i,j)*RAC(i,j));%%% convert to kg /yr
                             
                end
        end
            
end


%%%%%NEXT
% Nyear_s = 2;
% Nyear_f = 7;
% nDump_start = Nyear_s * 12;
% nDump_Finish = Nyear_f * 12;
Cavity_salt_wminus20 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Nyears = Nyear_f-Nyear_s;

for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath8w,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_wminus20(:,:,n) = (SALT);
  
end


SALT_season = zeros(size(Cavity_salt_wminus20,1),size(Cavity_salt_wminus20,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_wminus20(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_wminus20(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_wminus20(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_wminus20(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp8w = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));


salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_wminus20=0;
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_wminus20 = SALT_winter_tot_wminus20+ (salt_winter(i,j)*RAC(i,j));
                         
                end
        end
            
end



%%%%%NEXT
Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Nyear_s_2 = 9;
Nyear_f_2 = 18;
nDump_start2 = Nyear_s_2 * 12;
nDump_Finish2 = Nyear_f_2 * 12;

Nyear_s_3 = 19;
Nyear_f_3 = 27;
nDump_start3 = Nyear_s_3 * 12;
nDump_Finish3= Nyear_f_3 * 12;


Cavity_salt_wminus50 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_salt_wminus40 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_salt_wminus30 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_salt_wplus15 = NaN(Nx,Ny,nDump_Finish3-nDump_start3);
Cavity_salt_wminus5= NaN(Nx,Ny,nDump_Finish-nDump_start);
Nyears = Nyear_f-Nyear_s;

for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath4w,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_wminus50(:,:,n) = (SALT);
  
  SALT40 = rdmdsWrapper(fullfile(exppath13w,'/results/SFLUX'),dumpIters_3(n+nDump_start2));
  Cavity_salt_wminus40(:,:,n) = (SALT40); 
  
  SALT30 = rdmdsWrapper(fullfile(exppath12w,'/results/SFLUX'),dumpIters_3(n+nDump_start2));
  Cavity_salt_wminus30(:,:,n) = (SALT30);
  
  SALT5 = rdmdsWrapper(fullfile(exppath9w,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_wminus5(:,:,n) = (SALT5); 
end
for n=1:(nDump_Finish3-nDump_start3)
  SALT15 = rdmdsWrapper(fullfile(exppath6w,'/results/SFLUX'),dumpIters_4(n+nDump_start3+2));
  Cavity_salt_wplus15(:,:,n) = (SALT15);
end
Nyears2=7;
  


SALT_season = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);
SALT_season40 = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);
SALT_season30 = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);
SALT_season15 = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears2,4);
SALT_season5 = zeros(size(Cavity_salt_wminus5,1),size(Cavity_salt_wminus10,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_wminus50(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_wminus50(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_wminus50(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_wminus50(:,:,(n-1)*12+(9:11)),3);
 
  SALT_season40(:,:,n,1) = mean(Cavity_salt_wminus40(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season40(:,:,n,2) = mean(Cavity_salt_wminus40(:,:,(n-1)*12+(3:5)),3);
  SALT_season40(:,:,n,3) = mean(Cavity_salt_wminus40(:,:,(n-1)*12+(6:8)),3);
  SALT_season40(:,:,n,4) = mean(Cavity_salt_wminus40(:,:,(n-1)*12+(9:11)),3);
  
  SALT_season30(:,:,n,1) = mean(Cavity_salt_wminus30(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season30(:,:,n,2) = mean(Cavity_salt_wminus30(:,:,(n-1)*12+(3:5)),3);
  SALT_season30(:,:,n,3) = mean(Cavity_salt_wminus30(:,:,(n-1)*12+(6:8)),3);
  SALT_season30(:,:,n,4) = mean(Cavity_salt_wminus30(:,:,(n-1)*12+(9:11)),3); 
  
  SALT_season5(:,:,n,1) = mean(Cavity_salt_wminus5(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season5(:,:,n,2) = mean(Cavity_salt_wminus5(:,:,(n-1)*12+(3:5)),3);
  SALT_season5(:,:,n,3) = mean(Cavity_salt_wminus5(:,:,(n-1)*12+(6:8)),3);
  SALT_season5(:,:,n,4) = mean(Cavity_salt_wminus5(:,:,(n-1)*12+(9:11)),3); 
end



for n = 1:Nyears2
  SALT_season15(:,:,n,1) = mean(Cavity_salt_wplus15(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season15(:,:,n,2) = mean(Cavity_salt_wplus15(:,:,(n-1)*12+(3:5)),3);
  SALT_season15(:,:,n,3) = mean(Cavity_salt_wplus15(:,:,(n-1)*12+(6:8)),3);
  SALT_season15(:,:,n,4) = mean(Cavity_salt_wplus15(:,:,(n-1)*12+(9:11)),3);  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

M4 = squeeze(mean(SALT_season40,1));
N4 = squeeze(mean(M4,1));
M3 = squeeze(mean(SALT_season30,1));
N3 = squeeze(mean(M3,1));
M15 = squeeze(mean(SALT_season40,1));
N15 = squeeze(mean(M15,1));
M5 = squeeze(mean(SALT_season5,1));
N5 = squeeze(mean(M5,1));
err_exp4w = std(N(:,3));
err_exp13w = std(N4(:,3));
err_exp12w = std(N3(:,3));
err_exp6w = std(N15(:,3));
err_exp9w = std(N5(:,3));




SALT_seasonavg = squeeze(mean(SALT_season,3));
SALT_seasonavg40 = squeeze(mean(SALT_season40,3));
SALT_seasonavg30 = squeeze(mean(SALT_season30,3));
SALT_seasonavg15 = squeeze(mean(SALT_season15,3));
SALT_seasonavg5 = squeeze(mean(SALT_season5,3));


salt_winter = SALT_seasonavg(:,:,3);
salt_winter40 = SALT_seasonavg40(:,:,3);
salt_winter30 = SALT_seasonavg30(:,:,3);
salt_winter15 = SALT_seasonavg15(:,:,3);
salt_winter5 = SALT_seasonavg5(:,:,3);

%%%%%%% Finding average grid cell salinity
SALT_winter_tot_wminus50=0;
SALT_winter_tot_wminus40=0;
SALT_winter_tot_wminus30=0;
SALT_winter_tot_wplus15=0;
SALT_winter_tot_wminus5=0;

for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_wminus50 = SALT_winter_tot_wminus50+ (salt_winter(i,j)*RAC(i,j)); %%% convert to kg /yr
                                    SALT_winter_tot_wminus40 = SALT_winter_tot_wminus40+ (salt_winter40(i,j)*RAC(i,j)); %%% convert to kg /yr
                                    SALT_winter_tot_wminus30 = SALT_winter_tot_wminus30+ (salt_winter30(i,j)*RAC(i,j)); %%% convert to kg /yr
                                    SALT_winter_tot_wplus15 = SALT_winter_tot_wplus15+ (salt_winter15(i,j)*RAC(i,j)); %%% convert to kg /yr
                                    SALT_winter_tot_wminus5 = SALT_winter_tot_wminus5+ (salt_winter5(i,j)*RAC(i,j)); %%% convert to kg /yr

                           
                end
        end
            
end


%%%%%%%%%%%%%%%%%%%%
%%%Control

%%%%%NEXT
Nyear_s = 11;
Nyear_f = 20;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;
Cavity_salt_control = NaN(Nx,Ny,nDump_Finish-nDump_start);
Nyears = Nyear_f-Nyear_s;

for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath5,'/results/SFLUX'),dumpIters(n+nDump_start));
  Cavity_salt_control(:,:,n) = (SALT);
  
end


SALT_season = zeros(size(Cavity_salt_control,1),size(Cavity_salt_control,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_control(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_control(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_control(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_control(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp5 = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));


s_mean = 0;
salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_control=0;
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_control = SALT_winter_tot_control+ (salt_winter(i,j)*RAC(i,j));

                end
        end
            
end

%%%%%%%%%%%
%%%%%NEXT
Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;
Nyears = Nyear_f-Nyear_s;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_salt_cplus10 = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath2,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_cplus10(:,:,n) = (SALT);
  
end

SALT_season = zeros(size(Cavity_salt_cplus10,1),size(Cavity_salt_cplus10,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_cplus10(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_cplus10(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_cplus10(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_cplus10(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp2 = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));


salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_cplus10= 0
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_cplus10 = SALT_winter_tot_cplus10+ (salt_winter(i,j)*RAC(i,j));
            
                end
        end
            
end



%%%%%%%%%%%
%%%%%NEXT
Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;
Nyears = Nyear_f-Nyear_s;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_salt_cplus20 = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath3,'/results/SFLUX'),dumpIters_4(n+nDump_start));
  Cavity_salt_cplus20(:,:,n) = (SALT);
  
end

SALT_season = zeros(size(Cavity_salt_cplus20,1),size(Cavity_salt_cplus20,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_cplus20(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_cplus20(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_cplus20(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_cplus20(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp3 = std(N(:,3));

SALT_seasonavg = squeeze(mean(SALT_season,3));


salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_cplus20= 0
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_cplus20 = SALT_winter_tot_cplus20+ (salt_winter(i,j)*RAC(i,j));
           
                end
        end
            
end

%%%%%%%%%%%
%%%%%NEXT
Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Nyear_s2 = 18;
Nyear_f2 = 27;
nDump_start2 = Nyear_s2 * 12;
nDump_Finish2 = Nyear_f2 * 12;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_salt_cminus50 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_salt_cminus40 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_salt_cminus30 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_salt_cplus15 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_salt_cplus5 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_salt_cminus5 = NaN(Nx,Ny,nDump_Finish-nDump_start);

Nyears = Nyear_f-Nyear_s;

for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath1,'/results/SFLUX'),dumpIters_3(n+nDump_start2));
  Cavity_salt_cminus50(:,:,n) = (SALT);
  
  SALT40 = rdmdsWrapper(fullfile(exppath12,'/results/SFLUX'),dumpIters_2(n+nDump_start2));
  Cavity_salt_cminus40(:,:,n) = (SALT40); 
  
  SALT30 = rdmdsWrapper(fullfile(exppath11,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_cminus30(:,:,n) = (SALT30);
  
  SALT15 = rdmdsWrapper(fullfile(exppath13,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_cplus15(:,:,n) = (SALT15);
  
  SALTm5 = rdmdsWrapper(fullfile(exppath7,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_cminus5(:,:,n) = (SALTm5);
    
  SALTp5 = rdmdsWrapper(fullfile(exppath4,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_cplus5(:,:,n) = (SALTp5);
  
end


SALT_season = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);
SALT_season40 = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);
SALT_season30 = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);
SALT_season15 = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);
SALT_seasonp5 = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);
SALT_seasonm5 = zeros(size(Cavity_salt_wminus10,1),size(Cavity_salt_wminus10,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_cminus50(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_cminus50(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_cminus50(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_cminus50(:,:,(n-1)*12+(9:11)),3);
 
  SALT_season40(:,:,n,1) = mean(Cavity_salt_cminus40(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season40(:,:,n,2) = mean(Cavity_salt_cminus40(:,:,(n-1)*12+(3:5)),3);
  SALT_season40(:,:,n,3) = mean(Cavity_salt_cminus40(:,:,(n-1)*12+(6:8)),3);
  SALT_season40(:,:,n,4) = mean(Cavity_salt_cminus40(:,:,(n-1)*12+(9:11)),3);
  
  SALT_season30(:,:,n,1) = mean(Cavity_salt_cminus30(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season30(:,:,n,2) = mean(Cavity_salt_cminus30(:,:,(n-1)*12+(3:5)),3);
  SALT_season30(:,:,n,3) = mean(Cavity_salt_cminus30(:,:,(n-1)*12+(6:8)),3);
  SALT_season30(:,:,n,4) = mean(Cavity_salt_cminus30(:,:,(n-1)*12+(9:11)),3);  
  
  SALT_season15(:,:,n,1) = mean(Cavity_salt_cplus15(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season15(:,:,n,2) = mean(Cavity_salt_cplus15(:,:,(n-1)*12+(3:5)),3);
  SALT_season15(:,:,n,3) = mean(Cavity_salt_cplus15(:,:,(n-1)*12+(6:8)),3);
  SALT_season15(:,:,n,4) = mean(Cavity_salt_cplus15(:,:,(n-1)*12+(9:11)),3);  
  
  SALT_seasonm5(:,:,n,1) = mean(Cavity_salt_cminus5(:,:,(n-1)*12+[12 1 2]),3);
  SALT_seasonm5(:,:,n,2) = mean(Cavity_salt_cminus5(:,:,(n-1)*12+(3:5)),3);
  SALT_seasonm5(:,:,n,3) = mean(Cavity_salt_cminus5(:,:,(n-1)*12+(6:8)),3);
  SALT_seasonm5(:,:,n,4) = mean(Cavity_salt_cminus5(:,:,(n-1)*12+(9:11)),3);    
  
  SALT_seasonp5(:,:,n,1) = mean(Cavity_salt_cplus5(:,:,(n-1)*12+[12 1 2]),3);
  SALT_seasonp5(:,:,n,2) = mean(Cavity_salt_cplus5(:,:,(n-1)*12+(3:5)),3);
  SALT_seasonp5(:,:,n,3) = mean(Cavity_salt_cplus5(:,:,(n-1)*12+(6:8)),3);
  SALT_seasonp5(:,:,n,4) = mean(Cavity_salt_cplus5(:,:,(n-1)*12+(9:11)),3);   
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

M4 = squeeze(mean(SALT_season40,1));
N4 = squeeze(mean(M4,1));
M3 = squeeze(mean(SALT_season30,1));
N3 = squeeze(mean(M3,1));
M15 = squeeze(mean(SALT_season40,1));
N15 = squeeze(mean(M15,1));

Mm5 = squeeze(mean(SALT_seasonm5,1));
Nm5 = squeeze(mean(Mm5,1));
M5 = squeeze(mean(SALT_seasonp5,1));
N5 = squeeze(mean(M5,1));
err_exp1 = std(N(:,3));
err_exp12 = std(N4(:,3));
err_exp11 = std(N3(:,3));
err_exp13 = std(N15(:,3));
err_exp7 = std(Nm5(:,3));
err_exp4 = std(N5(:,3));



SALT_seasonavg = squeeze(mean(SALT_season,3));
SALT_seasonavg40 = squeeze(mean(SALT_season40,3));
SALT_seasonavg30 = squeeze(mean(SALT_season30,3));
SALT_seasonavg15 = squeeze(mean(SALT_season15,3));
SALT_seasonavgm5 = squeeze(mean(SALT_seasonm5,3));
SALT_seasonavgp5 = squeeze(mean(SALT_seasonp5,3));


salt_winter = SALT_seasonavg(:,:,3);
salt_winter40 = SALT_seasonavg40(:,:,3);
salt_winter30 = SALT_seasonavg30(:,:,3);
salt_winter15 = SALT_seasonavg15(:,:,3);
salt_winterp5 = SALT_seasonavgp5(:,:,3);
salt_winterm5 = SALT_seasonavgm5(:,:,3);

%%%%%%% Finding average grid cell salinity
SALT_winter_tot_cminus50=0;
SALT_winter_tot_cminus40=0;
SALT_winter_tot_cminus30=0;
SALT_winter_tot_cplus15=0;
SALT_winter_tot_cplus5=0;
SALT_winter_tot_cminus5=0;

for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_cminus50 = SALT_winter_tot_cminus50+ (salt_winter(i,j)*RAC(i,j)); %%% convert to kg /yr
                                    SALT_winter_tot_cminus40 = SALT_winter_tot_cminus40+ (salt_winter40(i,j)*RAC(i,j)); %%% convert to kg /yr
                                    SALT_winter_tot_cminus30 = SALT_winter_tot_cminus30+ (salt_winter30(i,j)*RAC(i,j)); %%% convert to kg /yr
                                    SALT_winter_tot_cplus15 = SALT_winter_tot_cplus15+ (salt_winter15(i,j)*RAC(i,j)); %%% convert to kg /yr
                                    SALT_winter_tot_cminus5 = SALT_winter_tot_cminus5+ (salt_winterm5(i,j)*RAC(i,j)); %%% convert to kg /yr
                                    SALT_winter_tot_cplus5 = SALT_winter_tot_cplus5+ (salt_winterp5(i,j)*RAC(i,j)); %%% convert to kg /yr

                end
        end
            
end


%%%%%%%%%%%
%%%%%NEXT
Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;
Nyears = Nyear_f-Nyear_s;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_salt_cminus20 = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath9,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_cminus20(:,:,n) = (SALT);
  
end

SALT_season = zeros(size(Cavity_salt_cminus20,1),size(Cavity_salt_cminus20,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_cminus20(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_cminus20(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_cminus20(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_cminus20(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp9 = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));


salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_cminus20= 0
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_cminus20 = SALT_winter_tot_cminus20+ (salt_winter(i,j)*RAC(i,j));
 
                end
        end
            
end


%%%%%%%%%%%
%%%%%NEXT
Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;
Nyears = Nyear_f-Nyear_s;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_salt_cminus10 = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SALT = rdmdsWrapper(fullfile(exppath10,'/results/SFLUX'),dumpIters_3(n+nDump_start));
  Cavity_salt_cminus10(:,:,n) = (SALT);
  
end

SALT_season = zeros(size(Cavity_salt_cminus10,1),size(Cavity_salt_cminus10,2),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,n,1) = mean(Cavity_salt_cminus10(:,:,(n-1)*12+[12 1 2]),3);
  SALT_season(:,:,n,2) = mean(Cavity_salt_cminus10(:,:,(n-1)*12+(3:5)),3);
  SALT_season(:,:,n,3) = mean(Cavity_salt_cminus10(:,:,(n-1)*12+(6:8)),3);
  SALT_season(:,:,n,4) = mean(Cavity_salt_cminus10(:,:,(n-1)*12+(9:11)),3);
  
end

M = squeeze(mean(SALT_season,1));
N = squeeze(mean(M,1));

err_exp10 = std(N(:,3));
SALT_seasonavg = squeeze(mean(SALT_season,3));

salt_winter = SALT_seasonavg(:,:,3);
%%%%%%% Finding average grid cell salinity
SALT_winter_tot_cminus10= 0
for i = 1:Nx
        for j = 1:Ny
           
                if start(i,j)==1


                                    SALT_winter_tot_cminus10 = SALT_winter_tot_cminus10+ (salt_winter(i,j)*RAC(i,j));
               
                end
        end
            
end


figure(1)
clf
    

datanewwarm = [SALT_winter_tot_wminus50,SALT_winter_tot_cminus50,SALT_winter_tot_wminus40,SALT_winter_tot_cminus40, SALT_winter_tot_wminus30,SALT_winter_tot_wminus20,SALT_winter_tot_wminus10,SALT_winter_tot_wminus5,SALT_winter_tot_1,SALT_winter_tot_wplus5,SALT_winter_tot_wplus10];
w_hssw=[.5,.5,.6,.6,.7,.8,.9,.95,1,1.05,1.1];



datanewcold= [ SALT_winter_tot_cminus30,SALT_winter_tot_cminus20,SALT_winter_tot_cminus10,SALT_winter_tot_cminus5,SALT_winter_tot_control,SALT_winter_tot_cplus5,SALT_winter_tot_cplus10,SALT_winter_tot_cplus15,SALT_winter_tot_wplus15,SALT_winter_tot_wplus20,SALT_winter_tot_cplus20];
w_chssw=[.7,.8,.9,.95,1,1.05,1.1,1.15,1.15,1.2,1.2];

combo = horzcat(datanewwarm,datanewcold);
combostrength=horzcat(w_hssw,w_chssw);

dataw = [SALT_winter_tot_wminus50,SALT_winter_tot_wminus40,SALT_winter_tot_wminus30,SALT_winter_tot_wminus20,SALT_winter_tot_wminus10, SALT_winter_tot_wminus5,SALT_winter_tot_1,SALT_winter_tot_wplus5,SALT_winter_tot_wplus10,SALT_winter_tot_wplus15,SALT_winter_tot_wplus20];
datac = [SALT_winter_tot_cminus50,SALT_winter_tot_cminus40, SALT_winter_tot_cminus30,SALT_winter_tot_cminus20,SALT_winter_tot_cminus10, SALT_winter_tot_cminus5,SALT_winter_tot_control,SALT_winter_tot_cplus5,SALT_winter_tot_cplus10,SALT_winter_tot_cplus15,SALT_winter_tot_cplus20];

hold on



%%%%plot warm 
% w = plot(w_s,dataw,'-or','LineWidth',2);
hold on
errorbar(w_hssw(1),datanewwarm(1),err_exp4w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_hssw(2),datanewwarm(2),err_exp13w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_hssw(3),datanewwarm(3),err_exp12w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_hssw(4),datanewwarm(4),err_exp8w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_hssw(5),datanewwarm(5),err_exp7w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_hssw(6),datanewwarm(6),err_exp9w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_hssw(7),datanewwarm(7),err_exp1w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_hssw(8),datanewwarm(8),err_exp11w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on

errorbar(w_hssw(9),datanewwarm(9),err_exp5w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on

errorbar(w_hssw(10),datanewwarm(10),err_exp6w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on

errorbar(w_hssw(11),datanewwarm(11),err_exp2w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on




%%%% plot cold
% c = plot(w_c,datac,'-ob','LineWidth',2);
hold on
errorbar(w_chssw(1),datanewcold(1),err_exp1,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_chssw(2),datanewcold(2),err_exp12,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');

hold on
errorbar(w_chssw(3),datanewcold(3),err_exp11,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_chssw(4),datanewcold(4),err_exp9,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_chssw(5),datanewcold(5),err_exp10,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_chssw(6),datanewcold(6),err_exp7,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_chssw(7),datanewcold(7),err_exp5,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_chssw(8),datanewcold(8),err_exp4,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_chssw(9),datanewcold(9),err_exp2,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_chssw(10),datanewcold(10),err_exp13,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_chssw(11),datanewcold(11),err_exp3,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on



% hold on
% errorbar(w_s(1),dataw(1),err_exp1,'-ks','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black');
% hold on
% errorbar(w_s(2),dataw(2),err_exp2,'-ks','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black');
% hold on
% errorbar(w_s(3),dataw(3),err_exp3,'-ks','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black');
% hold on
% errorbar(w_s(4),dataw(4),err_exp4,'-ks','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black');
% 



title('JJA Integrated Salt Flux, Ronne Polynya (g/s)','Fontsize',20,'interpreter','latex')
xlabel('Wind Perturbation ($\chi$)','Fontsize',16,'interpreter','latex')
ylabel('Integrated Salt Flux (g/s)','Fontsize',16,'interpreter','latex')

% ylim([.5e5 2e5]);





hhssw = fitlm(w_hssw,datanewwarm);
p=polyfit(w_hssw,datanewwarm,1);
f = polyval(p,w_hssw);
hold on
plot(w_hssw,f,'--r')

hc = fitlm(w_chssw,datanewcold);
pw=polyfit(w_chssw,datanewcold,1);
fw = polyval(pw,w_chssw);
hold on
plot(w_chssw,fw,'--b')

tots = fitlm(combostrength,combo);

% hc = fitlm(nnc,datac);
% pc=polyfit(nnc,datac,1);
% fc = polyval(pc,nnc);
% hold on
% plot(nnc,fc,'--b')
text(1.05,2.5e8,'R$^2$=.96','fontsize',18,'color','r','interpreter','latex');
text(1.05,1.5e8,'R$^2$=.99','fontsize',18,'color','b','interpreter','latex');
% text(1.05,1.5e8,'R$^2$=.98','fontsize',18,'color','b','interpreter','latex');



% text(.6,11e12,'$\sum_{p}$ = C$_{p}$V$_{wind}$ + $\sum_{p0}$','fontsize',15,'color','k','interpreter','latex');







%%%const = -1.4e8







