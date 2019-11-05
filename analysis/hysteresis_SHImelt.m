%%%% plot hysteresis shimelt

%%% Read experiment data
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
dumpIters4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters4 = dumpIters4(dumpIters4 >= nIter0_4);

deltaT2 = 400;
nDumps2 = round(nTimeSteps*deltaT2/dumpFreq);
dumpIters_orig = round((1:nDumps2)*dumpFreq/deltaT2); 
dumpIters2 = dumpIters_orig(dumpIters_orig >= nIter0);
nDumps2 = length(dumpIters2);


deltaT_3 = 300;
nIter0_3 = 1;
nDumps_3 = round(nTimeSteps*deltaT_3*10/dumpFreq);
dumpIters_3 = round((1:nDumps_3)*dumpFreq/deltaT_3); 
dumpIters_3 = dumpIters_3(dumpIters_3 >= nIter0_3);

%%% control

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

%%%%%%%%%%%%%% Warm Experiments First

% 
% tmin = 15*86400*360;
% tmax = 20*86400*360;
Nyear_s = 12;
Nyear_f = 21;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;





% Cavity_temp_control_warm = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
Cavity_shimelt_control_warm = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath1w,'/results/SHIfwFlx'),dumpIters4(n+nDump_start));
  Cavity_shimelt_control_warm(:,:,n) = (SHImelt);
  
  
end

% Cavity_temp_control_warm(:,:,:,1:nDump_start) = [];



SHImelt_control_warm=zeros(1,(nDump_Finish-nDump_start));
    
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                         if YC(i,j)<-75
            
                                  SHImelt_control_warm(1,h)= SHImelt_control_warm(1,h) + Cavity_shimelt_control_warm(i,j,h)*RAC(i,j);
                                     
                        end
                end
            end
    end
end
 
SHImelt_control_warm = -(SHImelt_control_warm*86400*365*1e-12);  % convrt to GT/year;  


err_exp1w = std(SHImelt_control_warm);

SHImelt_control_warm = mean(SHImelt_control_warm);

% %%%%%%%%% NEXT
% 
Nyear_s =19;
Nyear_f = 27;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;



% Cavity_temp_wplus5 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_shimelt_wplus15 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath10w,'/results/SHIfwFlx'),dumpIters4(n+nDump_start));
  Cavity_shimelt_wplus15(:,:,n) = (SHImelt);
  
  
end


SHImelt_control_wplus15=zeros(1,(nDump_Finish-nDump_start));
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                        if YC(i,j)<-75
                                  SHImelt_control_wplus15(1,h)= SHImelt_control_wplus15(1,h) + Cavity_shimelt_wplus15(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_wplus15 = -(SHImelt_control_wplus15*86400*365*1e-12);  % convrt to GT/year;  

err_exp6w = std(SHImelt_control_wplus15);

SHImelt_control_wplus15 = mean(SHImelt_control_wplus15);

% %%%%%%%%% NEXT
% 
Nyear_s = 11;
Nyear_f = 20;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;



% Cavity_temp_wplus5 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_shimelt_wplus5 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath11w,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_shimelt_wplus5(:,:,n) = (SHImelt);
  
  
end


SHImelt_control_wplus5=zeros(1,(nDump_Finish-nDump_start));
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                         if YC(i,j)<-75
                                  SHImelt_control_wplus5(1,h)= SHImelt_control_wplus5(1,h) + Cavity_shimelt_wplus5(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_wplus5 = -(SHImelt_control_wplus5*86400*365*1e-12);  % convrt to GT/year;  

err_exp11w = std(SHImelt_control_wplus5);

SHImelt_control_wplus5 = mean(SHImelt_control_wplus5);

%%%%%%NEXTTTT





% Cavity_temp_wplus5 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_shimelt_wminus5 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath9w,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_shimelt_wminus5(:,:,n) = (SHImelt);
  
  
end


SHImelt_control_wminus5=zeros(1,(nDump_Finish-nDump_start));
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                         if YC(i,j)<-75
            
                                  SHImelt_control_wminus5(1,h)= SHImelt_control_wminus5(1,h) + Cavity_shimelt_wminus5(i,j,h)*RAC(i,j);
                                     
                        end
                end
            end
    end
end
 
SHImelt_control_wminus5 = -(SHImelt_control_wminus5*86400*365*1e-12);  % convrt to GT/year;  

err_exp9w = std(SHImelt_control_wminus5);

SHImelt_control_wminus5 = mean(SHImelt_control_wminus5);

%%%%%%NEXTTTT

Nyear_s = 11;
Nyear_f = 20;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;


Cavity_shimelt_wminus10 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath7w,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_shimelt_wminus10(:,:,n) = (SHImelt);
  
  
end


SHImelt_control_wminus10=zeros(1,(nDump_Finish-nDump_start));
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                        if YC(i,j)<-75
                                  SHImelt_control_wminus10(1,h)= SHImelt_control_wminus10(1,h) + Cavity_shimelt_wminus10(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_wminus10 = -(SHImelt_control_wminus10*86400*365*1e-12);  % convrt to GT/year;  

err_exp7w = std(SHImelt_control_wminus10);

SHImelt_control_wminus10 = mean(SHImelt_control_wminus10);

%%%%%%NEXTTTT

Nyear_s = 11;
Nyear_f = 20;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;



% Cavity_temp_wplus5 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_shimelt_wminus20 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath8w,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_shimelt_wminus20(:,:,n) = (SHImelt);
  
  
end


SHImelt_control_wminus20=zeros(1,(nDump_Finish-nDump_start));
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                       if YC(i,j)<-75
            
                                  SHImelt_control_wminus20(1,h)= SHImelt_control_wminus20(1,h) + Cavity_shimelt_wminus20(i,j,h)*RAC(i,j);
                                        
                        
                       end
                end
            end
    end
end
 
SHImelt_control_wminus20 = -(SHImelt_control_wminus20*86400*365*1e-12);  % convrt to GT/year;  

err_exp8w = std(SHImelt_control_wminus20);

SHImelt_control_wminus20 = mean(SHImelt_control_wminus20);




%%%%%NEXT
Nyear_s = 11;
Nyear_f = 20;
nDump_start= Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Cavity_shimelt_wplus10 = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath5w,'/results/SHIfwFlx'),dumpIters4(n+nDump_start-1));
  Cavity_shimelt_wplus10(:,:,n) = (SHImelt);
  
end

% Cavity_temp_wplus10(:,:,:,1:nDump_start) = [];


SHImelt_control_wplus10=zeros(1,nDump_Finish-nDump_start);
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                        if YC(i,j)<-75            
                                  SHImelt_control_wplus10(1,h)= SHImelt_control_wplus10(1,h) + Cavity_shimelt_wplus10(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
 
SHImelt_control_wplus10 = -(SHImelt_control_wplus10*86400*365*1e-12);  % convrt to GT/year;  

err_exp5w = std(SHImelt_control_wplus10);

SHImelt_control_wplus10 = mean(SHImelt_control_wplus10);

%%%%%%%% NEXT
%%%%%NEXT
Nyear_s = 28;
Nyear_f = 37;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Cavity_SHImelt_wplus20 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath2w,'/results/SHIfwFlx'),dumpIters2(n+nDump_start-1));
  Cavity_SHImelt_wplus20(:,:,n) = (SHImelt);
  
end
% Cavity_temp_wplus20(1:nDump_start)=[];



SHImelt_control_wplus20=zeros(1,nDump_Finish-nDump_start);
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                        if YC(i,j)<-75            
                                  SHImelt_control_wplus20(1,h)= SHImelt_control_wplus20(1,h) + Cavity_SHImelt_wplus20(i,j,h)*RAC(i,j);
                                         
                        end
                end
            end
    end
end
 
SHImelt_control_wplus20 = -(SHImelt_control_wplus20*86400*365*1e-12);  % convrt to GT/year;  

err_exp2w = std(SHImelt_control_wplus20);

SHImelt_control_wplus20 = mean(SHImelt_control_wplus20);

%%%%%%NEXTTTT

Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;



% Cavity_temp_wplus5 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_shimelt_wminus50 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath4w,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_shimelt_wminus50(:,:,n) = (SHImelt);
  
  
end


SHImelt_control_wminus50=zeros(1,(nDump_Finish-nDump_start));
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                        if YC(i,j)<-75
                                  SHImelt_control_wminus50(1,h)= SHImelt_control_wminus50(1,h) + Cavity_shimelt_wminus50(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_wminus50 = -(SHImelt_control_wminus50*86400*365*1e-12);  % convrt to GT/year;  

err_exp4w = std(SHImelt_control_wminus50);

SHImelt_control_wminus50 = mean(SHImelt_control_wminus50);

%%%% Cold Experiments

Nyear_s = 9;
Nyear_f = 18;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Cavity_SHImelt_control = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath5,'/results/SHIfwFlx'),dumpIters(n+nDump_start-1));
  Cavity_SHImelt_control(:,:,n) = (SHImelt);
  
  
end

% Cavity_temp_control_warm(:,:,:,1:nDump_start) = [];



SHImelt_control=zeros(1,(nDump_Finish-nDump_start));
    
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                         if YC(i,j)<-75
                                  SHImelt_control(1,h)= SHImelt_control(1,h) + Cavity_SHImelt_control(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
 
SHImelt_control = -(SHImelt_control*86400*365*1e-12);  % convrt to GT/year;  


err_exp5 = std(SHImelt_control);

SHImelt_control = nanmean(SHImelt_control);

%%%%%%%%%%%

Nyear_s4 = 11;
Nyear_f4 = 20;
nDump_start4 = Nyear_s4 * 12;
nDump_Finish4 = Nyear_f4 * 12;


Cavity_SHImelt_minus5 = NaN(Nx,Ny,Nr,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath7,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start4));
  Cavity_SHImelt_minus5(:,:,n) = (SHImelt);
  
  
end

% Cavity_temp_wplus5(:,:,:,1:nDump_start) = [];




SHImelt_control_minus5=zeros(1,(nDump_Finish-nDump_start));
    
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                         if YC(i,j)<-75
                                  SHImelt_control_minus5(1,h)= SHImelt_control_minus5(1,h) + Cavity_SHImelt_minus5(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_minus5 = -(SHImelt_control_minus5*86400*365*1e-12);  % convrt to GT/year;  

err_exp7 = std(SHImelt_control_minus5);

SHImelt_control_minus5 = mean(SHImelt_control_minus5);


%%%%%NEXT
Nyear_s2 = 21;
Nyear_f2 = 30;
nDump_start2 = Nyear_s2 * 12;
nDump_Finish2 = Nyear_f2 * 12;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_SHImelt_minus50 = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath1,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start2-1));
  Cavity_SHImelt_minus50(:,:,n) = (SHImelt);
  
end

% Cavity_temp_wplus10(:,:,:,1:nDump_start) = [];


SHImelt_control_minus50=zeros(1,(nDump_Finish-nDump_start));
    
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                if XC(i,1)>-80 && XC(i,1)<-20
            
                        if YC(1,j) <-74
                                  SHImelt_control_minus50(1,h)= SHImelt_control_minus50(1,h) + Cavity_SHImelt_minus50(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_minus50 = -(SHImelt_control_minus50*86400*365*1e-12);  % convrt to GT/year;  

err_exp1 = std(SHImelt_control_minus50);

SHImelt_control_minus50 = nanmean(SHImelt_control_minus50);



%%%%%NEXT
Nyear_s = 11;
Nyear_f = 20;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_SHImelt_wminus40 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_SHImelt_wminus30 = NaN(Nx,Ny,nDump_Finish-nDump_start);


for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath13w,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start-1));
  Cavity_SHImelt_wminus40(:,:,n) = (SHImelt);
  SHImelt2 = rdmdsWrapper(fullfile(exppath12w,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start-1));
  Cavity_SHImelt_wminus30(:,:,n) = (SHImelt2);
end

% Cavity_temp_wplus10(:,:,:,1:nDump_start) = [];


SHImelt_control_wminus40=zeros(1,(nDump_Finish-nDump_start));
SHImelt_control_wminus30=zeros(1,(nDump_Finish-nDump_start));

 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                         if YC(i,j)<-75
                                  SHImelt_control_wminus40(1,h)= SHImelt_control_wminus40(1,h) + Cavity_SHImelt_wminus40(i,j,h)*RAC(i,j);
                                   SHImelt_control_wminus30(1,h)= SHImelt_control_wminus30(1,h) + Cavity_SHImelt_wminus30(i,j,h)*RAC(i,j);
                                       
                        end
                end
            end
    end
end
 
SHImelt_control_wminus40 = -(SHImelt_control_wminus40*86400*365*1e-12);  % convrt to GT/year;  
SHImelt_control_wminus30 = -(SHImelt_control_wminus30*86400*365*1e-12);  % convrt to GT/year;  

err_exp13w = std(SHImelt_control_wminus40);
err_exp12w = std(SHImelt_control_wminus40);

SHImelt_control_wminus40 = nanmean(SHImelt_control_wminus40);
SHImelt_control_wminus30 = nanmean(SHImelt_control_wminus30);

%%%%%NEXT
Nyear_s = 10;
Nyear_f = 19;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

% Cavity_temp_wplus10 = readIters(exppath,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Cavity_SHImelt_minus10 = NaN(Nx,Ny,nDump_Finish-nDump_start);
for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath10,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start-1));
  Cavity_SHImelt_minus10(:,:,n) = (SHImelt);
  
end

% Cavity_temp_wplus10(:,:,:,1:nDump_start) = [];


SHImelt_control_minus10=zeros(1,(nDump_Finish-nDump_start));
    
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                         if YC(i,j)<-75
                                  SHImelt_control_minus10(1,h)= SHImelt_control_minus10(1,h) + Cavity_SHImelt_minus10(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_minus10 = -(SHImelt_control_minus10*86400*365*1e-12);  % convrt to GT/year;  

err_exp10 = std(SHImelt_control_minus10);

SHImelt_control_minus10 = nanmean(SHImelt_control_minus10);

%%%%%%% NEXT
Nyear_s = 11;
Nyear_f = 20;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Cavity_SHImelt_c_minus15 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath6,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start-1));
  Cavity_SHImelt_c_minus15(:,:,n) = (SHImelt);
  
end



SHImelt_control_minus15=zeros(1,nDump_Finish-nDump_start);
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                         if YC(i,j)<-75
                                  SHImelt_control_minus15(1,h)= SHImelt_control_minus15(1,h) + Cavity_SHImelt_c_minus15(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end

SHImelt_control_minus15 = -(SHImelt_control_minus15*86400*365*1e-12);  % convrt to GT/year;  

err_exp6 = std(SHImelt_control_minus15);

SHImelt_control_minus15 = mean(SHImelt_control_minus15);

%%%%Next

Nyear_s = 11;
Nyear_f = 20;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;


Cavity_SHImelt_minus20 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath9,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_SHImelt_minus20(:,:,n) = (SHImelt);
  
end


SHImelt_control_minus20=zeros(1,nDump_Finish-nDump_start);
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                        if YC(i,j)<-75
                                  SHImelt_control_minus20(1,h)= SHImelt_control_minus20(1,h) + Cavity_SHImelt_minus20(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_minus20 = -(SHImelt_control_minus20*86400*365*1e-12);  % convrt to GT/year;  

err_exp9 = std(SHImelt_control_minus20);

SHImelt_control_minus20 = mean(SHImelt_control_minus20);

%%%%Next

Nyear_s = 11;
Nyear_f = 20;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;


Cavity_SHImelt_plus20 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath3,'/results/SHIfwFlx'),dumpIters4(n+nDump_start));
  Cavity_SHImelt_plus20(:,:,n) = (SHImelt);
  
end


SHImelt_control_plus20=zeros(1,nDump_Finish-nDump_start);
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                        if YC(i,j)<-75
                                  SHImelt_control_plus20(1,h)= SHImelt_control_plus20(1,h) + Cavity_SHImelt_plus20(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_plus20 = -(SHImelt_control_plus20*86400*365*1e-12);  % convrt to GT/year;  

err_exp3 = std(SHImelt_control_plus20);

SHImelt_control_plus20 = mean(SHImelt_control_plus20);

%%%%Next

Nyear_s = 11;
Nyear_f = 20;
nDump_start= Nyear_s * 12;
nDump_Finish = Nyear_f * 12;

Nyear_s2 = 26;
Nyear_f2 = 37;
nDump_start2= Nyear_s2 * 12;
nDump_Finish2 = Nyear_f2 * 12;

Cavity_SHImelt_plus10 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath2,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_SHImelt_plus10(:,:,n) = (SHImelt);
  
end


SHImelt_control_plus10=zeros(1,nDump_Finish-nDump_start);
 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                        if YC(i,j)<-75
                                  SHImelt_control_plus10(1,h)= SHImelt_control_plus10(1,h) + Cavity_SHImelt_plus10(i,j,h)*RAC(i,j);
                                        
                        end
                end
            end
    end
end
 
SHImelt_control_plus10 = -(SHImelt_control_plus10*86400*365*1e-12);  % convrt to GT/year;  

err_exp2 = std(SHImelt_control_plus10);

SHImelt_control_plus10 = mean(SHImelt_control_plus10);


Nyear_s = 11;
Nyear_f = 20;
nDump_start = Nyear_s * 12;
nDump_Finish = Nyear_f * 12;


Cavity_SHImelt_minus30 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_SHImelt_minus40 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_SHImelt_plus5 = NaN(Nx,Ny,nDump_Finish-nDump_start);
Cavity_SHImelt_plus15 = NaN(Nx,Ny,nDump_Finish-nDump_start);

for n=1:(nDump_Finish-nDump_start)
  
  SHImelt = rdmdsWrapper(fullfile(exppath11,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_SHImelt_minus30(:,:,n) = (SHImelt);
  SHImelt40 = rdmdsWrapper(fullfile(exppath12,'/results/SHIfwFlx'),dumpIters2(n+nDump_start2));
  Cavity_SHImelt_minus40(:,:,n) = (SHImelt40);
  SHImeltp5 = rdmdsWrapper(fullfile(exppath4,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_SHImelt_plus5(:,:,n) = (SHImeltp5);
  SHImeltp15 = rdmdsWrapper(fullfile(exppath13,'/results/SHIfwFlx'),dumpIters_3(n+nDump_start));
  Cavity_SHImelt_plus15(:,:,n) = (SHImeltp15); 
  
  
end


SHImelt_control_plus5=zeros(1,nDump_Finish-nDump_start);
SHImelt_control_minus40=zeros(1,nDump_Finish-nDump_start);
SHImelt_control_minus30=zeros(1,nDump_Finish-nDump_start);
SHImelt_control_plus15=zeros(1,nDump_Finish-nDump_start);

 for i=1:Nx
    for j=1:Ny
            for h = 1:(nDump_Finish-nDump_start)

                 if XC(i,j)>-80 && XC(i,j)<-20
                        if YC(i,j)<-75
                                  SHImelt_control_minus30(1,h)= SHImelt_control_minus30(1,h) + Cavity_SHImelt_minus30(i,j,h)*RAC(i,j);
                                  SHImelt_control_minus40(1,h)= SHImelt_control_minus40(1,h) + Cavity_SHImelt_minus40(i,j,h)*RAC(i,j);
                                  SHImelt_control_plus5(1,h)= SHImelt_control_plus5(1,h) + Cavity_SHImelt_plus5(i,j,h)*RAC(i,j);
                                  SHImelt_control_plus15(1,h)= SHImelt_control_plus15(1,h) + Cavity_SHImelt_plus15(i,j,h)*RAC(i,j);
                                     
                        end
                end
            end
    end
end
%  
SHImelt_control_minus30 = -(SHImelt_control_minus30*86400*365*1e-12);  % convrt to GT/year;  
SHImelt_control_minus40 = -(SHImelt_control_minus40*86400*365*1e-12);  % convrt to GT/year;  
SHImelt_control_plus5 = -(SHImelt_control_plus5*86400*365*1e-12);  % convrt to GT/year;  
SHImelt_control_plus15 = -(SHImelt_control_plus15*86400*365*1e-12);  % convrt to GT/year;  


err_exp11 = std(SHImelt_control_minus30);
err_exp12 = std(SHImelt_control_minus40);
err_exp4 = std(SHImelt_control_plus5);
err_exp13 = std(SHImelt_control_plus15);


SHImelt_control_minus30 = mean(SHImelt_control_minus30);
SHImelt_control_minus40 = mean(SHImelt_control_minus40);
SHImelt_control_plus5 = mean(SHImelt_control_plus5);
SHImelt_control_plus15 = mean(SHImelt_control_plus15);



figure(1)
clf

w_c = [.5,.6,.7,.8,.9,.95,1,1.05, 1.1,1.15,1.2];
w_s = [.5,.6,.7,.8,.9,.95,1,1.05 1.1,1.15,1.2];


dataw = [SHImelt_control_wminus50,SHImelt_control_wminus40,SHImelt_control_wminus30,SHImelt_control_wminus20,SHImelt_control_wminus10,SHImelt_control_wminus5, SHImelt_control_warm,SHImelt_control_wplus5,SHImelt_control_wplus10,SHImelt_control_wplus15,SHImelt_control_wplus20];
datac = [SHImelt_control_minus50,SHImelt_control_minus40, SHImelt_control_minus30,SHImelt_control_minus20,SHImelt_control_minus10, SHImelt_control_minus5,SHImelt_control,SHImelt_control_plus5,SHImelt_control_plus10,SHImelt_control_plus15,SHImelt_control_plus20];

hold on

%%%plot observed BML


%%%%plot warm 
w = plot(w_s,dataw,'-or','LineWidth',2);
hold on
errorbar(w_s(1),dataw(1),err_exp4w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(2),dataw(2),err_exp13w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(3),dataw(3),err_exp12w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(4),dataw(4),err_exp8w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(5),dataw(5),err_exp7w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(6),dataw(6),err_exp9w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(7),dataw(7),err_exp1w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(8),dataw(8),err_exp11w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on

errorbar(w_s(9),dataw(9),err_exp5w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on

errorbar(w_s(10),dataw(10),err_exp6w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on

errorbar(w_s(11),dataw(11),err_exp2w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on




%%%% plot cold
c = plot(w_c,datac,'-ob','LineWidth',2);
hold on
errorbar(w_c(1),datac(1),err_exp1,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(2),datac(2),err_exp12,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');

hold on
errorbar(w_c(3),datac(3),err_exp11,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(4),datac(4),err_exp9,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(5),datac(5),err_exp10,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(6),datac(6),err_exp7,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(7),datac(7),err_exp5,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(8),datac(8),err_exp4,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(9),datac(9),err_exp2,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(10),datac(10),err_exp13,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(11),datac(11),err_exp3,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on


t = plot(w_c(7),datac(7),'-sk','MarkerSize',15,'MarkerFaceColor','black');
hold on
q = plot(w_s(7),dataw(7),'-sk','MarkerSize',15,'MarkerFaceColor','black');



% xL=get(gca,'XLim');
% u=line(xL,[82 82],'color','black','linewidth',2);
% u.Marker='.';

% hold on
% axis([.5 1.5 0 1000]);
% hline(100,'k:')
% hline(900,'k:')
% vline(.8,'b:')
% vline(1.23,'r:');

set(gcf,'unit','pixel','position',[0 0 1600 1200])
set(gca,'fontsize',18);

title('Wind Perturbation vs Integrated Shelf Ice Melt (Gt/yr)','interpreter','latex','FontSize',20);
ylabel('Integrated Shelf Ice Melt (Gt/yr)','interpreter','latex','FontSize',16);
xlabel('Wind Perturbation ($\chi$)','interpreter','latex','FontSize',16);
hold on
x=.1:.1:1.7;
y=220;
plot(x,y*ones(size(x)),'k:','linewidth',2)
hold on
% % 
hold on
x=.1:.1:1.7;
y=1100;
plot(x,y*ones(size(x)),'k:','linewidth',2)
hold on



hold on

y = .67*ones(1,1500);
k9 = plot(y,1:1500,'k:','linewidth',2);
hold on
% % 
z = 1.12*ones(1,1500);
k3 = plot(z,1:1500,'k:','linewidth',2);
hold on





legend([w,c],{'Warm State Initialization','Cold State Initialization'},'interpreter','latex','location','northeast');
set(gca,'fontsize',16);



xlim([.5 1.25]);
ylim([0 1200]);

%%%plotting ronne
bathy(SHELFICEtopo-bathy==0)=NaN;
c1=bathy;

% %%%%%
% for i = 1:Nx
%     for j = 1:Ny
%                 if XC(i,j) >-58 && XC(i,j) <-48 && YC(i,j)<-74.5 && YC(i,j)>-76.5
%                          
%                                  bathy(i,j)=0;
%                     
%                              
%                          
%                 end
%     end
% end
% 
% 


 %%%%%
                          
                                 
                                 
                                 
                                 
                                 
