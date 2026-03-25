%%%%%%%%%%% Plotting Cavity Temperature vs freezing temperature at 1000m

%%%load experiment/setexpname
run ../newexp/defineGrid.m
loadexp


%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
deltaT = 440;
% nIter0 = 0;
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

deltaT4 = 200;
nIter04 = 1;
nDumps_4 = round(nTimeSteps*deltaT4/dumpFreq);
dumpIters4 = round((1:nDumps_4)*dumpFreq/deltaT4); 
dumpIters4 = dumpIters4(dumpIters4 >= nIter04);

deltaT2 = 400;
nDumps2 = round(nTimeSteps*deltaT2/dumpFreq);
dumpIters_orig = round((1:nDumps2)*dumpFreq/deltaT2); 
dumpIters2 = dumpIters_orig(dumpIters_orig >= nIter0);
nDumps2 = length(dumpIters2);


deltaT_3 = 300;
nIter0_3 = 1;
nDumps_3 = round(nTimeSteps*deltaT_3/dumpFreq);
dumpIters_3 = round((1:nDumps_3)*dumpFreq/deltaT_3); 
dumpIters_3 = dumpIters_3(dumpIters_3 >= nIter0_3);




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
exppath11 = '/data3/MITgcm_WS/experiments/s_hist41_c_minus30';
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

depth_level=16;





tmin = 9*86400*360;
tmax = 18*86400*360;


Temp_diff_control = readIters(exppath5,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
Salt_diff_control = readIters(exppath5,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);


% Temp_diff_cminus10 = readIters(exppath10,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_cminus10 = readIters(exppath10,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% 
% Temp_diff_cplus10 = readIters(exppath2,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_cplus10 = readIters(exppath2,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% Temp_diff_cminus20 = readIters(exppath9,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_cminus20 = readIters(exppath9,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% 
% Temp_diff_cplus15 =readIters(exppath13,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_cplus15 =readIters(exppath13,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% 
% Temp_diff_cminus5 = readIters(exppath7,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_cminus5 = readIters(exppath7,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% tmin = 18*86400*360;
% tmax = 27*86400*360;
% 
% Temp_diff_cminus50 =readIters(exppath1,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_cminus50 =readIters(exppath1,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

% tmin = 9*86400*360;
% tmax = 18*86400*360;
% 
% Temp_diff_cplus20 =  readIters(exppath3,'THETA',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_cplus20 =  readIters(exppath3,'SALT',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);
%  
% Temp_diff_cplus5 = readIters(exppath4,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_cplus5 = readIters(exppath4,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% tmin = 27*86400*360;
% tmax = 36*86400*360;
% Temp_diff_cmerid40 = readIters(exppath12,'THETA',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_cmerid40 = readIters(exppath12,'SALT',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
% 



% %%%%% Warm
% tmin = 9*86400*360;
% tmax = 18*86400*360;
% Temp_diff_wminus10 = readIters(exppath7w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_wminus10 = readIters(exppath7w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% Temp_diff_w_plus10 = readIters(exppath5w,'THETA',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_w_plus10 = readIters(exppath5w,'SALT',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);
% 
% 
% Temp_diff_w_minus20 = readIters(exppath8w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_w_minus20 = readIters(exppath8w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

% 
% Temp_diff_wminus5 = readIters(exppath9w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_wminus5 = readIters(exppath9w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% tmin = 27*86400*360;
% tmax = 36*86400*360;
% Temp_diff_w_plus20 = readIters(exppath2w,'THETA',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_w_plus20 = readIters(exppath2w,'SALT',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
% 
% tmin = 19*86400*360;
% tmax = 27*86400*360;
% Temp_diff_w_plus15 = readIters(exppath6w,'THETA',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_w_plus15 = readIters(exppath6w,'SALT',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);
% 
tmin = 9*86400*360;
tmax = 18*86400*360;
% Temp_diff_w_plus5 = readIters(exppath11w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_w_plus5 = readIters(exppath11w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% 
% Temp_diff_wmerid_minus50 = readIters(exppath4w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_wmerid_minus50 = readIters(exppath4w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% 
% Temp_diff_wmerid_minus40 = readIters(exppath13w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% Salt_diff_wmerid_minus40 = readIters(exppath13w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

Temp_diff_w_control = readIters(exppath1w,'THETA',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);
Salt_diff_w_control = readIters(exppath1w,'SALT',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);

Temp_diff_control_hssw= 0;
Temp_diff_control_cdw= 0;

% Temp_diff_cplus20_n=0;
% Temp_diff_cplus10_n= 0;
% Temp_diff_cplus5_n= 0;
% Temp_diff_cminus5_n= 0;
% Temp_diff_cminus10_n= 0;
% Temp_diff_cminus20_n= 0;
% Temp_diff_cminus50_n= 0;
% Temp_diff_cminus40_n= 0 ; 
% Temp_diff_cplus15_n= 0;
                                                      
Temp_diff_control_whssw= 0;
Temp_diff_control_wcdw= 0;

% Temp_diff_wplus20_n= 0;
% Temp_diff_wplus10_n= 0;
% Temp_diff_wplus5_n= 0;
% Temp_diff_wminus5_n= 0;
% Temp_diff_wminus10_n= 0;
% Temp_diff_wplus15_n= 0;
% Temp_diff_wminus20_n= 0;
% Temp_diff_wminus40_n= 0;
% Temp_diff_wminus50_n= 0;
                                                      
Salt_diff_control_hssw= 0;
Salt_diff_control_cdw= 0;
% 
% Salt_diff_cplus20_n=0;
% Salt_diff_cplus10_n= 0;
% Salt_diff_cplus5_n= 0;
% Salt_diff_cminus5_n= 0;
% Salt_diff_cminus10_n= 0;
% Salt_diff_cminus20_n= 0;
% Salt_diff_cminus50_n= 0;
% Salt_diff_cminus40_n= 0 ; 
% Salt_diff_cplus15_n= 0;
                                                      
Salt_diff_wcontrol_whssw= 0;
Salt_diff_wcontrol_wcdw= 0;

% Salt_diff_wplus20_n= 0;
% Salt_diff_wplus10_n= 0;
% Salt_diff_wplus5_n= 0;
% Salt_diff_wminus5_n= 0;
% Salt_diff_wminus10_n= 0;
% Salt_diff_wplus15_n= 0;
% Salt_diff_wminus20_n= 0;
% Salt_diff_wminus40_n= 0;
% Salt_diff_wminus50_n= 0;
%%%500mdepth = zz(40)

 %%% Integrate melt over the FRIS
for i=1:Nx
        for j=1:Ny
            for k = 1:Nr
                if SHELFICEtopo(i,j)==0 && bathy(i,j)<-500 && bathy(i,j)>-2000
                    if XC(i,j) <-25 &&XC(i,j)>-50

                       if Temp_diff_control_cdw < Temp_diff_control(i,j,k)  
                           Temp_diff_control_cdw=Temp_diff_control(i,j,k);
                           Salt_diff_control_cdw=Salt_diff_control(i,j,k);
                       end
                       if Salt_diff_control_hssw < Salt_diff_control(i,j,k)  && Temp_diff_control(i,j,k)<-1
                           Temp_diff_control_hssw=Temp_diff_control(i,j,k);
                           Salt_diff_control_hssw=Salt_diff_control(i,j,k);
                       end                       
                       
%                        if Temp_diff_cplus20_n < Temp_diff_cplus20(i,j,k)  
%                         Temp_diff_cplus20_n=Temp_diff_cplus20(i,j,k);
%                        end
%                        if Salt_diff_cplus20_n < Salt_diff_cplus20(i,j,k)  && Temp_diff_cplus20(i,j,k)<-1
%                         Salt_diff_cplus20_n=Salt_diff_cplus20(i,j,k);
%                        end           
%                        
%                        if Temp_diff_cplus10_n < Temp_diff_cplus10(i,j,k)  
%                         Temp_diff_cplus10_n=Temp_diff_cplus10(i,j,k);
%                        end
%                        if Salt_diff_cplus10_n < Salt_diff_cplus10(i,j,k)  && Temp_diff_cplus10(i,j,k)<-1
%                         Salt_diff_cplus10_n=Salt_diff_cplus10(i,j,k);
%                        end                       
%                        
%                        if Temp_diff_cplus5_n < Temp_diff_cplus5(i,j,k)  
%                         Temp_diff_cplus5_n=Temp_diff_cplus5(i,j,k);
%                        end
%                        if Salt_diff_cplus5_n < Salt_diff_cplus5(i,j,k)  && Temp_diff_cplus5(i,j,k)<-1
%                         Salt_diff_cplus5_n=Salt_diff_cplus5(i,j,k);
%                        end                       
%                        
%                        if Temp_diff_cminus5_n < Temp_diff_cminus5(i,j,k)  
%                         Temp_diff_cminus5_n=Temp_diff_cminus5(i,j,k);
%                        end
%                        if Salt_diff_cminus5_n < Salt_diff_cminus5(i,j,k)  && Temp_diff_cminus5(i,j,k)<-1
%                         Salt_diff_cminus5_n=Salt_diff_cminus5(i,j,k);
%                        end                      
%                        
%                        if Temp_diff_cminus10_n < Temp_diff_cminus10(i,j,k)  
%                         Temp_diff_cminus10_n=Temp_diff_cminus10(i,j,k);
%                        end
%                        if Salt_diff_cminus10_n < Salt_diff_cminus10(i,j,k)  && Temp_diff_cminus10(i,j,k)<-1
%                         Salt_diff_cminus10_n=Salt_diff_cminus10(i,j,k);
%                        end                      
%                        
%                        if Temp_diff_cminus20_n < Temp_diff_cminus20(i,j,k)  
%                         Temp_diff_cminus20_n=Temp_diff_cminus20(i,j,k);
%                        end
%                        if Salt_diff_cminus20_n < Salt_diff_cminus20(i,j,k)  && Temp_diff_cminus20(i,j,k)<-1
%                         Salt_diff_cminus20_n=Salt_diff_cminus20(i,j,k);
%                        end                      
%                        
%                        if Temp_diff_cminus50_n < Temp_diff_cminus50(i,j,k)  
%                         Temp_diff_cminus50_n=Temp_diff_cminus50(i,j,k);
%                        end
%                        if Salt_diff_cminus50_n < Salt_diff_cminus50(i,j,k)  && Temp_diff_cminus50(i,j,k)<-1
%                         Salt_diff_cminus50_n=Salt_diff_cminus50(i,j,k);
%                        end
%                        
%                        if Temp_diff_cminus40_n < Temp_diff_cmerid40(i,j,k)  
%                         Temp_diff_cminus40_n=Temp_diff_cmerid40(i,j,k);
%                        end
%                        if Salt_diff_cminus40_n < Salt_diff_cmerid40(i,j,k)  && Temp_diff_cmerid40(i,j,k)<-1
%                         Salt_diff_cminus40_n=Salt_diff_cmerid40(i,j,k);
%                        end                       
%                        
%                        if Temp_diff_cplus15_n < Temp_diff_cplus15(i,j,k)  
%                         Temp_diff_cplus15_n=Temp_diff_cplus15(i,j,k);
%                        end
%                         
%                        if Salt_diff_cplus15_n < Salt_diff_cplus15(i,j,k)  && Temp_diff_cplus15(i,j,k)<-1
%                         Salt_diff_cplus15_n=Salt_diff_cplus15(i,j,k);
%                        end                      

                       if Temp_diff_control_wcdw < Temp_diff_w_control(i,j,k)  
                        Temp_diff_control_wcdw=Temp_diff_w_control(i,j,k);
                        Salt_diff_wcontrol_wcdw=Salt_diff_w_control(i,j,k);
                  
                       end
                       if Salt_diff_wcontrol_whssw < Salt_diff_w_control(i,j,k)  && Temp_diff_w_control(i,j,k)<-1
                        Salt_diff_wcontrol_whssw=Salt_diff_w_control(i,j,k);
                        Temp_diff_control_whssw=Temp_diff_w_control(i,j,k);

                       end                      
                       
%                        if Temp_diff_wplus20_n < Temp_diff_w_plus20(i,j,k)  
%                         Temp_diff_wplus20_n=Temp_diff_w_plus20(i,j,k);
%                        end
%                        if Salt_diff_wplus20_n < Salt_diff_w_plus20(i,j,k)  && Temp_diff_w_plus20(i,j,k)<-1
%                         Salt_diff_wplus20_n=Salt_diff_w_plus20(i,j,k);
%                        end     
%                        
%                        
%                        if Temp_diff_wplus10_n < Temp_diff_w_plus10(i,j,k)  
%                         Temp_diff_wplus10_n=Temp_diff_w_plus10(i,j,k);
%                        end
%                        if Salt_diff_wplus10_n < Salt_diff_w_plus10(i,j,k)  && Temp_diff_w_plus10(i,j,k)<-1
%                         Salt_diff_wplus10_n=Salt_diff_w_plus10(i,j,k);
%                        end                      
%                        
%                        if Temp_diff_wplus5_n < Temp_diff_w_plus5(i,j,k)  
%                         Temp_diff_wplus5_n=Temp_diff_w_plus5(i,j,k);
%                        end
%                        if Salt_diff_wplus5_n < Salt_diff_w_plus5(i,j,k)  && Temp_diff_w_plus5(i,j,k)<-1
%                         Salt_diff_wplus5_n=Salt_diff_w_plus5(i,j,k);
%                        end                       
%                        
%                        if Temp_diff_wminus5_n < Temp_diff_wminus5(i,j,k)  
%                         Temp_diff_wminus5_n=Temp_diff_wminus5(i,j,k);
%                        end
%                        if Salt_diff_wminus5_n < Salt_diff_wminus5(i,j,k)  && Temp_diff_wminus5(i,j,k)<-1
%                         Salt_diff_wminus5_n=Salt_diff_wminus5(i,j,k);
%                        end                       
%                        
%                        if Temp_diff_wminus10_n < Temp_diff_wminus10(i,j,k)  
%                         Temp_diff_wminus10_n=Temp_diff_wminus10(i,j,k);
%                        end
%                        if Salt_diff_wminus10_n < Salt_diff_wminus10(i,j,k)  && Temp_diff_wminus10(i,j,k)<-1
%                         Salt_diff_wminus10_n=Salt_diff_wminus10(i,j,k);
%                        end                       
%                        
%                        if Temp_diff_wminus20_n < Temp_diff_w_minus20(i,j,k)  
%                         Temp_diff_wminus20_n=Temp_diff_w_minus20(i,j,k);
%                        end
%                        if Salt_diff_wminus20_n < Salt_diff_w_minus20(i,j,k)  && Temp_diff_w_minus20(i,j,k)<-1
%                         Salt_diff_wminus20_n=Salt_diff_w_minus20(i,j,k);
%                        end                       
%                        
%                        if Temp_diff_wminus50_n < Temp_diff_wmerid_minus50(i,j,k)  
%                         Temp_diff_wminus50_n=Temp_diff_wmerid_minus50(i,j,k);
%                        end
%                        if Salt_diff_wminus50_n < Salt_diff_wmerid_minus50(i,j,k)  && Temp_diff_wmerid_minus50(i,j,k)<-1
%                         Salt_diff_wminus50_n=Salt_diff_wmerid_minus50(i,j,k);
%                        end                       
%                        
%                        if Temp_diff_wminus40_n < Temp_diff_wmerid_minus40(i,j,k)  
%                         Temp_diff_wminus40_n=Temp_diff_wmerid_minus40(i,j,k);
%                        end
%                        if Salt_diff_wminus40_n < Salt_diff_wmerid_minus40(i,j,k)  && Temp_diff_wmerid_minus40(i,j,k)<-1
%                         Salt_diff_wminus40_n=Salt_diff_wmerid_minus40(i,j,k);
%                        end                      
%                        
%                        if Temp_diff_wplus15_n < Temp_diff_w_plus15(i,j,k)  
%                         Temp_diff_wplus15_n=Temp_diff_w_plus15(i,j,k);
%                        end
%                        if Salt_diff_wplus15_n < Salt_diff_w_plus15(i,j,k)  && Temp_diff_w_plus15(i,j,k)<-1
%                         Salt_diff_wplus15_n=Salt_diff_w_plus15(i,j,k);
%                        end                       
                       
                    end
            end
            end
        end
end

% salt_controlwarmtotal= [Salt_diff_wminus50_n,Salt_diff_wminus40_n,Salt_diff_wminus20_n,Salt_diff_wminus10_n,Salt_diff_wminus5_n,Salt_diff_wcontrol_n, Salt_diff_wplus5_n,Salt_diff_wplus10_n,Salt_diff_wplus15_n,Salt_diff_wplus20_n];
% temp_controlwarmtotal= [Temp_diff_wminus50_n,Temp_diff_wminus40_n,Temp_diff_wminus20_n,Temp_diff_wminus10_n,Temp_diff_wminus5_n,Temp_diff_wcontrol_nhssw, Temp_diff_wplus5_n,Temp_diff_wplus10_n,Temp_diff_wplus15_n,Temp_diff_wplus20_n];
% 
%          
% salt_control=[Salt_diff_cminus50_n,Salt_diff_cminus40_n,Salt_diff_cminus20_n,Salt_diff_cminus10_n,Salt_diff_cminus5_n,Salt_diff_control_n,Salt_diff_cplus5_n,Salt_diff_cplus10_n,Salt_diff_cplus15_n,Salt_diff_cplus20_n];
% temp_control=[Temp_diff_cminus50_n,Temp_diff_cminus40_n,Temp_diff_cminus20_n,Temp_diff_cminus10_n,Temp_diff_cminus5_n,Temp_diff_control_hssw,Temp_diff_cplus5_n,Temp_diff_cplus10_n,Temp_diff_cplus15_n,Temp_diff_cplus20_n];


 save CDWHSSWvals_343445 Salt_diff_wcontrol_whssw Temp_diff_control_whssw Temp_diff_control_wcdw Salt_diff_wcontrol_wcdw Salt_diff_control_hssw Temp_diff_control_hssw Temp_diff_control_cdw Salt_diff_control_cdw



