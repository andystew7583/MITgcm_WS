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



exppath1w = '/data3/MITgcm_WS/experiments/a_34_20boundary';
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




% %%%%% calculate freezing temperature
% tmin = 10*86400*360;
% tmax = 17*86400*360;
% Salt =readIters(exppath5,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
% 
% %%% Find average salinity at 1000m
% Salt(Salt==0) = NaN;
% Salt_avg = squeeze(nanmean(Salt(:,:,53),1));
% Salt_avg = squeeze(nanmean(Salt_avg,2));
% Pressure = (-rho0*g*zz(53))/Pa1dbar;
% temp_f = (.0901 - .0575*Salt_avg) - (7.61e-4 * Pressure);


tmin = 9*86400*360;
tmax = 18*86400*360;

%%%% then we plot the melt against the difference in the cavity temperature
% %%%% compared to freezing
% 
% nDumps = nDump_f-nDump_s;
% 
melttot_control =readIters(exppath5,'SHIfwFlx',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);
Temp_diff_control = readIters(exppath5,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);



melttot_c_meridminus10 = readIters(exppath10,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_cminus10 = readIters(exppath10,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

melttot_c_meridplus10 = readIters(exppath2,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_cplus10 = readIters(exppath2,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);


melttot_c_meridminus20= readIters(exppath9,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_cminus20 = readIters(exppath9,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);



melttot_c_meridplus15 = readIters(exppath13,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_cplus15 =readIters(exppath13,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);



melttot_c_meridminus5 = readIters(exppath7,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_cminus5 = readIters(exppath7,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

tmin = 18*86400*360;
tmax = 27*86400*360;

melttot_c_meridminus50 = readIters(exppath1,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_cminus50 =readIters(exppath1,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

tmin = 9*86400*360;
tmax = 18*86400*360;

melttot_c_meridplus20 =  readIters(exppath3,'SHIfwFlx',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,1);
Temp_diff_cplus20 =  readIters(exppath3,'THETA',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);
 

melttot_c_meridplus5 =  readIters(exppath4,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_cplus5 = readIters(exppath4,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

tmin = 27*86400*360;
tmax = 36*86400*360;
melttot_c_merid40=  readIters(exppath12,'SHIfwFlx',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,1);
Temp_diff_cmerid40 = readIters(exppath12,'THETA',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);




%%%%% Warm
tmin = 9*86400*360;
tmax = 18*86400*360;
melttot_w_meridminus10 = readIters(exppath7w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_wminus10 = readIters(exppath7w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

melttot_w_meridplus10 = readIters(exppath5w,'SHIfwFlx',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,1);
Temp_diff_w_plus10 = readIters(exppath5w,'THETA',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);


melttot_w_meridminus20= readIters(exppath8w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_w_minus20 = readIters(exppath8w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);


melttot_w_meridminus5 =  readIters(exppath9w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_wminus5 = readIters(exppath9w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

tmin = 27*86400*360;
tmax = 36*86400*360;
melttot_w_meridplus20 = readIters(exppath2w,'SHIfwFlx',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,1);
Temp_diff_w_plus20 = readIters(exppath2w,'THETA',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);

tmin = 19*86400*360;
tmax = 27*86400*360;
Temp_diff_w_plus15 = readIters(exppath6w,'THETA',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,Nr);
melttot_w_meridplus15 = readIters(exppath6w,'SHIfwFlx',dumpIters4,deltaT4,tmin,tmax,Nx,Ny,1);

tmin = 9*86400*360;
tmax = 18*86400*360;
melttot_w_meridplus5 = readIters(exppath11w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_w_plus5 = readIters(exppath11w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);


melttot_w_meridminus50 = readIters(exppath4w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_wmerid_minus50 = readIters(exppath4w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);


melttot_w_meridminus40 = readIters(exppath13w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Temp_diff_wmerid_minus40 = readIters(exppath13w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

melttot_w_control = readIters(exppath1w,'SHIfwFlx',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);
Temp_diff_w_control = readIters(exppath1w,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);


melttot_control_n=0;
melttot_c_meridminus10_n=0;
melttot_c_meridminus20_n=0;
melttot_c_meridminus5_n=0;
melttot_c_meridminus15_n=0;
melttot_c_meridminus50_n=0;
melttot_c_meridminus40_n=0;
melttot_c_meridplus15_n=0;
melttot_c_meridplus5_n=0;
melttot_c_meridplus10_n=0;
melttot_c_meridplus20_n=0;

melttot_w_meridminus10_n=0;
melttot_w_meridminus5_n=0;
melttot_w_meridminus20_n=0;
melttot_w_meridplus20_n=0;
melttot_w_meridplus10_n=0;
melttot_w_meridplus5_n=0;
melttot_w_meridplus15_n=0;
melttot_w_control_n=0;
melttot_w_meridminus40_n=0;
melttot_w_meridminus50_n=0;



Temp_diff_control_n= 0;
Temp_diff_cplus20_n=0;
Temp_diff_cplus10_n= 0;
Temp_diff_cplus5_n= 0;
Temp_diff_cminus5_n= 0;
Temp_diff_cminus10_n= 0;
Temp_diff_cminus20_n= 0;
Temp_diff_cminus50_n= 0;
Temp_diff_cminus40_n= 0 ;                                                
Temp_diff_cplus15_n= 0;
                                                      
Temp_diff_w_control_n= 0;
Temp_diff_w_plus20_n= 0;
Temp_diff_w_plus10_n= 0;
Temp_diff_w_plus5_n= 0;
Temp_diff_w_minus5_n= 0;
Temp_diff_w_minus10_n= 0;
Temp_diff_w_plus15_n= 0;
Temp_diff_w_minus20_n= 0;
Temp_diff_w_minus40_n= 0;
Temp_diff_w_minus50_n= 0;
                                                      
Volume = 0;

 %%% Integrate melt over the FRIS
T2=0;
for i=1:Nx
        for j=1:Ny
                if XC(i,1)>-80 && XC(i,1)<-20
            
                  if YC(1,j) <-75
                      if SHELFICEtopo(i,j)<0

                        
%                         melttot_w_meridplus2_n = melttot_w_meridplus2_n+melttot_w_meridplus2(i,j)*RAC(i,j);
               


                       if melttot_control(i,j)<0     
                        T2=Temp_diff_control(i,j,kmax(i,j));
  
                        melttot_control_n= melttot_control_n+melttot_control(i,j)*RAC(i,j);
                        Temp_diff_control_n= (Temp_diff_control_n + (T2-thetafrz)*RAC(i,j)*melttot_control(i,j));
                       end
                       if melttot_c_meridplus20(i,j)<0
                        Temp_diff_cplus20_n= ( Temp_diff_cplus20_n+ (Temp_diff_cplus20(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_c_meridplus20(i,j));
                        melttot_c_meridplus20_n = melttot_c_meridplus20_n+melttot_c_meridplus20(i,j)*RAC(i,j);

                       end
                       if melttot_c_meridplus10(i,j)<0
                        Temp_diff_cplus10_n= (Temp_diff_cplus10_n+ (Temp_diff_cplus10(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_c_meridplus10(i,j));
                        melttot_c_meridplus10_n = melttot_c_meridplus10_n+melttot_c_meridplus10(i,j)*RAC(i,j);

                       end
                       if melttot_c_meridplus5(i,j)<0
                        Temp_diff_cplus5_n= ( Temp_diff_cplus5_n+ (Temp_diff_cplus5(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_c_meridplus5(i,j));
                        melttot_c_meridplus5_n = melttot_c_meridplus5_n+melttot_c_meridplus5(i,j)*RAC(i,j);

                       end
                       if melttot_c_meridminus5(i,j)<0
                        Temp_diff_cminus5_n= ( Temp_diff_cminus5_n+ (Temp_diff_cminus5(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_c_meridminus5(i,j));
                        melttot_c_meridminus5_n = melttot_c_meridminus5_n+melttot_c_meridminus5(i,j)*RAC(i,j);

                       end
 
                       if melttot_c_meridminus10(i,j)<0
                        melttot_c_meridminus10_n = melttot_c_meridminus10_n+melttot_c_meridminus10(i,j)*RAC(i,j);
                   
                        Temp_diff_cminus10_n= (Temp_diff_cminus10_n+ (Temp_diff_cminus10(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_c_meridminus10(i,j));
                       end
                       if melttot_c_meridminus20(i,j)<0
                        melttot_c_meridminus20_n = melttot_c_meridminus20_n + melttot_c_meridminus20(i,j)*RAC(i,j);
      
                        Temp_diff_cminus20_n= (Temp_diff_cminus20_n+(Temp_diff_cminus20(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_c_meridminus20(i,j));
                       end
                       if melttot_c_meridminus50(i,j)<0
                        melttot_c_meridminus50_n = melttot_c_meridminus50_n+ melttot_c_meridminus50(i,j)*RAC(i,j);
        
                        Temp_diff_cminus50_n= (Temp_diff_cminus50_n+ (Temp_diff_cminus50(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_c_meridminus50(i,j));
                       end
                       if melttot_c_merid40(i,j)<0
                        melttot_c_meridminus40_n = melttot_c_meridminus40_n+ melttot_c_merid40(i,j)*RAC(i,j);
              
                        Temp_diff_cminus40_n= (Temp_diff_cminus40_n+(Temp_diff_cmerid40(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_c_merid40(i,j)); 
                       end
                       if melttot_c_meridplus15(i,j)<0
                        Temp_diff_cplus15_n= ( Temp_diff_cplus15_n+ (Temp_diff_cplus15(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_c_meridplus15(i,j));
                        melttot_c_meridplus15_n = melttot_c_meridplus15_n+melttot_c_meridplus15(i,j)*RAC(i,j);
                       end
                       
                       
                       
                       
                       if melttot_w_control(i,j)<0
                        melttot_w_control_n = melttot_w_control_n + melttot_w_control(i,j)*RAC(i,j);
                    
                        Temp_diff_w_control_n= (Temp_diff_w_control_n+(Temp_diff_w_control(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_control(i,j));
                       end
                       if melttot_w_meridplus20(i,j)<0
                        Temp_diff_w_plus20_n= ( Temp_diff_w_plus20_n+ (Temp_diff_w_plus20(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_meridplus20(i,j));
                        melttot_w_meridplus20_n = melttot_w_meridplus20_n+melttot_w_meridplus20(i,j)*RAC(i,j);
                       end
                     
                       if melttot_w_meridplus10(i,j)<0
                        Temp_diff_w_plus10_n= ( Temp_diff_w_plus10_n+ (Temp_diff_w_plus10(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_meridplus10(i,j));
                        melttot_w_meridplus10_n = melttot_w_meridplus10_n+melttot_w_meridplus10(i,j)*RAC(i,j);
                      
                       end
                       
                       if melttot_w_meridplus5(i,j)<0
                        Temp_diff_w_plus5_n= ( Temp_diff_w_plus5_n+ (Temp_diff_w_plus5(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_meridplus5(i,j));
                        melttot_w_meridplus5_n = melttot_w_meridplus5_n+melttot_w_meridplus5(i,j)*RAC(i,j);
                       end
                       if melttot_w_meridminus5(i,j)<0
                        Temp_diff_w_minus5_n= ( Temp_diff_w_minus5_n+ (Temp_diff_wminus5(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_meridminus5(i,j));
                        melttot_w_meridminus5_n = melttot_w_meridminus5_n+melttot_w_meridminus5(i,j)*RAC(i,j);
                       end
                       if melttot_w_meridminus10(i,j)<0
                        Temp_diff_w_minus10_n= ( Temp_diff_w_minus10_n+ (Temp_diff_wminus10(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_meridminus10(i,j));
                        melttot_w_meridminus10_n = melttot_w_meridminus10_n+melttot_w_meridminus10(i,j)*RAC(i,j);
                       end
                       if melttot_w_meridplus15(i,j)<0
                        Temp_diff_w_plus15_n= ( Temp_diff_w_plus15_n+ (Temp_diff_w_plus15(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_meridplus15(i,j));
                        melttot_w_meridplus15_n = melttot_w_meridplus15_n+melttot_w_meridplus15(i,j)*RAC(i,j);
                       end
                       
                       if melttot_w_meridminus20(i,j)<0
                        Temp_diff_w_minus20_n= ( Temp_diff_w_minus20_n+ (Temp_diff_w_minus20(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_meridminus20(i,j));
                        melttot_w_meridminus20_n = melttot_w_meridminus20_n+melttot_w_meridminus20(i,j)*RAC(i,j);
                             
                       end
                       if melttot_w_meridminus40(i,j)<0
                        Temp_diff_w_minus40_n= ( Temp_diff_w_minus40_n+ (Temp_diff_wmerid_minus40(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_meridminus40(i,j));
                        melttot_w_meridminus40_n = melttot_w_meridminus40_n+melttot_w_meridminus40(i,j)*RAC(i,j);
                       end
                       if melttot_w_meridminus50(i,j)<0
                        Temp_diff_w_minus50_n= ( Temp_diff_w_minus50_n+ (Temp_diff_wmerid_minus50(i,j,kmax(i,j))-thetafrz)*RAC(i,j)*melttot_w_meridminus50(i,j));
                        melttot_w_meridminus50_n = melttot_w_meridminus50_n+melttot_w_meridminus50(i,j)*RAC(i,j);
                       end
                      end
                  end
                end
        end
end


                

  

Temp_diff_control_n= Temp_diff_control_n./melttot_control_n;
Temp_diff_cplus20_n=Temp_diff_cplus20_n./melttot_c_meridplus20_n;
Temp_diff_cplus10_n= Temp_diff_cplus10_n./melttot_c_meridplus10_n;
Temp_diff_cplus5_n= Temp_diff_cplus5_n./melttot_c_meridplus5_n;
Temp_diff_cminus5_n=Temp_diff_cminus5_n./melttot_c_meridminus5_n;
Temp_diff_cminus10_n=Temp_diff_cminus10_n./melttot_c_meridminus10_n;
Temp_diff_cminus20_n= Temp_diff_cminus20_n./melttot_c_meridminus20_n;
Temp_diff_cminus50_n= Temp_diff_cminus50_n./melttot_c_meridminus50_n;
Temp_diff_cminus40_n= Temp_diff_cminus40_n./melttot_c_meridminus40_n;                                                
Temp_diff_cplus15_n= Temp_diff_cplus15_n./melttot_c_meridplus15_n;
                                                      
Temp_diff_w_control_n= Temp_diff_w_control_n./melttot_w_control_n;
Temp_diff_w_plus20_n=Temp_diff_w_plus20_n./melttot_w_meridplus20_n;
Temp_diff_w_plus10_n= Temp_diff_w_plus10_n./melttot_w_meridplus10_n;
Temp_diff_w_plus5_n= Temp_diff_w_plus5_n./melttot_w_meridplus5_n;
Temp_diff_w_minus5_n=Temp_diff_w_minus5_n./melttot_w_meridminus5_n;
Temp_diff_w_minus10_n=Temp_diff_w_minus10_n./melttot_w_meridminus10_n;
Temp_diff_w_minus20_n= Temp_diff_w_minus20_n./melttot_w_meridminus20_n;
Temp_diff_w_minus50_n= Temp_diff_w_minus50_n./melttot_w_meridminus50_n;
Temp_diff_w_minus40_n= Temp_diff_w_minus40_n./melttot_w_meridminus40_n;                                                
Temp_diff_w_plus15_n= Temp_diff_w_plus15_n./melttot_w_meridplus15_n;
                              


%%%convert to Salt Flux

melttot_w_control_n = (melttot_w_control_n)*34; 
melttot_w_meridplus5_n = (melttot_w_meridplus5_n)*34;  
melttot_w_meridplus10_n = (melttot_w_meridplus10_n)*34;  
melttot_w_meridplus15_n = (melttot_w_meridplus15_n)*34; 
melttot_w_meridplus20_n = (melttot_w_meridplus20_n)*34;  
melttot_w_meridminus5_n = (melttot_w_meridminus5_n)*34;   
melttot_w_meridminus10_n = (melttot_w_meridminus10_n)*34;    
melttot_w_meridminus20_n = (melttot_w_meridminus20_n)*34;    
melttot_w_meridminus40_n = (melttot_w_meridminus40_n)*34;    
melttot_w_meridminus50_n = (melttot_w_meridminus50_n)*34;   

melttot_control_n = (melttot_control_n)*34;   
melttot_c_meridplus5_n = (melttot_c_meridplus5_n)*34; 
melttot_c_meridplus10_n = (melttot_c_meridplus10_n)*34;  
melttot_c_meridplus15_n = (melttot_c_meridplus15_n)*34;  
melttot_c_meridplus20_n = (melttot_c_meridplus20_n)*34;  
melttot_c_meridminus5_n = (melttot_c_meridminus5_n)*34;   
melttot_c_meridminus10_n = (melttot_c_meridminus10_n)*34;  
melttot_c_meridminus20_n = (melttot_c_meridminus20_n)*34;  
melttot_c_meridminus50_n = (melttot_c_meridminus50_n)*34;    
melttot_c_meridminus40_n = (melttot_c_meridminus40_n)*34;    


%%%%make plot
totc= [Temp_diff_cminus20_n,Temp_diff_cminus10_n,Temp_diff_cminus5_n,Temp_diff_control_n,Temp_diff_cplus5_n,Temp_diff_cplus10_n,Temp_diff_cplus15_n,Temp_diff_w_plus15_n,Temp_diff_cplus20_n,Temp_diff_w_plus20_n];
                                                      
totmc=[melttot_c_meridminus20_n,melttot_c_meridminus10_n, melttot_c_meridminus5_n,melttot_control_n,melttot_c_meridplus5_n, melttot_c_meridplus10_n,melttot_c_meridplus15_n,melttot_c_meridplus15_n,melttot_c_meridplus20_n,melttot_w_meridplus20_n];

totw=[Temp_diff_w_minus50_n,Temp_diff_w_minus40_n,Temp_diff_cminus50_n,Temp_diff_cminus40_n,Temp_diff_w_minus20_n,Temp_diff_w_minus10_n,Temp_diff_w_minus5_n,Temp_diff_w_control_n,Temp_diff_w_plus5_n, Temp_diff_w_plus10_n];

wm=[melttot_w_meridminus50_n,melttot_w_meridminus40_n,melttot_c_meridminus50_n,melttot_c_meridminus40_n,melttot_w_meridminus20_n,melttot_w_meridminus10_n,melttot_w_meridminus5_n,melttot_w_control_n,melttot_w_meridplus5_n,melttot_w_meridplus10_n];
totc=sort(totc,'ascend');
totmc=sort(totmc,'descend');
totw=sort(totw,'ascend');
wm=sort(wm,'descend');

tot_temp = horzcat(totc,totw);
tot_temp = sort(tot_temp,'ascend');

tot_melt = horzcat(totmc,wm);
tot_melt = sort(tot_melt,'descend');

figure(1);
clf;
hold on


plot(Temp_diff_control_n,melttot_control_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_control_n+.001,melttot_control_n-.005,'C')

hold on
plot(Temp_diff_cplus5_n,melttot_c_meridplus5_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cplus5_n+.001,melttot_c_meridplus5_n-.005,'+5%')

hold on
plot(Temp_diff_cplus10_n,melttot_c_meridplus10_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cplus5_n+.001,melttot_c_meridplus5_n-.005,'+10%')

hold on
plot(Temp_diff_cplus15_n,melttot_c_meridplus15_n,'-bo','markerfacecolor','blue','markersize',10);

hold on
plot(Temp_diff_cplus10_n,melttot_c_meridplus10_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cplus10_n+.001,melttot_c_meridplus5_n-.005,'+10%')

hold on
plot(Temp_diff_cplus20_n,melttot_c_meridplus20_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cplus20_n+.001,melttot_c_meridplus20_n-.005,'+20%')
hold on

plot(Temp_diff_cminus50_n,melttot_c_meridminus50_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_cplus20_n+.001,melttot_c_meridplus20_n-.005,'+20%')

hold on
% plot(Temp_diff_cminus2_n,melttot_c_meridminus2_n,'-bo','markerfacecolor','blue');
% hold on
plot(Temp_diff_cminus5_n,melttot_c_meridminus5_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cminus5_n+.001,melttot_c_meridminus5_n-.005,'-5%')

hold on
plot(Temp_diff_cminus10_n,melttot_c_meridminus10_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cminus10_n+.001,melttot_c_meridminus10_n-.005,'-10%')

hold on
plot(Temp_diff_cminus20_n,melttot_c_meridminus20_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cminus20_n+.001,melttot_c_meridminus20_n-.005,'-5%')'markersize',10

hold on
plot(Temp_diff_cminus40_n ,melttot_c_meridminus40_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_cminus15_n+.001,melttot_c_meridminus15_n-.005,'-15%')

hold on




plot(Temp_diff_w_control_n,melttot_w_control_n,'-ro','markerfacecolor','red','markersize',10);
hold on
% text(Temp_diff_w_control_n+.001,melttot_w_control_n-.005,'Control')
% 
plot(Temp_diff_w_plus5_n,melttot_w_meridplus5_n,'-ro','markerfacecolor','red','markersize',10);
% % text(Temp_diff_w_plus5_n,melttot_w_meridplus5_n,'+5%')

hold on
plot(Temp_diff_w_plus10_n,melttot_w_meridplus10_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_w_plus10_n+.001,melttot_w_meridplus10_n+.001,'+10%')

hold on
% plot(Temp_diff_w_plus2_n,melttot_w_meridplus2_n,'-ro','markerfacecolor','red');
% text(Temp_diff_w_plus2_n+.001,melttot_w_meridplus2_n+.001,'+2%')
% hold on

plot(Temp_diff_w_plus20_n,melttot_w_meridplus20_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_w_plus20_n+.001,melttot_w_meridplus20_n+.001,'+20%')

hold on
plot(Temp_diff_w_minus5_n,melttot_w_meridminus5_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_w_minus5_n+.001,melttot_w_meridminus5_n+.001,'-5%')

hold on
plot(Temp_diff_w_minus10_n,melttot_w_meridminus10_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_w_minus10_n+.0001,melttot_w_meridminus10_n+.0001,'-10%')

hold on
plot(Temp_diff_w_minus20_n,melttot_w_meridminus20_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_w_minus20_n+.0001,melttot_w_meridminus20_n+.0001,'-20%')


hold on
plot(Temp_diff_w_plus15_n,melttot_w_meridplus15_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_w_plus15_n+.001,melttot_w_meridplus15_n+.001,'+15%')

hold on
plot(Temp_diff_w_minus40_n,melttot_w_meridminus40_n,'-ro','markerfacecolor','red','markersize',10);
% xlim([.15 .25]);

hold on
plot(Temp_diff_w_minus50_n,melttot_w_meridminus50_n,'-ro','markerfacecolor','red','markersize',10);






xlim([0 1.5])
r = fitlm(totc,totmc);
p = polyfit((totc),(totmc),1);
z = polyval(p,totc,1);
hold on
plot(totc,z,'--b','linewidth',3);

f = fitlm(totw,wm);
o = polyfit((totw),(wm),1);
h = polyval(o,totw);
hold on
plot(totw,h,'--r','linewidth',3);


totfit = fitlm(tot_temp,tot_melt);

% y = fitlm(tc,wc);
% o = polyfit((tc),(wc),1);
% h2 = polyval(o,tc);
% hold on
% plot(tc,h2,'--b','linewidth',2);


text(.2,-12e8,'R$^2$=.92','fontsize',19,'color','b','interpreter','latex')
text(.2,-11e8,'R$^2$=.98','fontsize',19,'color','r','interpreter','latex')
% text(1.1,.5e10,'R$^2$=.99','fontsize',19,'color','b','interpreter','latex')

%%%line of best fit
title('FRIS Freezing Bottom Temperature Anomaly($^\circ$C) vs Salt Flux from Freshwater $(g/s)$','FontSize',15,'interpreter','latex');
ylabel('Integrated Salt Flux from Freshwater $(g/s)$','FontSize',16,'interpreter','latex')
xlabel('Freezing Bottom Temperature Anomaly ($^o$C)','FontSize',16,'interpreter','latex');
ylim([-1.5e9 -1e8]);
set(gca,'fontsize',14);

%%int = 6E8;
