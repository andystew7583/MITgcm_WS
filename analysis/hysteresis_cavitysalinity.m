%%%% plot hysteresis

%%% Read experiment data
loadexp;


%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
deltaT = 440;
% nIter0 = 0;
nDumps = round(nTimeSteps*(deltaT/dumpFreq));
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);
 
deltaT_4 = 200;
nIter0_4 = 1;
nDumps_4 = round(nTimeSteps*(deltaT_4/dumpFreq));
dumpIters_4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters_4 = dumpIters_4(dumpIters_4 >= nIter0_4);



deltaT_5 = 200;
nIter0_5 = 2877120;
nDumps_5 = round(nTimeSteps*10*(deltaT_5/dumpFreq));
dumpIters_5 = round((1:nDumps_5)*dumpFreq/deltaT_5); 
dumpIters_5 = dumpIters_5(dumpIters_5 >= nIter0_5);

deltaT2 = 400;
nDumps2 = round(nTimeSteps*10*(deltaT2/dumpFreq));
dumpIters_2 = round((1:nDumps2)*dumpFreq/deltaT2); 
dumpIters2 = dumpIters_2(dumpIters_2 >= nIter0);
nDumps2 = length(dumpIters2);


deltaT_3 = 300;
nIter0_3 = 1;
nDumps_3 = round(nTimeSteps*10*(deltaT_3/dumpFreq));
dumpIters_3 = round((1:nDumps_3)*dumpFreq/deltaT_3); 
dumpIters_3 = dumpIters_3(dumpIters_3 >= nIter0_3);


%%%% Load Salinity and Temperature
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


%%%%Warm Experiments first

tmin = 9*86400*360;
tmax = 12*86400*360;
Cavity_temp_wcontrol = readIters(exppath1w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_wcontrol = readIters(exppath1w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

tmin = 13*86400*360;
tmax = 18*86400*360;
Cavity_temp_wcontrol2 = readIters(exppath1w,'SALT',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Melt_wcontrol2 = readIters(exppath1w,'SHIfwFlx',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,1);

tc=NaN(Nx,Ny,Nr);
mc = NaN(Nx,Nx,Nr);
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nr
            mc(i,j,k) = (Melt_wcontrol(i,j)+Melt_wcontrol2(i,j))/2;
            tc(i,j,k) =(Cavity_temp_wcontrol(i,j,k) + Cavity_temp_wcontrol2(i,j,k))/2;
        end
    end
end
Cavity_temp_wcontrol=tc;
Melt_wcontrol=mc;



tmin = 9*86400*360;
tmax = 18*86400*360;
Cavity_temp_wminus50 = readIters(exppath4w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_wminus50= readIters(exppath4w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_wminus20 = readIters(exppath8w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_wminus20= readIters(exppath8w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_wminus30 = readIters(exppath12w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_wminus30= readIters(exppath12w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_wminus5 = readIters(exppath9w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_wminus5= readIters(exppath9w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_wminus10 = readIters(exppath7w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_wminus10= readIters(exppath7w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_wplus5 = readIters(exppath11w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_wplus5= readIters(exppath11w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_wplus10 = readIters(exppath5w,'SALT',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Melt_wplus10= readIters(exppath5w,'SHIfwFlx',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,1);

tmin = 9*86400*360;
tmax = 18*86400*360;
Cavity_temp_wminus40 = readIters(exppath13w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_wminus40 = readIters(exppath13w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
tmin = 26*86400*360;
tmax = 35*86400*360;
Cavity_temp_wplus20 = readIters(exppath2w,'SALT',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
Melt_wplus20 = readIters(exppath2w,'SHIfwFlx',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,1);
tmin = 19*86400*360;
tmax = 28*86400*360;
Cavity_temp_wplus15 = readIters(exppath6w,'SALT',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Melt_wplus15= readIters(exppath6w,'SHIfwFlx',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,1);

Cavity_temp_wcontroln2=0;           
Cavity_temp_wminus10n2=0;
Cavity_temp_wminus20n2=0;
Cavity_temp_wminus5n2=0;
Cavity_temp_wminus50n2=0;
Cavity_temp_wminus40n2=0;
Cavity_temp_wminus30n2=0;
Cavity_temp_wplus20n2=0;
Cavity_temp_wplus10n2=0;
Cavity_temp_wplus5n2=0;
Cavity_temp_wplus15n2=0;
Meltw = 0;
Meltwplus5=0;
Meltwplus10=0;
Meltwplus20=0;
Meltwplus15=0;
Meltwminus5=0;
Meltwminus20=0;
Meltwminus10=0;
Meltwminus30=0;
Meltwminus50=0;
Meltwminus40=0;
for i = 1:Nx
    for j = 1:Ny
        if hFacC(i,j,kmax(i,j))>0
                 if XC(i,j)>-80 && XC(i,j)<-30
                         if YC(i,j)<-75
                            if  SHELFICEtopo(i,j)<0  
                             if Melt_wcontrol(i,j)<0    
                                Cavity_temp_wcontroln2 = Cavity_temp_wcontroln2+Cavity_temp_wcontrol(i,j,kmax(i,j))*RAC(i,j)*Melt_wcontrol(i,j); 
                                Meltw = Meltw + RAC(i,j)*Melt_wcontrol(i,j);
                             end
                             if Melt_wminus10(i,j)<0   
                                Meltwminus10=Meltwminus10+RAC(i,j)*Melt_wminus10(i,j);
                             
                                Cavity_temp_wminus10n2 = Cavity_temp_wminus10n2+Cavity_temp_wminus10(i,j,kmax(i,j))*RAC(i,j)*Melt_wminus10(i,j);
                             end
                             if Melt_wminus20(i,j)<0   
                              
                               Meltwminus20=Meltwminus20+RAC(i,j)*Melt_wminus20(i,j);
                              Cavity_temp_wminus20n2 = Cavity_temp_wminus20n2+Cavity_temp_wminus20(i,j,kmax(i,j))*RAC(i,j)*Melt_wminus20(i,j);
                             end
                             if Melt_wminus5(i,j)<0   
                              
                               Meltwminus5=Meltwminus5+RAC(i,j)*Melt_wminus5(i,j);
                               Cavity_temp_wminus5n2 = Cavity_temp_wminus5n2+Cavity_temp_wminus5(i,j,kmax(i,j))*RAC(i,j)*Melt_wminus5(i,j);
                             end
                             if Melt_wminus50(i,j)<0   
                                Meltwminus50=Meltwminus50+RAC(i,j)*Melt_wminus50(i,j);
                         
                                Cavity_temp_wminus50n2 = Cavity_temp_wminus50n2+Cavity_temp_wminus50(i,j,kmax(i,j))*RAC(i,j)*Melt_wminus50(i,j);
                             end
                             if Melt_wminus40(i,j)<0   
                                 Meltwminus40=Meltwminus40+RAC(i,j)*Melt_wminus40(i,j);
                           
                                Cavity_temp_wminus40n2 = Cavity_temp_wminus40n2+Cavity_temp_wminus40(i,j,kmax(i,j))*RAC(i,j)*Melt_wminus40(i,j);
                             end
                             if Melt_wminus30(i,j)<0   
                                Meltwminus30=Meltwminus30+RAC(i,j)*Melt_wminus30(i,j);
                           
                                Cavity_temp_wminus30n2 = Cavity_temp_wminus30n2 +Cavity_temp_wminus30(i,j,kmax(i,j))*RAC(i,j)*Melt_wminus30(i,j);
                             end
                             if Melt_wplus20(i,j)<0   
                                 Meltwplus20=Meltwplus20+RAC(i,j)*Melt_wplus20(i,j);
                             
                                Cavity_temp_wplus20n2 = Cavity_temp_wplus20n2+Cavity_temp_wplus20(i,j,kmax(i,j))*RAC(i,j)*Melt_wplus20(i,j);
                             end
                             if Melt_wplus10(i,j)<0   
                                Meltwplus10=Meltwplus10+RAC(i,j)*Melt_wplus10(i,j);
                              
                                Cavity_temp_wplus10n2 = Cavity_temp_wplus10n2+ Cavity_temp_wplus10(i,j,kmax(i,j))*RAC(i,j)*Melt_wplus10(i,j);
                             end
                             if Melt_wplus5(i,j)<0   
                                Meltwplus5=Meltwplus5+RAC(i,j)*Melt_wplus5(i,j);
                               
                                Cavity_temp_wplus5n2 = Cavity_temp_wplus5n2+Cavity_temp_wplus5(i,j,kmax(i,j))*RAC(i,j)*Melt_wplus5(i,j);
                             end
                             if Melt_wplus15(i,j)<0   
                               
                                Cavity_temp_wplus15n2 = Cavity_temp_wplus15n2+Cavity_temp_wplus15(i,j,kmax(i,j))*RAC(i,j)*Melt_wplus15(i,j);
                                                         
                                Meltwplus15=Meltwplus15+RAC(i,j)*Melt_wplus15(i,j);
                             end
                            end
                           
                         end
                 end
        end
    end
end




Cavity_salt_wminus10n=(Cavity_temp_wminus10n2/Meltwminus10);
Cavity_salt_wminus5n=(Cavity_temp_wminus5n2/Meltwminus5);
Cavity_salt_wminus20n=(Cavity_temp_wminus20n2/Meltwminus20);
Cavity_salt_wminus30n=(Cavity_temp_wminus30n2/Meltwminus30);
Cavity_salt_wminus40n=(Cavity_temp_wminus40n2/Meltwminus40);
Cavity_salt_wminus50n=(Cavity_temp_wminus50n2/Meltwminus50);
Cavity_salt_wcontroln=Cavity_temp_wcontroln2/Meltw;
Cavity_salt_wplus10n=Cavity_temp_wplus10n2/Meltwplus10;
Cavity_salt_wplus5n=Cavity_temp_wplus5n2/Meltwplus5;
Cavity_salt_wplus20n=Cavity_temp_wplus20n2/Meltwplus20;
Cavity_salt_wplus15n=Cavity_temp_wplus15n2/Meltwplus15;



% Cavity_temp_wminus10n=reshape(Cavity_temp_wminus10n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wminus5n=reshape(Cavity_temp_wminus5n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wminus20n=reshape(Cavity_temp_wminus20n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wminus40n=reshape(Cavity_temp_wminus40n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wminus50n=reshape(Cavity_temp_wminus50n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wcontroln=reshape(Cavity_temp_wcontroln2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wplus10n=reshape(Cavity_temp_wplus10n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wplus5n=reshape(Cavity_temp_wplus5n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wplus20n=reshape(Cavity_temp_wplus20n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wplus15n=reshape(Cavity_temp_wplus15n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));
% Cavity_temp_wminus30n=reshape(Cavity_temp_wminus30n2,1,size(Cavity_temp_wcontroln2,1)*size(Cavity_temp_wcontroln2,2));


 
 
 

% 
% Cavity_temp_wcontroln =  Cavity_temp_wcontroln(~Cavity_temp_wcontroln==0);
% temp_wcontrol_q3 = median(Cavity_temp_wcontroln(find(Cavity_temp_wcontroln>median(Cavity_temp_wcontroln))));
% % temp_wcontrol_q3 = min(min(min(Cavity_temp_wcontroln)));
% 
% 
% 
% Cavity_temp_wminus20n =  Cavity_temp_wminus20n(~Cavity_temp_wminus20n==0);
% s_wminus20_q3 = median(Cavity_temp_wminus20n(find(Cavity_temp_wminus20n>median(Cavity_temp_wminus20n))));
% % s_wminus20_q3 = min(min(min(Cavity_temp_wminus20n)));
% 
% 
% 
% Cavity_temp_wminus10n =  Cavity_temp_wminus10n(~Cavity_temp_wminus10n==0);
% s_wminus10_q3 = median(Cavity_temp_wminus10n(find(Cavity_temp_wminus10n>median(Cavity_temp_wminus10n))));
% % s_wminus10_q3 = min(min(min(Cavity_temp_wminus10n)));
% 
% Cavity_temp_wplus15n =  Cavity_temp_wplus15n(~Cavity_temp_wplus15n==0);
% s_wplus15_q3 = median(Cavity_temp_wplus15n(find(Cavity_temp_wplus15n>median(Cavity_temp_wplus15n))));
% % s_wplus15_q3 = min(min(min(Cavity_temp_wplus15n)));
% 
% 
% 
% Cavity_temp_wminus5n =  Cavity_temp_wminus5n(~Cavity_temp_wminus5n==0);
% s_wminus5_q3 = median(Cavity_temp_wminus5n(find(Cavity_temp_wminus5n>median(Cavity_temp_wminus5n))));
% % s_wminus5_q3 = min(min(min(Cavity_temp_wminus5n)));
% 
% 
% Cavity_temp_wminus40n =  Cavity_temp_wminus40n(~Cavity_temp_wminus40n==0);
% s_wminus40_q3 = median(Cavity_temp_wminus40n(find(Cavity_temp_wminus40n>median(Cavity_temp_wminus40n))));
% % s_wminus40_q3 = min(min(min(Cavity_temp_wminus40n)));
% 
% Cavity_temp_wminus50n =  Cavity_temp_wminus50n(~Cavity_temp_wminus50n==0);
% s_wminus50_q3 = median(Cavity_temp_wminus50n(find(Cavity_temp_wminus50n>median(Cavity_temp_wminus50n))));
% % s_wminus50_q3 = min(min(min(Cavity_temp_wminus50n)));
% 
% Cavity_temp_wplus10n =  Cavity_temp_wplus10n(~Cavity_temp_wplus10n==0);
% s_wplus10_q3 = median(Cavity_temp_wplus10n(find(Cavity_temp_wplus10n>median(Cavity_temp_wplus10n))));
% % s_wplus10_q3 = min(min(min(Cavity_temp_wplus10n)));
% 
% Cavity_temp_wplus20n =  Cavity_temp_wplus20n(~Cavity_temp_wplus20n==0);
% s_wplus20_q3 = median(Cavity_temp_wplus20n(find(Cavity_temp_wplus20n>median(Cavity_temp_wplus20n))));
% % s_wplus20_q3 = min(min(min((Cavity_temp_wplus20n))));
% 
% Cavity_temp_wminus30n =  Cavity_temp_wminus30n(~Cavity_temp_wminus30n==0);
% s_wminus30_q3 = median(Cavity_temp_wminus30n(find(Cavity_temp_wminus30n>median(Cavity_temp_wminus30n))));
% % s_wminus30_q3 = min(min(min(Cavity_temp_wminus30n)));
% 
% Cavity_temp_wplus5n =  Cavity_temp_wplus5n(~Cavity_temp_wplus5n==0);
% s_wplus5_q3 = median(Cavity_temp_wplus5n(find(Cavity_temp_wplus5n>median(Cavity_temp_wplus5n))));
% s_wplus5_q3 = min(min(min(Cavity_temp_wplus5n)));




% err_exp1w = std(Cavity_temp_wcontroln2(~Cavity_temp_wcontroln2==0));
% err_exp7w = std(Cavity_temp_wminus10n2(~Cavity_temp_wminus10n2==0));
% err_exp9w = std(Cavity_temp_wminus20n2(~Cavity_temp_wminus20n2==0));
% err_exp6w = std(Cavity_temp_wminus5n2(~Cavity_temp_wminus5n2==0));
% err_exp8w = std(Cavity_temp_wminus50n2(~Cavity_temp_wminus50n2==0));
% err_exp2w = std(Cavity_temp_wminus40n2(~Cavity_temp_wminus40n2==0));
% err_exp5w = std(Cavity_temp_wminus30n2(~Cavity_temp_wminus30n2==0));
% err_exp4w = std(Cavity_temp_wplus20n2(~Cavity_temp_wplus20n2==0));
% err_exp13w = std(Cavity_temp_wplus10n2(~Cavity_temp_wplus10n2==0));
% err_exp12w = std(Cavity_temp_wplus5n2(~Cavity_temp_wplus5n2==0));
% err_exp11w = std(Cavity_temp_wplus15n2(~Cavity_temp_wminus10n2==0));



%%%%cold


%%%%%%%Control experiments
tmin = 10*86400*360;
tmax = 19*86400*365;
Cavity_temp_cplus20 = readIters(exppath3,'SALT',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
Melt_cplus20=readIters(exppath3,'SHIfwFlx',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,1);



tmin = 11*86400*360;
tmax = 20*86400*365;
Cavity_temp_control = readIters(exppath5,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
Melt_control= readIters(exppath5,'SHIfwFlx',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);
tmin = 9*86400*360;
tmax = 18*86400*365;
Cavity_temp_cminus20 = readIters(exppath9,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_cminus20=readIters(exppath9,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_cminus5 = readIters(exppath7,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_cminus5= readIters(exppath7,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_cminus10 = readIters(exppath10,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_cminus10= readIters(exppath10,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_cplus5 = readIters(exppath4,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_cplus5= readIters(exppath4,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_cplus10 = readIters(exppath2,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_cplus10= readIters(exppath2,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Cavity_temp_cplus15 = readIters(exppath13,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_cplus15= readIters(exppath13,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

tmin = 18*86400*360;
tmax = 27*86400*365;
Cavity_temp_cminus50 = readIters(exppath1,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_cminus50= readIters(exppath1,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
tmin = 27*86400*360;
tmax = 36*86400*365;
Cavity_temp_cminus40 = readIters(exppath12,'SALT',dumpIters_2,deltaT2,tmin,tmax,Nx,Ny,Nr);
Melt_cminus40= readIters(exppath12,'SHIfwFlx',dumpIters_2,deltaT2,tmin,tmax,Nx,Ny,1);

tmin = 7*86400*360;
tmax = 16*86400*365;
Cavity_temp_cminus30 = readIters(exppath11,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
Melt_cminus30= readIters(exppath11,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

% Cavity_temp_controln2=zeros(Nx,Ny);           
% Cavity_temp_cminus10n2=zeros(Nx,Ny);
% Cavity_temp_cminus20n2=zeros(Nx,Ny);
% Cavity_temp_cminus5n2=zeros(Nx,Ny);
% Cavity_temp_cminus50n2=zeros(Nx,Ny);
% Cavity_temp_cminus40n2=zeros(Nx,Ny);
% Cavity_temp_cminus30n2=zeros(Nx,Ny);
% Cavity_temp_cplus20n2=zeros(Nx,Ny);
% Cavity_temp_cplus10n2=zeros(Nx,Ny);
% Cavity_temp_cplus5n2=zeros(Nx,Ny);
% Cavity_temp_cplus15n2=zeros(Nx,Ny);
% 
% for i = 1:Nx
%     for j = 1:Ny
%                     if XC(i,j)>-80 && XC(i,j)<-60
% 
%                          if YC(i,j)<-75
%                              if  SHELFICEtopo(i,j)<0                    
%                                 Cavity_temp_controln2(i,j) = Cavity_temp_control(i,j,kmax(i,j));             
%                                 Cavity_temp_cminus10n2(i,j) = Cavity_temp_cminus10(i,j,kmax(i,j));
%                                 Cavity_temp_cminus20n2(i,j) = Cavity_temp_cminus20(i,j,kmax(i,j));
%                                 Cavity_temp_cminus5n2(i,j) = Cavity_temp_cminus5(i,j,kmax(i,j));
%                                 Cavity_temp_cminus50n2(i,j) = Cavity_temp_cminus50(i,j,kmax(i,j));
%                                 Cavity_temp_cminus40n2(i,j) = Cavity_temp_cminus40(i,j,kmax(i,j));
%                                 Cavity_temp_cminus30n2(i,j) = Cavity_temp_cminus30(i,j,kmax(i,j));
%                                 Cavity_temp_cplus20n2(i,j) = Cavity_temp_cplus20(i,j,kmax(i,j));
%                                 Cavity_temp_cplus10n2(i,j) = Cavity_temp_cplus10(i,j,kmax(i,j));
%                                 Cavity_temp_cplus5n2(i,j) = Cavity_temp_cplus5(i,j,kmax(i,j));
%                                 Cavity_temp_cplus15n2(i,j) = Cavity_temp_cplus15(i,j,kmax(i,j));
%              
%                   
%                              end
%                          end
%                     end
%                  
% 
%         
%     end
% end
% 
% 
% Cavity_temp_cminus10n=reshape(Cavity_temp_cminus10n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_cminus5n=reshape(Cavity_temp_cminus5n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_cminus20n=reshape(Cavity_temp_cminus20n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_cminus40n=reshape(Cavity_temp_cminus40n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_cminus50n=reshape(Cavity_temp_cminus50n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_controln=reshape(Cavity_temp_controln2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_cplus10n=reshape(Cavity_temp_cplus10n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_cplus5n=reshape(Cavity_temp_cplus5n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_cplus20n=reshape(Cavity_temp_cplus20n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_cminus30n=reshape(Cavity_temp_cminus30n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% Cavity_temp_cplus15n=reshape(Cavity_temp_cplus15n2,1,size(Cavity_temp_controln2,1)*size(Cavity_temp_controln2,2));
% 
%          
% Cavity_temp_cplus15n =  Cavity_temp_cplus15n(~Cavity_temp_cplus15n==0);
% Cavity_temp_cplus5n =  Cavity_temp_cplus5n(~Cavity_temp_cplus5n==0);
% Cavity_temp_cplus20n =  Cavity_temp_cplus20n(~Cavity_temp_cplus20n==0);
% Cavity_temp_cminus10n =  Cavity_temp_cminus10n(~Cavity_temp_cminus10n==0);
% Cavity_temp_cminus20n =  Cavity_temp_cminus20n(~Cavity_temp_cminus20n==0);
% Cavity_temp_cminus50n =  Cavity_temp_cminus50n(~Cavity_temp_cminus50n==0);
% Cavity_temp_cminus30n =  Cavity_temp_cminus30n(~Cavity_temp_cminus30n==0);
% Cavity_temp_cminus5n =  Cavity_temp_cminus5n(~Cavity_temp_cminus5n==0);
% Cavity_temp_cminus40n =  Cavity_temp_cminus40n(~Cavity_temp_cminus40n==0);
% Cavity_temp_controln =  Cavity_temp_controln(~Cavity_temp_controln==0);
% Cavity_temp_cplus10n =  Cavity_temp_cplus10n(~Cavity_temp_cplus10n==0);
% 
% 
% temp_control_q3 = median(Cavity_temp_controln(find(Cavity_temp_controln>median(Cavity_temp_controln))));
% 
% s_cminus20_q3 = median(Cavity_temp_cminus20n(find(Cavity_temp_cminus20n>median(Cavity_temp_cminus20n))));
% 
% s_cminus10_q3 = median(Cavity_temp_cminus10n(find(Cavity_temp_cminus10n>median(Cavity_temp_cminus10n))));
% 
% s_cplus15_q3 = median(Cavity_temp_cplus15n(find(Cavity_temp_cplus15n>median(Cavity_temp_cplus15n))));
% 
% s_cminus5_q3 = median(Cavity_temp_cminus5n(find(Cavity_temp_cminus5n>median(Cavity_temp_cminus5n))));
% 
% s_cminus40_q3 = median(Cavity_temp_cminus40n(find(Cavity_temp_cminus40n>median(Cavity_temp_cminus40n))));
% 
% s_cminus50_q3 = median(Cavity_temp_cminus50n(find(Cavity_temp_cminus50n>median(Cavity_temp_cminus50n))));
% 
% s_cplus10_q3 = median(Cavity_temp_cplus10n(find(Cavity_temp_cplus10n>median(Cavity_temp_cplus10n))));
% 
% s_cplus20_q3 = median(Cavity_temp_cplus20n(find(Cavity_temp_cplus20n>median(Cavity_temp_cplus20n))));
% 
% s_cminus30_q3 = median(Cavity_temp_cminus30n(find(Cavity_temp_cminus30n>median(Cavity_temp_cminus30n))));
% 
% s_cplus5_q3 = median(Cavity_temp_cplus5n(find(Cavity_temp_cplus5n>median(Cavity_temp_cplus5n))));
% 
% 
% err_exp1 = std(Cavity_temp_cminus50n2(Cavity_temp_cminus50n2==0))
% err_exp5 = std(Cavity_temp_controln2(Cavity_temp_controln2==0))
% err_exp10 = std(Cavity_temp_cminus10n2(Cavity_temp_cminus10n2==0))
% err_exp7 = std(Cavity_temp_cminus5n2(Cavity_temp_cminus5n2==0))
% % err_exp6 = std(min_cminus15);
% err_exp9 = std(Cavity_temp_cminus30n2(Cavity_temp_cminus30n2==0))
% err_exp3 = std(Cavity_temp_cplus20n2(Cavity_temp_cplus20n2==0))
% err_exp2 = std(Cavity_temp_cplus10n2(Cavity_temp_cplus10n2==0))
% err_exp4 = std(Cavity_temp_cplus5n2(Cavity_temp_cplus5n2==0))
% err_exp12 = std(Cavity_temp_cminus40n2(Cavity_temp_cminus40n2==0))
% err_exp13 = std(Cavity_temp_cplus15n2(Cavity_temp_cplus15n2==0))

Cavity_temp_controln2=0;           
Cavity_temp_cminus10n2=0;
Cavity_temp_cminus20n2=0;
Cavity_temp_cminus5n2=0;
Cavity_temp_cminus50n2=0;
Cavity_temp_cminus40n2=0;
Cavity_temp_cminus30n2=0;
Cavity_temp_cplus20n2=0;
Cavity_temp_cplus10n2=0;
Cavity_temp_cplus5n2=0;
Cavity_temp_cplus15n2=0;
Melt = 0;
Meltcplus5=0;
Meltcplus10=0;
Meltcplus20=0;
Meltcplus15=0;
Meltcminus5=0;
Meltcminus20=0;
Meltcminus10=0;
Meltcminus30=0;
Meltcminus50=0;
Meltcminus40=0;
for i = 1:Nx
    for j = 1:Ny
            if hFacC(i,j,kmax(i,j))>0
                 if XC(i,j)>-80 && XC(i,j)<-30
                         if YC(i,j)<-75
                            if  SHELFICEtopo(i,j)<0  
                             if Melt_control(i,j)<0    
                                Cavity_temp_controln2 = Cavity_temp_controln2+Cavity_temp_control(i,j,kmax(i,j))*RAC(i,j)*Melt_control(i,j); 
                                Melt = Melt + RAC(i,j)*Melt_control(i,j);
                             end
                             if Melt_cminus10(i,j)<0   
                                Meltcminus10=Meltcminus10+RAC(i,j)*Melt_cminus10(i,j);
                             
                                Cavity_temp_cminus10n2 = Cavity_temp_cminus10n2+Cavity_temp_cminus10(i,j,kmax(i,j))*RAC(i,j)*Melt_cminus10(i,j);
                             end
                             if Melt_cminus20(i,j)<0   
                              
                               Meltcminus20=Meltcminus20+RAC(i,j)*Melt_cminus20(i,j);
                              Cavity_temp_cminus20n2 = Cavity_temp_cminus20n2+Cavity_temp_cminus20(i,j,kmax(i,j))*RAC(i,j)*Melt_cminus20(i,j);
                             end
                             if Melt_cminus5(i,j)<0   
                              
                               Meltcminus5=Meltcminus5+RAC(i,j)*Melt_cminus5(i,j);
                               Cavity_temp_cminus5n2 = Cavity_temp_cminus5n2+Cavity_temp_cminus5(i,j,kmax(i,j))*RAC(i,j)*Melt_cminus5(i,j);
                             end
                             if Melt_cminus50(i,j)<0   
                                Meltcminus50=Meltcminus50+RAC(i,j)*Melt_cminus50(i,j);
                         
                                Cavity_temp_cminus50n2 = Cavity_temp_cminus50n2+Cavity_temp_cminus50(i,j,kmax(i,j))*RAC(i,j)*Melt_cminus50(i,j);
                             end
                             if Melt_cminus40(i,j)<0   
                                 Meltcminus40=Meltcminus40+RAC(i,j)*Melt_cminus40(i,j);
                           
                                Cavity_temp_cminus40n2 = Cavity_temp_cminus40n2+Cavity_temp_cminus40(i,j,kmax(i,j))*RAC(i,j)*Melt_cminus40(i,j);
                             end
                             if Melt_cminus30(i,j)<0   
                                Meltcminus30=Meltcminus30+RAC(i,j)*Melt_cminus30(i,j);
                           
                                Cavity_temp_cminus30n2 = Cavity_temp_cminus30n2 +Cavity_temp_cminus30(i,j,kmax(i,j))*RAC(i,j)*Melt_cminus30(i,j);
                             end
                             if Melt_cplus20(i,j)<0   
                                 Meltcplus20=Meltcplus20+RAC(i,j)*Melt_cplus20(i,j);
                             
                                Cavity_temp_cplus20n2 = Cavity_temp_cplus20n2+Cavity_temp_cplus20(i,j,kmax(i,j))*RAC(i,j)*Melt_cplus20(i,j);
                             end
                             if Melt_wplus10(i,j)<0   
                                Meltcplus10=Meltcplus10+RAC(i,j)*Melt_cplus10(i,j);
                              
                                Cavity_temp_cplus10n2 = Cavity_temp_cplus10n2+ Cavity_temp_cplus10(i,j,kmax(i,j))*RAC(i,j)*Melt_cplus10(i,j);
                             end
                             if Melt_cplus5(i,j)<0   
                                Meltcplus5=Meltcplus5+RAC(i,j)*Melt_cplus5(i,j);
                               
                                Cavity_temp_cplus5n2 = Cavity_temp_cplus5n2+Cavity_temp_cplus5(i,j,kmax(i,j))*RAC(i,j)*Melt_cplus5(i,j);
                             end
                             if Melt_cplus15(i,j)<0   
                               
                                Cavity_temp_cplus15n2 = Cavity_temp_cplus15n2+Cavity_temp_cplus15(i,j,kmax(i,j))*RAC(i,j)*Melt_cplus15(i,j);
                                                         
                                Meltcplus15=Meltcplus15+RAC(i,j)*Melt_cplus15(i,j);
                             end
                            end
                           
                         end
                
                 end
            end
        end
    
        
    
end
                           
                             
       



Cavity_salt_cminus10n=(Cavity_temp_cminus10n2/Meltcminus10);
Cavity_salt_cminus5n=(Cavity_temp_cminus5n2/Meltcminus5);
Cavity_salt_cminus20n=(Cavity_temp_cminus20n2/Meltcminus20);
Cavity_salt_cminus30n=(Cavity_temp_cminus30n2/Meltcminus30);
Cavity_salt_cminus40n=(Cavity_temp_cminus40n2/Meltcminus40);
Cavity_salt_cminus50n=(Cavity_temp_cminus50n2/Meltcminus50);
Cavity_salt_controln=Cavity_temp_controln2/Melt;
Cavity_salt_cplus10n=Cavity_temp_cplus10n2/Meltcplus10;
Cavity_salt_cplus5n=Cavity_temp_cplus5n2/Meltcplus5;
Cavity_salt_cplus20n=Cavity_temp_cplus20n2/Meltcplus20;
Cavity_salt_cplus15n=Cavity_temp_cplus15n2/Meltcplus15;


figure(1)
clf

 
c_s = [.5,.6,.7,.8,.9,.95,1,1.05, 1.1,1.15, 1.2];
w_s = [.5,.6,.7,.8,.9,.95,1,1.05,1.1,1.15,1.2];

% dataw = [s_wminus50_q3,s_wminus40_q3,s_wminus30_q3,s_wminus20_q3,s_wminus10_q3,s_wminus5_q3,temp_wcontrol_q3,s_wplus5_q3,s_wplus10_q3,s_wplus15_q3,s_wplus20_q3];
% datac = [s_cminus50_q3,s_cminus40_q3,s_cminus30_q3,s_cminus20_q3,s_cminus10_q3,s_cminus5_q3,temp_control_q3,s_cplus5_q3,s_cplus10_q3,s_cplus15_q3,s_cplus20_q3];

dataw = [Cavity_salt_wminus50n,Cavity_salt_wminus40n,Cavity_salt_wminus30n,Cavity_salt_wminus20n,Cavity_salt_wminus10n,Cavity_salt_wminus5n,Cavity_salt_wcontroln,Cavity_salt_wplus5n,Cavity_salt_wplus10n,Cavity_salt_wplus15n,Cavity_salt_wplus20n];
datac = [Cavity_salt_cminus50n,Cavity_salt_cminus40n,Cavity_salt_cminus30n,Cavity_salt_cminus20n,Cavity_salt_cminus10n,Cavity_salt_cminus5n,Cavity_salt_controln,Cavity_salt_cplus5n,Cavity_salt_cplus10n,Cavity_salt_cplus15n,Cavity_salt_cplus20n];


%%%%plot warm 
p1 = plot(w_s,dataw,'-or','LineWidth',2,'MarkerSize',10,'Markerfacecolor','r');
hold on
plot(w_s(7),dataw(7),'-sk','MarkerSize',15,'Markerfacecolor','k');
hold on
%%%% plot cold
p2 = plot(c_s,datac,'-ob','LineWidth',2,'MarkerSize',10,'Markerfacecolor','b');
hold on
plot(c_s(7),datac(7),'-sk','MarkerSize',15,'Markerfacecolor','k');
% hold on
% hold on
% errorbar(w_s(1),dataw(1),err_exp4w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% errorbar(w_s(2),dataw(2),err_exp13w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% 
% hold on
% errorbar(w_s(3),dataw(3),err_exp12w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% errorbar(w_s(4),dataw(4),err_exp8w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% errorbar(w_s(5),dataw(5),err_exp7w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% errorbar(w_s(6),dataw(6),err_exp9w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% errorbar(w_s(7),dataw(7),err_exp1w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% errorbar(w_s(8),dataw(8),err_exp11w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% 
% errorbar(w_s(9),dataw(9),err_exp5w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% errorbar(w_s(10),dataw(10),err_exp6w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% errorbar(w_s(11),dataw(11),err_exp2w,'-rs','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on
% 
% 
% errorbar(c_s(1),datac(1),err_exp1,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on
% errorbar(c_s(2),datac(2),err_exp12,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on
% errorbar(c_s(3),datac(3),err_exp11,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on
% errorbar(c_s(4),datac(4),err_exp9,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on
% errorbar(c_s(5),datac(5),err_exp10,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on
% errorbar(c_s(6),datac(6),err_exp7,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on

% errorbar(c_s(7),datac(7),err_exp5,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on
% errorbar(c_s(8),datac(8),err_exp4,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on
% 
% errorbar(c_s(9),datac(9),err_exp2,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on
% errorbar(c_s(10),datac(10),err_exp13,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on
% errorbar(c_s(11),datac(11),err_exp3,'-bs','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% hold on



set(gcf,'unit','pixel','position',[0 0 1600 1200])
title('Wind Perturbation ($\chi$) vs Melt-Weighted Bottom Cavity Salinity (psu)','interpreter','latex','FontSize',20);
ylabel('Melt-Weighted Bottom Cavity Salinity (psu)','interpreter','latex','FontSize',16)
xlabel('Wind Perturbation ($\chi$)','interpreter','latex','FontSize',16);
hold on

% y = 1*ones(1,1500);
% k2 = plot(y,1:1500,'m:','linewidth',2);
% hold on


% Cp = 4e8;
% Cm = 9.6e8;
% Sigma_p0 = -1.4e8;
% Sigma_m0 = 5.1e7;
% alpha = 6.6e-5;
% beta = 7.6e-4;
% Psi = 2.2e6;
% S0 = 34;
% rho0 = 1000;
% Saasw = 34.2;
% Scdw = 34.35;
% Tcdw = -1;
% Thssw = -2.08;
% Tf = -2.38;
% V_ref = 2.7;
%  
%  
% V_lo = (1/Cp) * ( rho0*Psi*( Scdw - Saasw + (alpha/beta)*(Thssw-Tcdw)) ...
%                 + Cm*(Thssw-Tf) - Sigma_p0 + Sigma_m0);
% V_hi = (1/Cp) * ( rho0*Psi*( Scdw - Saasw + (alpha/beta)*(Thssw-Tcdw)) ...
%                 + Cm*(Tcdw-Tf) - Sigma_p0 + Sigma_m0);
%               
% (1/Cp) * ( rho0*Psi*( Scdw - Saasw) ) 
% (1/Cp) * ( rho0*Psi*(alpha/beta)*(Thssw-Tcdw) ) 
% (1/Cp) * Cm*(Thssw-Tf) 
% (1/Cp) * (- Sigma_p0 + Sigma_m0)
%  
% V_lo/V_ref
% V_hi/V_ref
%  
% Vwind = V_ref*[.5:.1:1.2];
% Shssw = Saasw + (Sigma_p0 + Cp*Vwind - Sigma_m0 - Cm*(Thssw-Tf)) / (rho0*Psi)
%               
% hold on
% m=plot([.5:.1:1.2],Shssw,'--k','linewidth',3);
% 
% 
% ylim([34.25 34.6]);
% legend([p1,p2],{'Warm State Initialization','Cold State Initialization'},'interpreter','latex','fontsize',15,'location','northwest')
% set(gca,'FontSize',17');
% 
