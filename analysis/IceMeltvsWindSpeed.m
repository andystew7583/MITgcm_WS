%%%%%%%Checking to see if V_i scales with V_meridional

%%%%%Goal: look at average ice speeds over SW weddell and compare to the
%%%%%Wind Speeds


%%%%%% First we check the magnitude of the wind vectors with time




%%% Frequency of diagnostic output

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


inputpath1 = '/data3/MITgcm_WS/experiments/s_hist38_c_minus50/input';
inputpath2 = '/data3/MITgcm_WS/experiments/s_hist31_c_meridplus10/input';
inputpath3 = '/data3/MITgcm_WS/experiments/s_hist1_meridplus20/input';
inputpath4 = '/data3/MITgcm_WS/experiments/s_hist32_c_meridplus5/input';
inputpath5 = '/data3/MITgcm_WS/experiments/a_3445_20boundary/input';
inputpath6 = '/data3/MITgcm_WS/experiments/s_hist37_c_minus15_2/input';
inputpath7 = '/data3/MITgcm_WS/experiments/s_hist34_cminus5_2/input';
inputpath8 = '/data3/MITgcm_WS/experiments/s_hist19_cminus2_2/input';
inputpath9 = '/data3/MITgcm_WS/experiments/s_hist36_c_minus20_2/input';
inputpath10 = '/data3/MITgcm_WS/experiments/s_hist35_cminus10_2/input';
inputpath11='/data3/MITgcm_WS/experiments/s_hist41_c_minus30/input';
inputpath12='/data3/MITgcm_WS/experiments/s_hist41_c_minus40/input';
inputpath13 = '/data3/MITgcm_WS/experiments/s_hist45_c_meridplus15/input';

inputpath1w = '/data3/MITgcm_WS/experiments/a_34_20boundary/input';
inputpath2w = '/data3/MITgcm_WS/experiments/s_hist4_warm_meridplus20/input';
inputpath3w = '/data3/MITgcm_WS/experiments/s_hist18_warm_windsplus2/input';
inputpath4w = '/data3/MITgcm_WS/experiments/w_hist42_wminus50/input';
inputpath5w = '/data3/MITgcm_WS/experiments/s_hist40_w_plus10/input';
inputpath6w = '/data3/MITgcm_WS/experiments/s_hist24_w_mplus15_2/input';
inputpath7w = '/data3/MITgcm_WS/experiments/s_hist32_w_meridminus10/input';
inputpath8w = '/data3/MITgcm_WS/experiments/s_hist32_w_meridminus20/input';
inputpath9w = '/data3/MITgcm_WS/experiments/s_hist32_w_meridminus5/input';
inputpath10w = '/data3/MITgcm_WS/experiments/s_hist24_w_mplus15_2/input';
inputpath11w = '/data3/MITgcm_WS/experiments/s_hist43_w_meridplus5/input';
inputpath12w = '/data3/MITgcm_WS/experiments/s_hist44_w_meridminus30/input';
inputpath13w =  '/data3/MITgcm_WS/experiments/s_hist46_w_minus40/input';




% Nyears = (nDump_f-nDump_s)/12;

%%%% load Sea Ice Velocity Data
tmin = 9*86400*360;
tmax = 18*86400*360;
Melt_c = readIters(exppath5,'SHIfwFlx',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);

Melt_cminus5 = readIters(exppath7,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);




Melt_cminus20 = readIters(exppath9,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


Melt_cminus15 = readIters(exppath6,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);




Melt_cplus5 = readIters(exppath4,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


Melt_cplus10 = readIters(exppath2,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Melt_cplus15 = readIters(exppath13,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);







tmin = 7*86400*360;
tmax = 16*86400*360;
Melt_cminus30 = readIters(exppath11,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

tmin = 9*86400*360;
tmax = 18*86400*360;
Melt_cminus10 = readIters(exppath10,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Melt_cplus20 = readIters(exppath3,'SHIfwFlx',dumpIters4,deltaT_4,tmin,tmax,Nx,Ny,1);




%%%%%% Warm experiments
tmin = 6*86400*360;
tmax = 7*86400*360;

Melt_w = readIters(exppath1w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


Melt_wminus5 = readIters(exppath9w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Melt_wminus10 = readIters(exppath7w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

Melt_wminus20 = readIters(exppath8w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


Melt_wplus10 = readIters(exppath5w,'SHIfwFlx',dumpIters4,deltaT_4,tmin,tmax,Nx,Ny,1);


tmin = 28*86400*360;
tmax = 37*86400*360;

Melt_wplus20 = readIters(exppath2w,'SHIfwFlx',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,1);


tmin = 19*86400*360;
tmax = 28*86400*360;

Melt_cminus50 = readIters(exppath1,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Melt_cminus40 = readIters(exppath12,'SHIfwFlx',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,1);

tmin = 19*86400*360;
tmax = 27*86400*360;
Melt_wplus15 = readIters(exppath10w,'SHIfwFlx',dumpIters4,deltaT_4,tmin,tmax,Nx,Ny,1);


tmin = 9*86400*360;
tmax = 18*86400*360;


Melt_wplus5 = readIters(exppath11w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Melt_wminus50 = readIters(exppath4w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

tmin = 7*86400*360;
tmax = 16*86400*360;

Melt_wminus30 = readIters(exppath12w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);
Melt_wminus40 = readIters(exppath13w,'SHIfwFlx',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


% tmax = 12*86400*360;




days=3287;
%%%%%cold
% 
Vwind_c = NaN(Nx,Ny,3287);
Vwind_cminus50=NaN(Nx,Ny,3287);
% Vwind_cminus2 = NaN(Nx,Ny,3287);
Vwind_cminus5= NaN(Nx,Ny,3287);
Vwind_cminus10= NaN(Nx,Ny,3287);
Vwind_cminus15= NaN(Nx,Ny,3287);
Vwind_cminus20= NaN(Nx,Ny,3287);
Vwind_cplus5= NaN(Nx,Ny,3287);
Vwind_cplus10= NaN(Nx,Ny,3287);
Vwind_cplus20= NaN(Nx,Ny,3287);
Vwind_cminus30= NaN(Nx,Ny,3287);
Vwind_cminus40= NaN(Nx,Ny,3287);
Vwind_cplus15= NaN(Nx,Ny,3287);

%%%%%warm
Vwind_wminus50= NaN(Nx,Ny,3287);
Vwind_wminus10= NaN(Nx,Ny,3287);
Vwind_wplus15= NaN(Nx,Ny,3287);
Vwind_wminus20= NaN(Nx,Ny,3287);
Vwind_wminus5= NaN(Nx,Ny,3287);
Vwind_wplus2= NaN(Nx,Ny,3287);
Vwind_wplus10= NaN(Nx,Ny,3287);
Vwind_wplus20= NaN(Nx,Ny,3287);
Vwind_wc = NaN(Nx,Ny,3287);
Vwind_wminus30= NaN(Nx,Ny,3287);
Vwind_wminus40= NaN(Nx,Ny,3287);
Vwind_wplus5= NaN(Nx,Ny,3287);




fidv = fopen(fullfile(inputpath1,'vwindfile.bin'),'r','b');
fid2v = fopen(fullfile(inputpath2,'vwindfile.bin'),'r','b');
fid3v = fopen(fullfile(inputpath3,'vwindfile.bin'),'r','b');
fid4v = fopen(fullfile(inputpath4,'vwindfile.bin'),'r','b');
fid5v = fopen(fullfile(inputpath5,'vwindfile.bin'),'r','b');
fid6v = fopen(fullfile(inputpath6,'vwindfile.bin'),'r','b');
fid7v = fopen(fullfile(inputpath7,'vwindfile.bin'),'r','b');
fid8v = fopen(fullfile(inputpath8,'vwindfile.bin'),'r','b');
fid9v = fopen(fullfile(inputpath9,'vwindfile.bin'),'r','b');
fid10v= fopen(fullfile(inputpath10,'vwindfile.bin'),'r','b');
fid11v= fopen(fullfile(inputpath11,'vwindfile.bin'),'r','b');
fid12v= fopen(fullfile(inputpath12,'vwindfile.bin'),'r','b');
fid13v= fopen(fullfile(inputpath13,'vwindfile.bin'),'r','b');



fid1wv = fopen(fullfile(inputpath1w,'vwindfile.bin'),'r','b');
fid2wv = fopen(fullfile(inputpath2w,'vwindfile.bin'),'r','b');
fid3wv = fopen(fullfile(inputpath3w,'vwindfile.bin'),'r','b');
fid4wv = fopen(fullfile(inputpath4w,'vwindfile.bin'),'r','b');
fid5wv = fopen(fullfile(inputpath5w,'vwindfile.bin'),'r','b');
fid6wv = fopen(fullfile(inputpath6w,'vwindfile.bin'),'r','b');
fid7wv = fopen(fullfile(inputpath7w,'vwindfile.bin'),'r','b');
fid8wv = fopen(fullfile(inputpath8w,'vwindfile.bin'),'r','b');
fid9wv = fopen(fullfile(inputpath9w,'vwindfile.bin'),'r','b');
fid10wv = fopen(fullfile(inputpath10w,'vwindfile.bin'),'r','b');
fid11wv = fopen(fullfile(inputpath11w,'vwindfile.bin'),'r','b');
fid12wv = fopen(fullfile(inputpath12w,'vwindfile.bin'),'r','b');
fid13wv = fopen(fullfile(inputpath13w,'vwindfile.bin'),'r','b');




for k=1:days

% 
  Vwind_c(:,:,k) = fread(fid5v,[Nx Ny],'real*8');
  Vwind_cminus30(:,:,k) = fread(fid11v,[Nx Ny],'real*8');
  Vwind_cminus5(:,:,k) = fread(fid7v,[Nx Ny],'real*8');
  Vwind_cminus10(:,:,k) = fread(fid10v,[Nx Ny],'real*8');
  Vwind_cminus20(:,:,k) = fread(fid9v,[Nx Ny],'real*8');
  Vwind_cminus15(:,:,k) = fread(fid6v,[Nx Ny],'real*8');
  Vwind_cminus50(:,:,k) = fread(fidv,[Nx Ny],'real*8');
  Vwind_cplus5(:,:,k) = fread(fid4v,[Nx Ny],'real*8');
  Vwind_cplus10(:,:,k) = fread(fid2v,[Nx Ny],'real*8');
  Vwind_cplus20(:,:,k) = fread(fid3v,[Nx Ny],'real*8');
  Vwind_cminus40(:,:,k) = fread(fid12v,[Nx Ny],'real*8');
  Vwind_cplus15(:,:,k) = fread(fid13v,[Nx Ny],'real*8');

  
  Vwind_wc(:,:,k) = fread(fid1wv,[Nx Ny],'real*8');
  Vwind_wplus5(:,:,k) = fread(fid11wv,[Nx Ny],'real*8');
  Vwind_wplus10(:,:,k) = fread(fid5wv,[Nx Ny],'real*8');
  Vwind_wplus20(:,:,k) = fread(fid2wv,[Nx Ny],'real*8');
  Vwind_wplus15(:,:,k) = fread(fid6wv,[Nx Ny],'real*8');
  Vwind_wminus10(:,:,k) = fread(fid7wv,[Nx Ny],'real*8');
  Vwind_wminus20(:,:,k) = fread(fid8wv,[Nx Ny],'real*8');
  Vwind_wminus50(:,:,k) = fread(fid4wv,[Nx Ny],'real*8'); 
  Vwind_wminus30(:,:,k) = fread(fid12wv,[Nx Ny],'real*8'); 
  Vwind_wminus5(:,:,k) = fread(fid9wv,[Nx Ny],'real*8'); 
  Vwind_wminus40(:,:,k) = fread(fid13wv,[Nx Ny],'real*8'); 
 
end
fclose(fidv);
fclose(fid2v);
fclose(fid4v);
fclose(fid3v);
fclose(fid5v);
fclose(fid6v);
fclose(fid7v);
fclose(fid9v);
fclose(fid10v);
fclose(fid11v);
fclose(fid12v);
fclose(fid13v);


fclose(fid1wv);
fclose(fid2wv);
fclose(fid4wv);
fclose(fid5wv);
fclose(fid6wv);
fclose(fid7wv);
fclose(fid8wv);
fclose(fid9wv);
fclose(fid10wv);
fclose(fid11wv);
fclose(fid12wv);
fclose(fid13wv);






Vwind_c = squeeze(mean(Vwind_c,3));
Vwind_cminus50 = squeeze(mean(Vwind_cminus50,3));
Vwind_cminus5= squeeze(mean(Vwind_cminus5,3));
Vwind_cminus10= squeeze(mean(Vwind_cminus10,3));
Vwind_cminus15= squeeze(mean(Vwind_cminus15,3));
Vwind_cminus20= squeeze(mean(Vwind_cminus20,3));
Vwind_cplus5= squeeze(mean(Vwind_cplus5,3));
Vwind_cplus10= squeeze(mean(Vwind_cplus10,3));
Vwind_cplus20= squeeze(mean(Vwind_cplus20,3));
Vwind_cplus15= squeeze(mean(Vwind_cplus15,3));
Vwind_cminus30= squeeze(mean(Vwind_cminus30,3));
Vwind_cminus40= squeeze(mean(Vwind_cminus40,3));





Vwind_wminus5= squeeze(mean(Vwind_wminus5,3));
Vwind_wminus10= squeeze(mean(Vwind_wminus10,3));
Vwind_wplus15= squeeze(mean(Vwind_wplus15,3));
Vwind_wminus20= squeeze(mean(Vwind_wminus20,3));
Vwind_wplus5= squeeze(mean(Vwind_wplus5,3));
Vwind_wminus40= squeeze(mean(Vwind_wminus40,3));
Vwind_wplus10= squeeze(mean(Vwind_wplus10,3));
Vwind_wplus20= squeeze(mean(Vwind_wplus20,3));
Vwind_wc = squeeze(mean(Vwind_wc,3));
Vwind_wminus50= squeeze(mean(Vwind_wminus50,3));
Vwind_wminus30= squeeze(mean(Vwind_wminus30,3));


Vol =0;
    
Melt_cn=0;
Melt_cminus5n=0;
Melt_cminus20n=0;
Melt_cminus40n=0;
Melt_cminus30n=0;
Melt_cminus10n=0;
Melt_cplus5n=0;
Melt_cplus10n=0;
Melt_cplus15n=0;
Melt_cplus20n=0;
Melt_cminus50n =0;

Melt_wn=0;
Melt_wminus5n=0;
Melt_wminus20n=0;
Melt_wminus40n=0;
Melt_wminus30n=0;
Melt_wminus10n=0;
Melt_wplus5n=0;
Melt_wplus10n=0;
Melt_wplus15n=0;
Melt_wplus20n=0;
Melt_wminus50n=0;
U_c = 0;
U_cminus30 = 0;
U_cminus40 = 0;
U_cminus50= 0;
U_cminus= 0;
U_cminus10=0;
U_cminus20=0;
U_cplus5=0;
U_cplus10=0;
U_cplus15=0;
U_cplus20=0;


U_w = 0;
U_wminus30 = 0;
U_wminus40 = 0;
U_wminus50= 0;
U_wminus= 0;
U_wminus10=0;
U_wminus20=0;
U_wplus5=0;
U_wplus10=0;
U_wplus15=0;
U_wplus20=0;
U_wminus5 = 0;
U_cminus5 = 0;





 %%% Integrate melt over the FRIS
    for i=1:Nx
        for j=1:Ny
            if XC(i,j)>-80 && XC(i,j) <-20
                if YC(i,j) <-75
                 if SHELFICEtopo(i,j)< 0

                        Melt_cn= Melt_cn+ Melt_c(i,j)*RAC(i,j);
                        Melt_cminus5n= Melt_cminus5n + Melt_cminus5(i,j)*RAC(i,j);
                        
                        Melt_cminus30n= Melt_cminus30n+Melt_cminus30(i,j)*RAC(i,j);
                        Melt_cminus40n= Melt_cminus40n+Melt_cminus40(i,j)*RAC(i,j);               
                        Melt_cminus10n= Melt_cminus10n+Melt_cminus10(i,j)*RAC(i,j);
                        Melt_cminus20n= Melt_cminus20n+Melt_cminus20(i,j)*RAC(i,j);
                        Melt_cminus50n= Melt_cminus50n+Melt_cminus50(i,j)*RAC(i,j);               
                        Melt_cplus5n= Melt_cplus5n+Melt_cplus5(i,j)*RAC(i,j);
                        Melt_cplus10n= Melt_cplus10n+Melt_cplus10(i,j)*RAC(i,j);
                        Melt_cplus20n= Melt_cplus20n+Melt_cplus20(i,j)*RAC(i,j);
                        Melt_cplus15n= Melt_cplus15n+Melt_cplus15(i,j)*RAC(i,j);

                        Melt_wn= Melt_wn+Melt_w(i,j)*RAC(i,j);
                        Melt_wminus5n= Melt_wminus5n+Melt_wminus5(i,j)*RAC(i,j);
                        Melt_wminus10n= Melt_wminus10n+Melt_wminus10(i,j)*RAC(i,j);
                        Melt_wminus20n= Melt_wminus20n+Melt_wminus20(i,j)*RAC(i,j);
                        Melt_wplus15n= Melt_wplus15n+Melt_wplus15(i,j)*RAC(i,j);               
                        Melt_wplus5n= Melt_wplus5n+Melt_wplus5(i,j)*RAC(i,j);
                        Melt_wplus10n= Melt_wplus10n+Melt_wplus10(i,j)*RAC(i,j);
                        Melt_wplus20n= Melt_wplus20n+Melt_wplus20(i,j)*RAC(i,j);                        
                        Melt_wminus30n= Melt_wminus30n+Melt_wminus30(i,j)*RAC(i,j);
                        Melt_wminus40n= Melt_wminus40n+Melt_wminus40(i,j)*RAC(i,j);                          
                        Melt_wminus50n= Melt_wminus50n+Melt_wminus50(i,j)*RAC(i,j);                          
                        
                 end
                end
            end
                
        end
    end
    
   for i=1:Nx
        for j=1:Ny            
                 if XC(i,j) >-60 && XC(i,j) <-50
                         if YC(i,j) > -76 && YC(i,j) <-74
                                 if  hFacC(i,j,1) >0

                                    U_c = U_c + (Vwind_c(i,j)) * hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                                                     
                                    
                                    U_cminus20 = U_cminus20 + (Vwind_cminus20(i,j)) * hFacC(i,j,1)*RAC(i,j)*delR(1);                                   
                                    
                                    U_cminus10 = U_cminus10 + (Vwind_cminus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                     
                                    
                                    U_cminus5 = U_cminus5 + (Vwind_cminus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    
                                    U_cminus50 = U_cminus50 + (Vwind_cminus50(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                                                      
                                    
                                    U_cplus5 = U_cplus5 +(Vwind_cplus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    
                                    U_cplus10 = U_cplus10+(Vwind_cplus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                                                       
                                    U_cplus20 = U_cplus20 +(Vwind_cplus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                                    
                                    U_cplus15 = U_cplus15+(Vwind_cplus15(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                                                 
                                    U_cminus40 = U_cminus40 +(Vwind_cminus40(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);   
                                    
                                    U_cminus30 = U_cminus30 +(Vwind_cminus30(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    
                                    %Warm
                                    
                                    
                                    U_wplus15 = U_wplus15+ (Vwind_wplus15(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                       
                                    
                                    U_wminus20 = U_wminus20+(Vwind_wminus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                   
                                    
                                    U_wminus10 = U_wminus10+(Vwind_wminus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                     
%                                    
                                    U_wminus5 = U_wminus5+(Vwind_wminus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 

                                    U_w =U_w +(Vwind_wc(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 


                                    
                                    U_wminus50 = U_wminus50+(Vwind_wminus50(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    
                                    U_wplus10 =U_wplus10 +(Vwind_wplus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                                                       
                                    U_wplus20 = U_wplus20+ (Vwind_wplus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                                                          
                                    U_wplus5 = U_wplus5 +(Vwind_wplus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                                                                                           
                                                                             
                                    U_wminus30 = U_wminus30 +(Vwind_wminus30(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                     
                            
                                    U_wminus40 = U_wminus40+(Vwind_wminus40(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    Vol = Vol+ hFacC(i,j,1)*RAC(i,j)*delR(1);
                         
                         
                         
                         
                         
                         
                                 end
                         
                         end
                end
        
        end
    end
  
U_c = U_c/Vol;
U_cminus30 = U_cminus30/Vol;
U_cminus40 = U_cminus40/Vol;
U_cminus50= U_cminus50/Vol;
U_cminus10=U_cminus10/Vol;
U_cminus20=U_cminus20/Vol;
U_cplus5=U_cplus5/Vol;
U_cminus5 = U_cminus5/Vol;
U_cplus10=U_cplus10/Vol;
U_cplus15=U_cplus15/Vol;
U_cplus20=U_cplus20/Vol;


U_w = U_w/Vol;
U_wminus30 = U_wminus30/Vol;
U_wminus40 = U_wminus40/Vol;
U_wminus50= U_wminus50/Vol;
U_wminus10=U_wminus10/Vol;
U_wminus20=U_wminus20/Vol;
U_wplus5=U_wplus5/Vol;
U_wminus5 = U_wminus5/Vol;
U_wplus10=U_wplus10/Vol;
U_wplus15=U_wplus15/Vol;
U_wplus20=U_wplus20/Vol;





Melt_cn = -(Melt_cn*86400*365*1e-12);  % convrt to GT/year;  
Melt_cminus5n = -(Melt_cminus5n*86400*365*1e-12);  % convrt to GT/year;  
Melt_cminus40n = -(Melt_cminus40n*86400*365*1e-12);  % convrt to GT/year;  
Melt_cminus30n = -(Melt_cminus30n*86400*365*1e-12);  % convrt to GT/year;  
Melt_cminus10n = -(Melt_cminus10n*86400*365*1e-12);  % convrt to GT/year;  
Melt_cminus20n = -(Melt_cminus20n*86400*365*1e-12);  % convrt to GT/year;  
Melt_cminus50n = -(Melt_cminus50n*86400*365*1e-12);  % convrt to GT/year;  
Melt_cplus5n = -(Melt_cplus5n*86400*365*1e-12);  % convrt to GT/year;  
Melt_cplus10n = -(Melt_cplus10n*86400*365*1e-12);  % convrt to GT/year;  
Melt_cplus20n = -(Melt_cplus20n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wn = -(Melt_wn*86400*365*1e-12);  % convrt to GT/year;  
Melt_wminus30n = -(Melt_wminus30n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wplus5n = -(Melt_wplus5n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wplus10n = -(Melt_wplus10n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wplus15n = -(Melt_wplus15n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wplus20n = -(Melt_wplus20n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wminus5n = -(Melt_wminus5n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wminus10n = -(Melt_wminus10n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wminus20n = -(Melt_wminus20n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wminus40n = -(Melt_wminus40n*86400*365*1e-12);  % convrt to GT/year;  
Melt_wminus50n = -(Melt_wminus50n*86400*365*1e-12);  % convrt to GT/year;  

%%%%make plot

figure(1);
clf;
hold on


plot(U_c,Melt_cn,'-ok','markerfacecolor','b','markersize',10);
hold on
plot(U_cminus50,Melt_cminus50n,'-ok','markerfacecolor','b','markersize',10);
hold on

plot(U_cminus20,Melt_cminus20n,'-ok','markerfacecolor','b','markersize',10);
hold on

plot(U_cminus5,Melt_cminus5n,'-ok','markerfacecolor','b','markersize',10);
hold on

plot(U_cminus10,Melt_cminus10n,'-ok','markerfacecolor','b','markersize',10);
hold on

plot(U_cminus30,Melt_cminus30n,'-ok','markerfacecolor','b','markersize',10);
hold on

plot(U_cminus40,Melt_cminus40n,'-ok','markerfacecolor','b','markersize',10);
hold on

plot(U_cplus5,Melt_cplus5n,'-ok','markerfacecolor','b','markersize',10);
hold on

plot(U_cplus10,Melt_cplus10n,'-ok','markerfacecolor','b','markersize',10);
hold on

plot(U_cplus20,Melt_cplus20n,'-ok','markerfacecolor','b','markersize',10);
hold on

plot(U_w,Melt_wn,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wminus5,Melt_wminus5n,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wminus20,Melt_wminus20n,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wminus10,Melt_wminus10n,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wminus30,Melt_wminus30n,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wplus5,Melt_wplus5n,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wplus10,Melt_wplus10n,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wplus20,Melt_wplus20n,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wplus15,Melt_wplus15n,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wminus40,Melt_wminus40n,'-or','markerfacecolor','r','markersize',10);
hold on

plot(U_wminus50,Melt_wminus50n,'-or','markerfacecolor','r','markersize',10);
hold on



xlabel('Mean Offshore V$_{Wind}$($m/s$) ','Interpreter','latex','FontSize',16)

ylabel('Mean Integrated Shelf Ice Melt','Interpreter','latex','FontSize',16)

title('Avg V$_{wind}$ Speed ($m/s$) over Ronne Polynya vs Integrated FRIS Melt (Gt)','FontSize',20,'interpreter','latex');


ws_w = [U_wminus50,U_wminus40,U_wminus30,U_wminus20,U_wminus10,U_wminus5,U_w,U_wplus5,U_wplus10,U_wplus20];
Is_w = [Melt_wminus50n,Melt_wminus40n,Melt_wminus30n,Melt_wminus20n,Melt_wminus10n,Melt_wminus5n,Melt_wn,Melt_wplus5n,Melt_wplus10n,Melt_wplus20n];

h = fitlm(Is_w,ws_w);

ws_c = [U_cminus50,U_cminus40,U_cminus30,U_cminus20,U_cminus10,U_cminus5,U_c,U_cplus5,U_cplus10,U_cplus20];
Is_c = [Melt_cminus50n,Melt_cminus40n,Melt_cminus30n,Melt_cminus20n,Melt_cminus10n,Melt_cminus5n,Melt_cn,Melt_cplus5n,Melt_cplus10n,Melt_cplus20n];

z = fitlm(Is_c,ws_c);
