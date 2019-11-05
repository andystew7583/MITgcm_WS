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




exppath1 = '/data3/MITgcm_WS/experiments/s_hist38_c_minus50';
exppath2 = '/data3/MITgcm_WS/experiments/s_hist31_c_meridplus10';
exppath3 = '/data3/MITgcm_WS/experiments/s_hist1_meridplus20';
exppath4 = '/data3/MITgcm_WS/experiments/s_hist32_c_meridplus5';
exppath5 = '/data3/MITgcm_WS/experiments/a_3445_20boundary';
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

inputpath1w = '/data3/MITgcm_WS/experiments/n_342/input';
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


% Nyears = (nDump_f-nDump_s)/12;

%%%% load Sea Ice Velocity Data
tmin = 9*86400*360;
tmax = 18*86400*360;
SIVice_c = readIters(exppath5,'SIvice',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,1);
tmin = 9*86400*360;
tmax = 18*86400*360;
SIVice_w = readIters(exppath1w,'SIvice',dumpIters4,deltaT_4,tmin,tmax,Nx,Ny,1);


%%%% load Sea Ice Velocity Data
tmin = 9*86400*360;
tmax = 18*86400*360;
SIVice_cminus5 = readIters(exppath7,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


SIVice_cminus10 = readIters(exppath10,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


SIVice_cminus20 = readIters(exppath9,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


SIVice_cminus15 = readIters(exppath6,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);




SIVice_cplus5 = readIters(exppath4,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


SIVice_cplus10 = readIters(exppath2,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

SIVice_cplus15 = readIters(exppath13,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


SIVice_cplus20 = readIters(exppath3,'SIvice',dumpIters4,deltaT_4,tmin,tmax,Nx,Ny,1);




tmin = 9*86400*360;
tmax = 18*86400*360;
SIVice_cminus30 = readIters(exppath11,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);






%%%%%% Warm experiments
tmin = 9*86400*360;
tmax = 18*86400*360;


SIVice_wplus5 = readIters(exppath11w,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

SIVice_wminus50 = readIters(exppath4w,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

SIVice_wminus5 = readIters(exppath9w,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

SIVice_wminus10 = readIters(exppath7w,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

SIVice_wminus20 = readIters(exppath8w,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


SIVice_wplus10 = readIters(exppath5w,'SIvice',dumpIters4,deltaT_4,tmin,tmax,Nx,Ny,1);

tmin = 28*86400*360;
tmax = 37*86400*360;
SIVice_wplus20 = readIters(exppath2w,'SIvice',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,1);

tmin = 19*86400*360;
tmax = 28*86400*360;
SIVice_cminus40 = readIters(exppath12,'SIvice',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,1);

SIVice_cminus50 = readIters(exppath1,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);


tmin = 19*86400*360;
tmax = 27*86400*360;

SIVice_wplus15 = readIters(exppath10w,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

tmin = 9*86400*360;
tmax = 18*86400*360;

SIVice_wminus30 = readIters(exppath12w,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);

SIVice_wminus40 = readIters(exppath13w,'SIvice',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,1);



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


I_c = 0;
U_c = 0;
U_cminus30 = 0;
U_cminus40 = 0;
U_cminus50= 0;
U_cminus10=0;
U_cminus20=0;
U_cplus5=0;
U_cplus10=0;
U_cplus15=0;
U_cplus20=0;
I_cminus30 = 0;
I_cminus40 = 0;
I_cminus50= 0;
I_cminus10=0;
I_cminus20=0;
I_cplus5=0;
I_cplus10=0;
I_cplus15=0;
I_cplus20=0;


I_w = 0;
U_w = 0;
U_wminus30 = 0;
U_wminus40 = 0;
U_wminus50= 0;
U_wminus10=0;
U_wminus20=0;
U_wplus5=0;
U_wplus10=0;
U_wplus15=0;
U_wplus20=0;
I_wminus30 = 0;
I_wminus40 = 0;
I_wminus50= 0;
I_wminus10=0;
I_wminus20=0;
I_wplus5=0;
I_wplus10=0;
I_wplus15=0;
I_wplus20=0;
U_wminus5 = 0;
I_wminus5 = 0;
U_cminus5 = 0;
I_cminus5 = 0;


Vol = 0;
Vol2 = 0;
for i =1:Nx
%%%%% Finding average grid cell salinity
    for j = 1:Ny
                       if start(i,j)==1

                                    I_c =  I_c + (SIVice_c(i,j) * hFacC(i,j,1)*RAC(i,j)*delR(1)); 
                                
                                    I_cminus10 = I_cminus10 +  (SIVice_cminus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);


                                    
                                    I_cminus20 =  I_cminus20+ (SIVice_cminus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    
                                    
                                    I_cminus5 = I_cminus5 + (SIVice_cminus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);
                                    
                                    I_cminus50 =  I_cminus50+ (SIVice_cminus50(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);
                                                                      
                                    
                                    I_cplus5 =  I_cplus5+ (SIVice_cplus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    
                                    I_cplus10 =  I_cplus10+(SIVice_cplus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                                                       
                                    I_cplus20 = I_cplus20 +(SIVice_cplus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                                    

                                    I_cplus15 = I_cplus15+(SIVice_cplus15(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    
                                    
                                    I_cminus40 =  I_cminus40+(SIVice_cminus40(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    
                                    I_cminus30=  I_cminus30 +(SIVice_cminus30(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 

                                    %Warm
                                    
                                    
                                    I_wplus15 =  I_wplus15+(SIVice_wplus15(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);   
                                    
                                    I_wminus20 = I_wminus20+(SIVice_wminus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);
                                    
                                    I_wminus10=  I_wminus10+(SIVice_wminus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
%                                    
                                    I_wminus5 = I_wminus5 +(SIVice_wminus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 

                                    I_w =  I_w+(SIVice_w(i,j)* hFacC(i,j,1)*RAC(i,j)*delR(1)); 


                                    
                                    I_wminus50 = I_wminus50 +(SIVice_wminus50(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    
                                    I_wplus10 =I_wplus10 +(SIVice_wplus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                                                       
                                    I_wplus20 = I_wplus20 + (SIVice_wplus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                                                          
                                    I_wplus5 = I_wplus5 +(SIVice_wplus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                                                             
                                    I_wminus30=  I_wminus30+(SIVice_wminus30(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                            
                                    I_wminus40=  I_wminus40+(SIVice_wminus40(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    
                                    Vol = Vol + hFacC(i,j,1)*RAC(i,j)*delR(1);
                                    


                                    U_c = U_c + (Vwind_c(i,j)) * hFacC(i,j,1)*RAC(i,j)*delR(1); 
    
                                    U_cminus30 = U_cminus30 +(Vwind_cminus30(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);  

                                
                                    U_cminus20 = U_cminus20 + (Vwind_cminus20(i,j)) * hFacC(i,j,1)*RAC(i,j)*delR(1);                                   
                                    U_cminus10 = U_cminus10 + (Vwind_cminus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                     
                                    U_cminus5 = U_cminus5 + (Vwind_cminus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    U_cminus50 = U_cminus50 + (Vwind_cminus50(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    U_cplus5 = U_cplus5 +(Vwind_cplus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    U_cplus10 = U_cplus10+(Vwind_cplus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    U_cplus20 = U_cplus20 +(Vwind_cplus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    U_cplus15 = U_cplus15+(Vwind_cplus15(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    U_wplus15 = U_wplus15+ (Vwind_wplus15(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                       
                                    U_wminus20 = U_wminus20+(Vwind_wminus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                   
                                    U_cminus40 = U_cminus40 +(Vwind_cminus40(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);   
                                    U_wminus10 = U_wminus10+(Vwind_wminus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                     
                                    U_wminus50 = U_wminus50+(Vwind_wminus50(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    U_w =U_w +(Vwind_wc(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    U_wminus5 = U_wminus5+(Vwind_wminus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1); 
                                    U_wplus10 =U_wplus10 +(Vwind_wplus10(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    U_wplus20 = U_wplus20+ (Vwind_wplus20(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                    
                                    U_wplus5 = U_wplus5 +(Vwind_wplus5(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                                                                                           
                                    U_wminus40 = U_wminus40+(Vwind_wminus40(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                          
                                    U_wminus30 = U_wminus30 +(Vwind_wminus30(i,j))* hFacC(i,j,1)*RAC(i,j)*delR(1);                                     

                                    Vol2 = Vol2 + hFacC(i,j,1)*RAC(i,j)*delR(1);
                        end
                                    
                              
      end
end

   


U_c = U_c/Vol2;
U_cminus30 = U_cminus30/Vol2;
U_cminus40 = U_cminus40/Vol2;
U_cminus50= U_cminus50/Vol2;
U_cminus10=U_cminus10/Vol2;
U_cminus20=U_cminus20/Vol2;
U_cplus5=U_cplus5/Vol2;
U_cminus5 = U_cminus5/Vol2;
U_cplus10=U_cplus10/Vol2;
U_cplus15=U_cplus15/Vol2;
U_cplus20=U_cplus20/Vol2;
I_c = I_c/Vol;
I_cminus30 = I_cminus30/Vol;
I_cminus40 = I_cminus40/Vol;
I_cminus50= I_cminus50/Vol;
I_cminus5= I_cminus5/Vol;
I_cminus10=I_cminus10/Vol;
I_cminus20=I_cminus20/Vol;
I_cplus5=I_cplus5/Vol;
I_cplus10=I_cplus10/Vol;
I_cplus15=I_cplus15/Vol;
I_cplus20=I_cplus20/Vol;


I_w = I_w/Vol;
U_w = U_w/Vol2;
U_wminus30 = U_wminus30/Vol2;
U_wminus40 = U_wminus40/Vol2;
U_wminus50= U_wminus50/Vol2;
U_wminus10=U_wminus10/Vol2;
U_wminus20=U_wminus20/Vol2;
U_wplus5=U_wplus5/Vol2;
U_wminus5 = U_wminus5/Vol2;
U_wplus10=U_wplus10/Vol2;
U_wplus15=U_wplus15/Vol2;
U_wplus20=U_wplus20/Vol2;
I_wminus30 = I_wminus30/Vol;
I_wminus40 = I_wminus40/Vol;
I_wminus50= I_wminus50/Vol;
I_wminus5= I_wminus5/Vol;
I_wminus10=I_wminus10/Vol;
I_wminus20=I_wminus20/Vol;
I_wplus5=I_wplus5/Vol;
I_wplus10=I_wplus10/Vol;
I_wplus15=I_wplus15/Vol;
I_wplus20=I_wplus20/Vol;

figure(1);
clf;
hold on
alpha = [.5,.6,.7,.8,.9,.95,1,1.05,1.1,1.15 ,1.2];

plot(alpha(7),U_c,'-ob','markerfacecolor','blue','markersize',10);
hold on
% text(I_c,U_c,'C');


plot(alpha(1),U_cminus50,'-ob','markerfacecolor','blue','markersize',10);
% text(U_cminus50+.001,I_cminus50,'-50%');

hold on

plot(alpha(2),U_cminus40,'-ob','markerfacecolor','blue','markersize',10);
% text(U_cminus50+.001,I_cminus50,'-50%');

hold on
plot(alpha(3),U_cminus30,'-ob','markerfacecolor','blue','markersize',10);
% text(U_cminus50+.001,I_cminus50,'-50%');

hold on
plot(alpha(10),U_cplus15,'-ob','markerfacecolor','blue','markersize',10);
% text(U_cminus50+.001,I_cminus50,'-50%');

hold on
plot(alpha(4),U_cminus20,'-ob','markerfacecolor','blue','markersize',10);
% text(I_cminus20,U_cminus20,'-20%');

hold on

plot(alpha(6),U_cminus5,'-ob','markerfacecolor','blue','markersize',10);
% text(I_cminus5,U_cminus5,'-5%');

hold on

plot(alpha(5),U_cminus10,'-ob','markerfacecolor','blue','markersize',10);
% text(I_cminus10,U_cminus10,'-10%');
hold on

plot(alpha(8),U_cplus5,'-ob','markerfacecolor','blue','markersize',10);
% text(U_cminus50+.001,I_cminus50,'-50%');
hold on

plot(alpha(9),U_cplus10,'-ob','markerfacecolor','blue','markersize',10);
% text(U_cminus50+.001,I_cminus50,'-50%');



hold on

plot(alpha(11),U_cplus20,'-ob','markerfacecolor','blue','markersize',10);
% text(I_cplus20,U_cplus20,'+20%');

hold on




%%%%%% Warm Experiments
plot(alpha(7),U_w,'-or','markerfacecolor','red','markersize',10);
hold on
% text(I_w,U_w,'C');


plot(alpha(6),U_wminus5,'-or','markerfacecolor','red','markersize',10);
% text(I_wminus5,U_wminus5,'-5%');

hold on

plot(alpha(4),U_wminus20,'-or','markerfacecolor','r','markersize',10);
hold on
% text(I_wminus20,U_wminus20,'-20%');


plot(alpha(5),U_wminus10,'-or','markerfacecolor','r','markersize',10);
% text(I_wminus10,U_wminus10,'-10%');

hold on
plot(alpha(3),U_wminus30,'-or','markerfacecolor','r','markersize',10);
% text(I_wminus10,U_wminus10,'-10%');

hold on
% plot(U_wplus2,I_wplus2,'-or','markerfacecolor','r');
% hold on
% 
plot(alpha(8),U_wplus5,'-or','markerfacecolor','r','markersize',10);
% text(I_wplus5,U_wplus5,'+5%');

hold on

plot(alpha(10),U_wplus15,'-or','markerfacecolor','r','markersize',10);

hold on
plot(alpha(9),U_wplus10,'-or','markerfacecolor','r','markersize',10);
% text(I_wplus10,U_wplus10,'+10%');

hold on

plot(alpha(11),U_wplus20,'-or','markerfacecolor','r','markersize',10);
hold on
% text(I_wplus20,U_wplus20,'+20%');

plot(alpha(1),U_wminus50,'-or','markerfacecolor','r','markersize',10);
hold on
% text(I_wplus20,U_wplus20,'+20%');

hold on
plot(alpha(2),U_wminus40,'-or','markerfacecolor','r','markersize',10);
hold on


hold on
% text(I_wplus20,U_wplus20,'+20%');

hold on
% plot(U_wplus15,I_wplus15,'-or','markerfacecolor','r');
% text(U_wplus15+.001,I_wplus15,'+15%');
% 
% xlim([.025 .05]);



title('Mean Offshore $V_{wind}$ vs Wind Perturbation ($\chi$) over Ronne Polynya','Interpreter','latex','FontSize',16);
xlabel('Wind Perturbation ($\chi$) ','Interpreter','latex','FontSize',13)

ylabel('Mean Offshore V$_{Wind}$($m/s$)','Interpreter','latex','FontSize',13)



% ws_w = [U_wminus50,U_wminus40,U_wminus30,U_wminus20,U_wminus10,U_wminus5,U_w,U_wplus5,U_wplus10,U_wplus20];
% Is_w = [I_wminus50,I_wminus40,I_wminus30,I_wminus20,I_wminus10,I_wminus5,I_w,I_wplus5,I_wplus10,I_wplus20];
% 
% p = polyfit(log(Is_w),log(ws_w),1);
% z = polyval(p,log(Is_w));
% hold on
% loglog(Is_w,exp(z));
% ws_c = [U_cminus50,U_cminus40,U_cminus30,U_cminus20,U_cminus10,U_cminus5,U_c,U_cplus5,U_cplus10,U_cplus20];
% Is_c = [I_cminus50,I_cminus40,I_cminus30,I_cminus20,I_cminus10,I_cminus5,I_c,I_cplus5,I_cplus10,I_cplus20];
% 
% z = fitlm(Is_c,ws_c);







