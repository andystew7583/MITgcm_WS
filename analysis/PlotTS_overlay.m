%%%%%% Plot T/S 

% Script to plot Temperature and salinity

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


%%%% Load Salinity and Temperature

% 
% tmin = 12*86400*360;
% tmax = 19*86400*360;
% %
% %
% thetawc = readIters(exppath1w,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
% saltwc = readIters(exppath1w,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
% tmin = 10*86400*360;
% tmax = 19*86400*360;
% thetawcplus10 = readIters(exppath5w,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
% saltwcplus10 = readIters(exppath5w,'SALT',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
% 
% tmin = 28*86400*360;
% tmax = 37*86400*360;
% thetawcplus20 = readIters(exppath2w,'THETA',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
% saltwcplus20 = readIters(exppath2w,'SALT',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
% tmin = 19*86400*360;
% tmax = 27*86400*360;
% thetawcplus15 = readIters(exppath10w,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
% saltwcplus15 = readIters(exppath10w,'SALT',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
% 
% tmin = 10*86400*360;
% tmax = 19*86400*360;
% thetawcplus5 = readIters(exppath11w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% saltwcplus5 = readIters(exppath11w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% thetawcminus5 = readIters(exppath9w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% saltwcminus5 = readIters(exppath9w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% thetawcminus20 = readIters(exppath8w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% saltwcminus20 = readIters(exppath8w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% thetawcminus30 = readIters(exppath12w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% saltwcminus30 = readIters(exppath12w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% thetawcminus10 = readIters(exppath7w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% saltwcminus10 = readIters(exppath7w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% 
% tmin = 9*86400*360;
% tmax = 18*86400*360;
% thetawcminus50 = readIters(exppath4w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% saltwcminus50 = readIters(exppath4w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% 
% tmin = 9*86400*360;
% tmax = 18*86400*360;
% thetawcminus40 = readIters(exppath13w,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
% saltwcminus40 = readIters(exppath13w,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);


tmin = 9*86400*360;
tmax = 18*86400*360;
thetawc = readIters(exppath5,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
saltwc = readIters(exppath5,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);



thetawcplus10 = readIters(exppath2,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
saltwcplus10 = readIters(exppath2,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

thetawcplus20 = readIters(exppath3,'THETA',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);
saltwcplus20 = readIters(exppath3,'SALT',dumpIters_4,deltaT_4,tmin,tmax,Nx,Ny,Nr);


thetawcplus15 = readIters(exppath13,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
saltwcplus15 = readIters(exppath13,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);


thetawcplus5 = readIters(exppath4,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
saltwcplus5= readIters(exppath4,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

thetawcminus5 = readIters(exppath7,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
saltwcminus5 = readIters(exppath7,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

thetawcminus20 = readIters(exppath9,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
saltwcminus20 = readIters(exppath9,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

tmin = 9*86400*360;
tmax = 16*86400*360;
thetawcminus30 = readIters(exppath11,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
saltwcminus30 = readIters(exppath11,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

tmin = 9*86400*360;
tmax = 18*86400*360;
thetawcminus10 = readIters(exppath10,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
saltwcminus10 = readIters(exppath10,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);


tmin = 21*86400*360;
tmax = 29*86400*360;
thetawcminus50 = readIters(exppath1,'THETA',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);
saltwcminus50 = readIters(exppath1,'SALT',dumpIters_3,deltaT_3,tmin,tmax,Nx,Ny,Nr);

tmin = 21*86400*360;
tmax = 30*86400*360;
thetawcminus40 = readIters(exppath12,'THETA',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
saltwcminus40 = readIters(exppath12,'SALT',dumpIters_2,deltaT2,tmin,tmax,Nx,Ny,Nr);

                thetawc_n = zeros(Nx,Ny,Nr);
                thetawcplus20_n=zeros(Nx,Ny,Nr);
                thetawcplus10_n=zeros(Nx,Ny,Nr);
                thetawcplus5_n=zeros(Nx,Ny,Nr);
                thetawcminus5_n=zeros(Nx,Ny,Nr);
                thetawcplus15_n=zeros(Nx,Ny,Nr);
                thetawcminus10_n=zeros(Nx,Ny,Nr);
                thetawcminus20_n=zeros(Nx,Ny,Nr);
                thetawcminus30_n=zeros(Nx,Ny,Nr);
                thetawcminus40_n=zeros(Nx,Ny,Nr);
                thetawcminus50_n=zeros(Nx,Ny,Nr);
                
                saltwc_n= zeros(Nx,Ny,Nr);
                saltwcplus20_n=zeros(Nx,Ny,Nr);
                saltwcplus10_n=zeros(Nx,Ny,Nr);
                saltwcplus5_n=zeros(Nx,Ny,Nr);
                saltwcminus5_n=zeros(Nx,Ny,Nr);
                saltwcplus15_n=zeros(Nx,Ny,Nr);
                saltwcminus10_n=zeros(Nx,Ny,Nr);
                saltwcminus20_n=zeros(Nx,Ny,Nr);
                saltwcminus30_n=zeros(Nx,Ny,Nr);
                saltwcminus40_n=zeros(Nx,Ny,Nr);
                saltwcminus50_n=zeros(Nx,Ny,Nr);

%%% define where the cavity is
for i = 1:Nx
    for j = 1:Ny
       for k = 1:Nr
        if XC(i,j)>-80 && XC(i,j)<-30 && YC(i,j)<-75
            if SHELFICEtopo(i,j)<0
                
                thetawc_n(i,j,k) = thetawc(i,j,k);
                thetawcplus20_n(i,j,k)=thetawcplus20(i,j,k);
                thetawcplus10_n(i,j,k)=thetawcplus10(i,j,k);
                thetawcplus5_n(i,j,k)=thetawcplus5(i,j,k);
                thetawcminus5_n(i,j,k)=thetawcminus5(i,j,k);
                thetawcplus15_n(i,j,k)=thetawcplus15(i,j,k);
                thetawcminus10_n(i,j,k)=thetawcminus10(i,j,k);
                thetawcminus20_n(i,j,k)=thetawcminus20(i,j,k);
                thetawcminus30_n(i,j,k)=thetawcminus30(i,j,k);
                thetawcminus40_n(i,j,k)=thetawcminus40(i,j,k);
                thetawcminus50_n(i,j,k)=thetawcminus50(i,j,k);
                
                saltwc_n(i,j,k) = saltwc(i,j,k);
                saltwcplus20_n(i,j,k)=saltwcplus20(i,j,k);
                saltwcplus10_n(i,j,k)=saltwcplus10(i,j,k);
                saltwcplus5_n(i,j,k)=saltwcplus5(i,j,k);
                saltwcminus5_n(i,j,k)=saltwcminus5(i,j,k);
                saltwcplus15_n(i,j,k)=saltwcplus15(i,j,k);
                saltwcminus10_n(i,j,k)=saltwcminus10(i,j,k);
                saltwcminus20_n(i,j,k)=saltwcminus20(i,j,k);
                saltwcminus30_n(i,j,k)=saltwcminus30(i,j,k);
                saltwcminus40_n(i,j,k)=saltwcminus40(i,j,k);
                saltwcminus50_n(i,j,k)=saltwcminus50(i,j,k);
              
            end
        end
       end
    end
end


%%%%%%%%%%%%%%%%%




Nx_n = size(thetawc_n,1);
Ny_n = size(thetawc_n,2);
    
thetawc_n = reshape(thetawc_n,[1],Nx_n*(Ny_n)*(Nr));
saltwc_n = reshape(saltwc_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcminus5_n = reshape(thetawcminus5_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcminus5_n = reshape(saltwcminus5_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcminus50_n = reshape(thetawcminus50_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcminus50_n = reshape(saltwcminus50_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcminus30_n = reshape(thetawcminus30_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcminus30_n = reshape(saltwcminus30_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcminus40_n = reshape(thetawcminus40_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcminus40_n = reshape(saltwcminus40_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcminus20_n = reshape(thetawcminus20_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcminus20_n = reshape(saltwcminus20_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcplus5_n = reshape(thetawcplus5_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcplus5_n = reshape(saltwcplus5_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcplus10_n = reshape(thetawcplus10_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcplus10_n = reshape(saltwcplus10_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcplus15_n = reshape(thetawcplus15_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcplus15_n = reshape(saltwcplus15_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcplus20_n = reshape(thetawcplus20_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcplus20_n = reshape(saltwcplus20_n,[1],(Nx_n)*(Ny_n)*(Nr));

thetawcminus10_n = reshape(thetawcminus10_n,[1],Nx_n*(Ny_n)*(Nr));
saltwcminus10_n = reshape(saltwcminus10_n,[1],(Nx_n)*(Ny_n)*(Nr));

saltwcminus50_n=saltwcminus50_n(~saltwcminus50_n==0);
thetawcminus50_n=thetawcminus50_n(~thetawcminus50_n==0);
saltwcminus40_n=saltwcminus40_n(~saltwcminus40_n==0);
thetawcminus40_n=thetawcminus40_n(~thetawcminus40_n==0);
saltwcminus30_n=saltwcminus30_n(~saltwcminus30_n==0);
thetawcminus30_n=thetawcminus30_n(~thetawcminus30_n==0);
saltwcplus20_n=saltwcplus20_n(~saltwcplus20_n==0);
thetawcplus20_n=thetawcplus20_n(~thetawcplus20_n==0);
saltwcminus10_n=saltwcminus10_n(~saltwcminus10_n==0);
thetawcminus10_n=thetawcminus10_n(~thetawcminus10_n==0);
saltwc_n=saltwc_n(~saltwc_n==0);
thetawc_n=thetawc_n(~thetawc_n==0);
saltwcplus10_n=saltwcplus10_n(~saltwcplus10_n==0);
thetawcplus10_n=thetawcplus10_n(~thetawcplus10_n==0);
saltwcplus15_n=saltwcplus15_n(~saltwcplus15_n==0);
thetawcplus15_n=thetawcplus15_n(~thetawcplus15_n==0);

figure(1)
clf
cmap = lines((11));
set(gca,'FontSize',15);
hold on

s1=scatter(saltwcminus50_n(1:5:end),thetawcminus50_n(1:5:end),5,cmap(10,:),'filled');
hold on
s2=scatter(saltwcminus40_n(1:5:end),thetawcminus40_n(1:5:end),5,cmap(9,:),'filled');
hold on
s3=scatter(saltwcminus30_n(1:5:end),thetawcminus30_n(1:5:end),5,cmap(11,:),'filled');
hold on
% s4=scatter(saltwcminus20,thetawcminus20,10,cmap(2,:),'filled');
hold on;
s5=scatter(saltwcminus10_n(1:5:end),thetawcminus10_n(1:5:end),5,cmap(3,:),'filled');
hold on
% s6=scatter(saltwcminus5,thetawcminus5,10,cmap(4,:),'filled');
hold on
s7=scatter(saltwc_n(1:5:end),thetawc_n(1:5:end),5,cmap(1,:),'filled');
hold on;
% s8=scatter(saltwcplus5,thetawcplus5,10,cmap(7,:),'filled');
hold on
s9=scatter(saltwcplus10_n( 1:5:end),thetawcplus10_n(1:5:end),5,cmap(6,:),'filled');
hold on
s10=scatter(saltwcplus15_n(1:5:end),thetawcplus15_n(1:5:end),5,cmap(8,:),'filled');
hold on
s11=scatter(saltwcplus20_n(1:5:end),thetawcplus20_n(1:5:end),5,cmap(5,:),'filled');
hold on

hold on
title('Temperature/Salinity in FRIS Cavity, ``Warm" FRIS Initialization Experiments', 'interpreter','latex','Fontsize',17);
xlabel('Salinity (psu)','interpreter','latex');
ylabel('Temperature ($^\circ$C)','interpreter','latex');
xlim([33.7 35]);
ylim([-3 .5]);



thetawc_n(thetawc_n==0)=[];
saltwc_n(saltwc_n==0)=[];
%%% Get ranges of T/S
ss_max = max(max(saltwc_n))+2;
pt_max = max(max(thetawc_n)+2);
pt_min = min(min(thetawc_n)-1);
ss_min = min(min(saltwc_n))-1;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/100;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/100;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% Calculate potential density
pd = densmdjwf(SS_grid,PT_grid,500) - 1000;
hold on
[C,f] = contour(SS_grid,PT_grid,pd,20.0:.1:40.0,'EdgeColor','k');
clabel(C,f)
set(gca,'fontsize',17);

load starred
z = -1.8*ones(6,1);
k3 = plot(33:38,z,'k:','linewidth',2);
hold on

h = legend([s1,s2,s3,s5,s7,s9,s10,s11,k3],{'-50\%','-40\%','-30\%','-10\%','Control','+10\%','+15\%','+20\%','Freezing Temperature($^\circ$C)'},'interpreter','latex');
set(h,'FontSize',13);
