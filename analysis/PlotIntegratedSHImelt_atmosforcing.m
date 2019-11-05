%%%%%%% plot integrated shelf ice melt for array of atmos perturbation
%%%%%%% experiments

%%% Read experiment data

setExpname
loadexp
diagfreq = diag_frequency(end);

%%% Frequency of diagnostic output


dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters_orig = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters_orig(dumpIters_orig >= nIter0);
nDumps = length(dumpIters);

%%% for experiments where timestep, nIter0 is different

deltaT_2 = 400;
nDumps_2 = round(nTimeSteps*deltaT_2/dumpFreq);
nIter0_2 = 11;
dumpIters_2_orig = round((1:nDumps_2)*dumpFreq/deltaT_2);
dumpIters_2 = dumpIters_2_orig(dumpIters_2_orig >= nIter0_2);
Ny_2 = 224;
nDumps_2 = length(dumpIters_2);

deltaT_3 = 300;
nDumps_3 = round(nTimeSteps*deltaT_3/dumpFreq);
nIter0_3 = 1;
dumpIters_3 = round((1:nDumps_3)*dumpFreq/deltaT_3);
dumpIters_3 = dumpIters_3(dumpIters_3 >= nIter0_3);
% 
% nDumps_3 = length(dumpIters_3);
% 
deltaT_4 = 200;
nDumps_4 = round(nTimeSteps*deltaT_4/dumpFreq);
nIter0_4 = 1;
dumpIters_4 = round((1:nDumps_4)*dumpFreq/deltaT_4);
dumpIters_4 = dumpIters_4(dumpIters_4 >= nIter0_4);

nDumps_4 = length(dumpIters_4);


exppath1 =  '/data3/MITgcm_WS/experiments/z_control_zonalminus50';
exppath2 = '/data3/MITgcm_WS/experiments/z_minus20_control';
exppath7 = '/data3/MITgcm_WS/experiments/t_warm_tempminus5_h2';
exppath4 = '/data3/MITgcm_WS/experiments/s_hist8_nozonalwind';
exppath5 = '/data3/MITgcm_WS/experiments/s_hist6_warm_coldminus30';
exppath6 = '/data3/MITgcm_WS/experiments/s_hist4_warm_meridplus20';
exppath3 = '/data3/MITgcm_WS/experiments/s_hist38_c_minus50';
exppath8 = '/data3/MITgcm_WS/experiments/a_3445_20boundary';
exppath9 = '/data3/MITgcm_WS/experiments/s_hist11_control_30warmer';
exppath10 = '/data3/MITgcm_WS/experiments/z_warm_zonalplus20';

tt = zeros(1,nDumps);
melttot_exp1 = zeros(1,nDumps);
melttot_exp2 = zeros(1,nDumps_2);
melttot_exp3 = zeros(1,nDumps_3);
melttot_exp4 = zeros(1,nDumps_2);
melttot_exp5 = zeros(1,nDumps_2);
melttot_exp6 = zeros(1,nDumps_2);
melttot_exp7 = zeros(1,nDumps_2);
melttot_exp9 = zeros(1,nDumps_2);
melttot_exp10 = zeros(1,nDumps_2);

%%% To calculate time averages
melt_mean = zeros(Nx,Ny);
mean_cntr = 0;
meltlen=0
nDumps = 12*35;
for n=1:nDumps
  tt(n) =  dumpIters(n)*deltaT/86400/365;
  tt(n)
  
  %%% Attempt to load melt ave per month

%   exp1 = rdmdsWrapper(fullfile(exppath1,'/results/SHIfwFlx'),dumpIters_3(n));      
%     
%   if (isempty(exp1))
%         exp1=zeros(Nx,Ny);
%   end
%   
%   exp2 = rdmdsWrapper(fullfile(exppath2,'/results/SHIfwFlx'),dumpIters_3(n));      
%   if (isempty(exp2))
%         exp2=zeros(Nx,Ny);
%   end
  exp3 = rdmdsWrapper(fullfile(exppath3,'/results/SHIfwFlx'),dumpIters_3(n));      
  if (isempty(exp3))
        exp3=zeros(Nx,Ny);
  end
  
  exp4 = rdmdsWrapper(fullfile(exppath4,'/results/SHIfwFlx'),dumpIters_2(n));      
  if (isempty(exp4))
        exp4=zeros(Nx,Ny);
  end
  
  exp5 = rdmdsWrapper(fullfile(exppath5,'/results/SHIfwFlx'),dumpIters_2(n));      
  if (isempty(exp5))
        exp5=zeros(Nx,Ny);
  end
  
  exp6 = rdmdsWrapper(fullfile(exppath6,'/results/SHIfwFlx'),dumpIters_2(n));      
  if (isempty(exp6))
        exp6=zeros(Nx,Ny);
  end
%   exp7 = rdmdsWrapper(fullfile(exppath7,'/results/SHIfwFlx'),dumpIters_3(n));      
%   if (isempty(exp7))
%         exp7=zeros(Nx,Ny);
%   end
%   
%   exp9 = rdmdsWrapper(fullfile(exppath9,'/results/SHIfwFlx'),dumpIters_2(n));      
%   if (isempty(exp9))
%         exp9=zeros(Nx,Ny);
%   end
%   
%   exp10 = rdmdsWrapper(fullfile(exppath10,'/results/SHIfwFlx'),dumpIters_4(n));      
%   if (isempty(exp10))
%         exp10=zeros(Nx,Ny);
%   end
  %%% Time-averages
  
  
  %%% Integrate melt over the FRIS
%   melttot_exp1(n) = 0;
%   melttot_exp2(n)=0;
  melttot_exp3(n)=0;
  melttot_exp4(n)=0;
  melttot_exp5(n)=0;
  melttot_exp6(n)=0;
%   melttot_exp7(n)=0;
%   melttot_exp9(n)=0;
%   melttot_exp10(n)=0;

  for i=1:Nx
    for j=1:Ny
        if XC(i,1)>-80 && XC(i,1)<-20
                        if YC(1,j) <-74  
%                             melttot_exp1(n) = melttot_exp1(n) + exp1(i,j)*RAC(i,j);
%                             melttot_exp2(n) = melttot_exp2(n) + exp2(i,j)*RAC(i,j);
                            melttot_exp3(n) = melttot_exp3(n) + exp3(i,j)*RAC(i,j);
                            melttot_exp4(n) = melttot_exp4(n) + exp4(i,j)*RAC(i,j);
                            melttot_exp5(n) = melttot_exp5(n) + exp5(i,j)*RAC(i,j);
                            melttot_exp6(n) = melttot_exp6(n) + exp6(i,j)*RAC(i,j);
%                             melttot_exp7(n) = melttot_exp7(n) + exp7(i,j)*RAC(i,j);
%                             melttot_exp10(n) = melttot_exp10(n) + exp10(i,j)*RAC(i,j);
%                             melttot_exp9(n) = melttot_exp9(n) + exp9(i,j)*RAC(i,j);
% 

                        end
        
        end
    end
  end
  meltlen = meltlen + 1;

  
  
 end

 
 
% melttot_exp1(melttot_exp1==0)=NaN;
% 
% melttot_exp1 = (melttot_exp1*86400*365*1e-12);  % convrt to GT/year;  
% 
% 
% melttot_exp2(melttot_exp2==0)=NaN;
% 
% melttot_exp2 = (melttot_exp2*86400*365*1e-12); 

melttot_exp3(melttot_exp3==0)=NaN;

melttot_exp3 = (melttot_exp3*86400*365*1e-12); 
% 

melttot_exp4(melttot_exp4==0)=NaN;

melttot_exp4 = (melttot_exp4*86400*365*1e-12); 

melttot_exp5(melttot_exp5==0)=NaN;

melttot_exp5 = (melttot_exp5*86400*365*1e-12); 


melttot_exp6(melttot_exp6==0)=NaN;

melttot_exp6 = (melttot_exp6*86400*365*1e-12); 

% melttot_exp7(melttot_exp7==0)=NaN;
% 
% melttot_exp7 = (melttot_exp7*86400*365*1e-12); 
% 
% melttot_exp9(melttot_exp9==0)=NaN;
% 
% melttot_exp9 = (melttot_exp9*86400*365*1e-12); 
% melttot_exp10(melttot_exp10==0)=NaN;
% 
% melttot_exp10 = (melttot_exp10*86400*365*1e-12); 


melt_mean = melt_mean/mean_cntr;

figure(1);
clf
c =(jet(100));
set(gca,'FontSize',16);
% t1 = plot(tt(1:meltlen),-melttot_exp1(1:meltlen),'color',c(5,:),'Linewidth',1.5);
% hold on 
% t2 = plot(tt(1:meltlen),-melttot_exp2(1:meltlen),'color',c(20,:),'Linewidth',1.5);
hold on
t3 = plot(tt(1:meltlen),-melttot_exp3(1:meltlen),'color',c(10,:),'Linewidth',.5);
hold on
t4 = plot(tt(1:meltlen),-melttot_exp4(1:meltlen),'color',c(30,:),'Linewidth',.5);
hold on
t5 = plot(tt(1:meltlen),-melttot_exp5(1:meltlen),'color',c(75,:),'Linewidth',.5);
hold on
t6 = plot(tt(1:meltlen),-melttot_exp6(1:meltlen),'color',c(90,:),'Linewidth',.5);
hold on
% t7 = plot(tt(1:meltlen),-melttot_exp7(1:meltlen),'color',c(90,:),'Linewidth',1.5);
% hold on
% % t9 = plot(tt(1:meltlen),-melttot_exp9(1:meltlen),'color',c(10,:),'Linewidth',1.5);
% hold on
% % t10 = plot(tt(1:meltlen),-melttot_exp10(1:meltlen),'color',c(110,:),'Linewidth',1.5);


axis([0 25 0 950]);

xlabel('t (years)','interpreter','latex');
ylabel('Integrated Shelf Ice Melt (Gt/yr)','interpreter','latex');
title('Integrated FRIS Melt (Gt)', 'interpreter','latex','fontsize',16);
% h1 = line([0 0],[1 2000]);
% h2 = line([2 2],[1 2000]);
% c = patch([0 2 2 0],[1 1 2000 2000],[.5 .5 .5]);
% 
% c.FaceAlpha=.3;
% axis([0 25 0 1200]);

z = 124*ones(40,1);
k3 = plot(0:39,z,'k:','linewidth',2);
hold on
l = legend([t4,t3,t5,t6,k3],{'REF, No Zonal Wind','REF, Merid -50$\%$','FRESH, Temp -30$^\circ$C','FRESH, Merid +20$\%$','Observed'},'interpreter','latex','location','northwest','fontsize',12);

