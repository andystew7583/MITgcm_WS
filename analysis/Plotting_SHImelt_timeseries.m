%%%%%%% plot Melt Time Series

%%% Read experiment data
setExpname
loadexp;

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(end);

%%% Frequency of diagnostic output
% nIter0 =1568180;
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters_orig = round((1:nDumps)*dumpFreq/deltaT); 
dumpIters = dumpIters_orig(dumpIters_orig >= nIter0);
nDumps = length(dumpIters);
% 

deltaT_2 = 400;
nIter0_2 = 0;
nDumps_2 = round(nTimeSteps*deltaT_2/dumpFreq);
dumpIters_2 = round((1:nDumps_2)*dumpFreq/deltaT_2); 
dumpIters_2 = dumpIters_2(dumpIters_2 >= nIter0_2);
nDumps_2 = length(dumpIters_2);

deltaT_3 = 300;
nIter0_3 = 1;
nDumps_3 = round(nTimeSteps*deltaT_3/dumpFreq);
dumpIters_3 = round((1:nDumps_3)*dumpFreq/deltaT_3); 
dumpIters_3 = dumpIters_3(dumpIters_3 >= nIter0_3);


deltaT_4 = 440;
nIter0_4 = 0;
nDumps_4 = round(nTimeSteps*deltaT_4/dumpFreq);
dumpIters_4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters_4 = dumpIters_4(dumpIters_4 >= nIter0_4);

deltaT_5 = 200;
nIter0_5 = 0;
nDumps_5 = round(nTimeSteps*deltaT_5/dumpFreq);
dumpIters_5 = round((1:nDumps_5)*dumpFreq/deltaT_5); 
dumpIters_5 = dumpIters_5(dumpIters_5 >= nIter0_5);

% nDumps = 120;
% nDumps_2 = 120;
tt = zeros(1,nDumps);
tt_2 = zeros(1,nDumps);
melttot_exp1 = zeros(1,nDumps);
melttot_exp2 = zeros(1,nDumps+1);
melttot_exp3 = zeros(1,nDumps+1);


%%% To calculate time averages
melt_mean = zeros(Nx,Ny);
mean_cntr = 0;
meltlen = 0;
meltlen_3 = 0;
meltlen_2 = 0;
p = 1;
for n=1:nDumps
%    
% % if n<=19
  tt(n) =  dumpIters(n)*deltaT/86400/365;
  tt(n)
%   
  melt_exp1 = rdmdsWrapper(fullfile(exppath,'/results/SHIfwFlx'),dumpIters(n));  
  if (isempty(melt_exp1))
        break 
  end
% 
% 
% 
%   
  %%% Time-averages
  if (n > 0)
    melt_mean = melt_mean + melt_exp1;
    mean_cntr = mean_cntr + 1;
  end  
  
  %%% Integrate melt over the FRIS
melttot_exp1(n)  =0;
for i=1:Nx
    for j=1:Ny
        if XC(i,1)>-80 && XC(i,1)<-20
            
            if YC(1,j) <-74

                melttot_exp1(n)=melttot_exp1(n)+ melt_exp1(i,j)*RAC(i,j);
            end
        end
            
        
    end
 end
  
 meltlen=meltlen+1;




end

 
%  
% for n=2:nDumps+1
%   
% %   exppath_5 = '/data3/MITgcm_WS/experiments/htc5_800m_cont';       
% %   exppath_6 = '/data3/MITgcm_WS/experiments/htc5_800m';       
% 
%   tt_3(n) =  dumpIters_orig(n)*deltaT/86400/365;
%   tt_3(n)
%   exppath_2 = '/data3/MITgcm_WS/experiments/n_349';
%      if n<=24
%             melt_exp2 = rdmdsWrapper(fullfile(exppath_2,'/results/SHIfwFlx'),dumpIters_3(n));      
%              if (isempty(melt_exp2))
%                      break;
%              end
%      
%    
% 
% 
%             
%           
% 
%              
%      elseif (n>=25) 
%          
%            melt_exp2 = rdmdsWrapper(fullfile(exppath_2,'/results/SHIfwFlx'),dumpIters_5(n-5));      
%              if (isempty(melt_exp2))
%                      break;
%              end
%      end
%      
% 
%    
% % 
%      
%   %%% Integrate melt over the FRIS
%   
% % melttot_exp2(n)  =0;
% % for i=1:Nx
% %     for j=1:Ny
% %         if XC(i,1)>-80 && XC(i,1)<-20
% %             
% %             if YC(1,j) <-74
% % 
% %                  melttot_exp2(n)=melttot_exp2(n)+ melt_exp2(i,j)*RAC(i,j);
% %             end
% %         else
% %             melttot_exp2(n) = melttot_exp2(n);
% %               
% %         end
% %     end
% % end
% %  meltlen_2 = meltlen_2+1;
% % end 
% % %  
% % % % 
% deltaT_3 = 300;
% nIter0_3 = 6480;
% nDumps_3 = round(nTimeSteps*deltaT_3/dumpFreq);
% dumpIters_3 = round((1:nDumps_3)*dumpFreq/deltaT_3); 
% dumpIters_3 = dumpIters_3(dumpIters_3 >= nIter0_3);
%  for n=2:nDumps+1
% 
%          
%   tt_2(n) =  dumpIters_orig(n)*deltaT/86400/365;
%   tt_2(n)
%   exppath_3 = '/data3/MITgcm_WS/experiments/a_34_20boundary';
%      if n<=4
%             melt_exp3 = rdmdsWrapper(fullfile(exppath_3,'/results/SHIfwFlx'),dumpIters_2(n));      
%              if (isempty(melt_exp3))
%                      break;
%              end
%      elseif (n>=5)&&(n<=51 )
%      
%             
%             melt_exp3 = rdmdsWrapper(fullfile(exppath_3,'/results/SHIfwFlx'),dumpIters_3(n-2));      
% 
%      elseif n>51 
%            tt_2(n) = dumpIters_4(n)*deltaT_4/86400/365;
%            melt_exp3 = rdmdsWrapper(fullfile(exppath_3,'/results/SHIfwFlx'),dumpIters_4(n+20));      
%              if (isempty(melt_exp3))
%                      
%                  break;
%              end
%      end
%  
%   
%  melttot_exp3(n) = 0;
% for i=1:Nx
%     for j=1:Ny
%             if XC(i,1)>-80 && XC(i,1)<-20
%                 
%                 if YC(1,j) <-74       
%                    
%             
%                     melttot_exp3(n) = melttot_exp3(n) + melt_exp3(i,j)*RAC(i,j);
%                 end
%             else
%                     melttot_exp3(n) = melttot_exp3(n);
%                        
%             end
%             
%      end
% end
% 
% meltlen_3=meltlen_3+1;
%   
%  
% end
%  
% 
melttot_exp1 = (melttot_exp1*86400*365*1e-12);  % convrt to GT/year;  
% 
% melttot_exp2 = (melttot_exp2*86400*365*1e-12); 
% % % 
melttot_exp3 = (melttot_exp3*86400*365*1e-12); 
% % %
% % 
% % melttot_exp3(melttot_exp3==0)=NaN;
% % % melttot_exp2(melttot_exp2==0)=NaN;
% % 
% % melt_mean = melt_mean/mean_cntr;
% % 
figure(1);
clf;
hold on
c1=plot(tt(1:meltlen),-melttot_exp1(1:meltlen),'b','linewidth',1);
hold on 
% n=plot(tt_3(1:meltlen_2),-melttot_exp2(1:meltlen_2),'g');
hold on
% w=plot(tt_2(1:meltlen_3),-melttot_exp3(1:meltlen_3),'r');

% 
% 
% % 
xlabel('t (years)','interpreter','latex');
ylabel('Integrated Shelf Ice Melt (Gt/yr)','interpreter','latex');
title('Integrated FRIS Melt', 'interpreter','latex');
% % 
% % hold on
% % 

%%%shade up to 1 years
h1 = line([0 0],[1 2000]);
h2 = line([2 2],[1 2000]);
c = patch([0 2 2 0],[1 1 2000 2000],[.5 .5 .5]);

c.FaceAlpha=.3;
axis([0 25 0 2000]);

z = 155*ones(40,1);
k3 = plot(0:39,z,'k:','linewidth',2);
hold on
% % 
legend([c1,k3],{'34.45psu','Observed Melt'},'location','northeast','interpreter','latex','fontsize',13)
% 



%          
%          
%          