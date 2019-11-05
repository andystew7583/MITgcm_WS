%%% This needs to be set to ensure we are using the correct output
%%% frequency
loadexp;
diagfreq = diag_frequency(end);

%%% Frequency of diagnostic output

nIter0 = 0;
%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters_orig = round((1:nDumps)*dumpFreq/deltaT); 
dumpIters = dumpIters_orig(dumpIters_orig >= nIter0);


%%% Load velocity

loadexp;
%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
deltaT = 440;
% nIter0 = 0;
nDumps = round(nTimeSteps*(deltaT/dumpFreq));
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);
 
deltaT_2 = 400;
nIter0_2 = 0;
nDumps_2 = round(nTimeSteps*deltaT_2/dumpFreq);
dumpIters_2 = round((1:nDumps_2)*dumpFreq/deltaT_2); 
dumpIters_2 = dumpIters_2(dumpIters_2 >= nIter0_2);
nDumps_2 = length(dumpIters_2);



deltaT_3 = 300;
nIter0_3 = 6480;
nDumps_3 = round(nTimeSteps*deltaT_3/dumpFreq);
dumpIters_3 = round((1:nDumps_3)*dumpFreq/deltaT_3); 
dumpIters_3 = dumpIters_3(dumpIters_3 >= nIter0_3);

deltaT_4 = 200;
nIter0_4 = 0;
nDumps_4 = round(nTimeSteps*deltaT_4/dumpFreq);
dumpIters_4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters_4 = dumpIters_4(dumpIters_4 >= nIter0_4);



deltaT_5 = 200;
nIter0_5 = 0;
nDumps_5 = round(nTimeSteps*deltaT_5/dumpFreq);
dumpIters_5 = round((1:nDumps_5)*dumpFreq/deltaT_5); 
dumpIters_5 = dumpIters_5(dumpIters_5 >= nIter0_5);

%%% Mesh grids for plotting
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
kn = ones(Nx,Ny);
kp= ones(Nx,Ny);
wn = 0.5*ones(Nx,Ny);
wp = 0.5*ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    if (sum(hFacC(i,j,:),3)==0)
      continue;
    end
    zmid = 0.5 * (SHELFICEtopo(i,j) + bathy(i,j));
    kmid = max(find(squeeze(zz)>zmid));
    if (isempty(kmid))
      continue;
    end
    kp(i,j) = kmid;
    kn(i,j) = kp(i,j) + 1;
    wp(i,j) = (zmid-zz(kn(i,j))) / (zz(kp(i,j))-zz(kn(i,j)));
    wn(i,j) = 1 - wp(i,j);
  end
end


tt = zeros(1,28*12);
tt2=zeros(1,28*12);
thetatot_exp1 = zeros(1,28*12);
thetatot_exp2 = zeros(1,28*12);
thetatot_exp3 = zeros(1,28*12);
salttot_exp1 = zeros(1,28*12);
salttot_exp2 = zeros(1,28*12);

melttot_exp1=zeros(1,28*12);
melttot_exp2=zeros(1,28*12);
melttot_exp3=zeros(1,28*12);


theta_cntr = 0;
thetalen = 0;
% thetalen2 = 0;
% thetalen3 = 0;
% for n=1:nDumps

start_dump = 1;
end_dump = 28*12;



for n=1:end_dump
  tt(n) =  dumpIters_3(n)*deltaT_3/86400/365;
  tt(n)
  
  %%% Attempt to load melt ave per month
  theta_exp1 = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));  
  salt_exp1 = rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n));  
  melt_exp1=rdmdsWrapper(fullfile(exppath,'/results/SHIfwFlx'),dumpIters(n));  
  if (isempty(melt_exp1))
        break     
  end
        

  for i=1:Nx
    for j=1:Ny
     if hFacC(i,j,kmax(i,j))>0
        if XC(i,j)>-80 && XC(i,j)<-20
            
            if YC(i,j) <-75
                if SHELFICEtopo(i,j) <0
                  if melt_exp1(i,j)<0
                    thetatot_exp1(n) = thetatot_exp1(n) + (theta_exp1(i,j,kmax(i,j))*RAC(i,j)*melt_exp1(i,j));
                    salttot_exp1(n) = salttot_exp1(n) + (salt_exp1(i,j,kmax(i,j))*RAC(i,j)*melt_exp1(i,j));
                    melttot_exp1(n)=melttot_exp1(n) + (melt_exp1(i,j)*RAC(i,j));
                  end
                end
            end
        end
     end
    end

  
  
  
  end
end

nIter0_3 = 6480;
nDumps_3 = round(nTimeSteps*deltaT_3/dumpFreq);
dumpIters_3 = round((1:nDumps_3)*dumpFreq/deltaT_3); 
dumpIters_3 = dumpIters_3(dumpIters_3 >= nIter0_3);
 for n=2:nDumps+1

         
  tt_2(n) =  dumpIters_orig(n)*deltaT/86400/365;
  tt_2(n);
  exppath_3 = '/data3/MITgcm_WS/experiments/a_34_20boundary';
     if n<=4
           tt_2(n) =  dumpIters_orig(n)*deltaT/86400/365;

            theta_exp2 = rdmdsWrapper(fullfile(exppath_3,'/results/THETA'),dumpIters_2(n));     
            salt_exp2 = rdmdsWrapper(fullfile(exppath_3,'/results/SALT'),dumpIters_2(n));  
            melt_exp2=rdmdsWrapper(fullfile(exppath_3,'/results/SHIfwFlx'),dumpIters_2(n));  

     elseif (n>=5)&&(n<=51 )

            theta_exp2 = rdmdsWrapper(fullfile(exppath_3,'/results/THETA'),dumpIters_3(n-2));      
            salt_exp2 = rdmdsWrapper(fullfile(exppath_3,'/results/SALT'),dumpIters_3(n-2));      
            melt_exp2 = rdmdsWrapper(fullfile(exppath_3,'/results/SHIfwFlx'),dumpIters_3(n-2));      

     elseif n>51 
           theta_exp2 = rdmdsWrapper(fullfile(exppath_3,'/results/THETA'),dumpIters_orig(n+20));    
           salt_exp2 = rdmdsWrapper(fullfile(exppath_3,'/results/SALT'),dumpIters_orig(n+20));      
           melt_exp2 = rdmdsWrapper(fullfile(exppath_3,'/results/SHIfwFlx'),dumpIters_orig(n+20));           
           
           
             if (isempty(theta_exp2))
                     
                 break;
             end
     end
 
  
  for i=1:Nx
    for j=1:Ny
     if hFacC(i,j,kmax(i,j))>0
        if XC(i,j)>-80 && XC(i,j)<-20
            
            if YC(i,j) <-75
                if SHELFICEtopo(i,j) <0
                  if melt_exp2(i,j)<0
                    thetatot_exp2(n) = thetatot_exp2(n) + (theta_exp2(i,j,kmax(i,j))*RAC(i,j)*melt_exp2(i,j));
                    salttot_exp2(n) = salttot_exp2(n) + (salt_exp2(i,j,kmax(i,j))*RAC(i,j)*melt_exp2(i,j));
                    
                    melttot_exp2(n)=melttot_exp2(n) + (melt_exp2(i,j)*RAC(i,j));
                  end
                end
            end
        end
     end
    end

  
  
  
  end
  
 
 end



% 
% for n=1:end_dump
% 
%   if n<=24
%             theta_exp3 = rdmdsWrapper(fullfile(exppath_3,'/results/THETA'),dumpIters_3(n));      
%             melt_exp3=rdmdsWrapper(fullfile(exppath_3,'/results/SHIfwFlx'),dumpIters_3(n));  
%             
%              if (isempty(theta_exp3))
%                      break;
%              end
%      
%  
%                  
%   elseif n>=25
%          
%            theta_exp3 = rdmdsWrapper(fullfile(exppath_3,'/results/THETA'),dumpIters_5(n-5)); 
%            melt_exp3=rdmdsWrapper(fullfile(exppath_3,'/results/SHIfwFlx'),dumpIters_5(n-5));  
%           
%              if (isempty(theta_exp3))
%                      break;
%              end
%              if (isempty(melt_exp3))
%                         break;
%              end
%   end
%   
%   for i=1:Nx
%     for j=1:Ny
%      if hFacC(i,j,kmax(i,j))>0
% 
%         if XC(i,j)>-80 && XC(i,j)<-20
%             
%             if YC(i,j) <-75
%                 if SHELFICEtopo(i,j) <0
%                    if melt_exp3(i,j)<0
%                     thetatot_exp3(n) = thetatot_exp3(n) + (theta_exp3(i,j,kmax(i,j))*RAC(i,j)*melt_exp3(i,j));
%                     melttot_exp3(n)=melttot_exp3(n) + (melt_exp3(i,j)*RAC(i,j));
% 
%                    end
%                 end
%             end
%         end
%      end
%     end
% 
%   
%   
%   
%   end
% end
%   
                   

thetatot_exp1 = thetatot_exp1./melttot_exp1;
% theta_mean_1 = thetatot_exp1/thetalen;
thetatot_exp2 = thetatot_exp2./melttot_exp2;
salttot_exp1 = salttot_exp1./melttot_exp1;
salttot_exp2 = salttot_exp2./melttot_exp2;
