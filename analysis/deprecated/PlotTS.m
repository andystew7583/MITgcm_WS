%%%%%% Plot T/S 

% Script to plot Temperature and salinity

%%% Read experiment data
loadexp;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
% deltaT = 440;
nIter0 = 0;
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%%% Load Salinity and Temperature

% tmin = 10;
% tmax = 19;
tmin = 9;
tmax = 18;
%%%number of months pertaining to these timeframes
d_s = tmin*12;
d_f = tmax*12;

theta_tot = NaN(Nx,Ny,Nr,d_f-d_s);
salt_tot = NaN(Nx,Ny,Nr,d_f-d_s);

for n=1:(d_f-d_s)
  
  SALT = rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n+d_s));
  salt_tot(:,:,:,n) = (SALT); %%%convert to kg/m2/s
  TEMP= rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n+d_s));
  theta_tot(:,:,:,n) = (TEMP); %%%convert to kg/m2/s
  %tt(n) is the day of the run
  
end

SALT_complete= squeeze(nanmean(salt_tot,4));
THETA_complete = squeeze(nanmean(theta_tot,4));
Nyears=9;

%%%finding seasonal averages
SALT_season = zeros(size(salt_tot,1),size(salt_tot,2),size(salt_tot,3),Nyears,4);
Temp_season = zeros(size(salt_tot,1),size(salt_tot,2),size(salt_tot,3),Nyears,4);

for n = 1:Nyears
  SALT_season(:,:,:,n,1) = mean(salt_tot(:,:,:,(n-1)*12+[12 1 2]),4);
  SALT_season(:,:,:,n,2) = mean(salt_tot(:,:,:,(n-1)*12+(3:5)),4);
  SALT_season(:,:,:,n,3) = mean(salt_tot(:,:,:,(n-1)*12+(6:8)),4);
  SALT_season(:,:,:,n,4) = mean(salt_tot(:,:,:,(n-1)*12+(9:11)),4);
  Temp_season(:,:,:,n,1) = mean(theta_tot(:,:,:,(n-1)*12+[12 1 2]),4);
  Temp_season(:,:,:,n,2) = mean(theta_tot(:,:,:,(n-1)*12+(3:5)),4);
  Temp_season(:,:,:,n,3) = mean(theta_tot(:,:,:,(n-1)*12+(6:8)),4);
  Temp_season(:,:,:,n,4) = mean(theta_tot(:,:,:,(n-1)*12+(9:11)),4);
end


SALT_seasonavg = squeeze(mean(SALT_season,4));
Temp_seasonavg = squeeze(mean(Temp_season,4));
% Temp_JJA = squeeze(mean(Temp_seasonavg,4));
% Salt_JJA = squeeze(mean(SALT_seasonavg,4));
Temp_JJA = (Temp_seasonavg(:,:,:,3));
Salt_JJA = (SALT_seasonavg(:,:,:,3));
  
%%%finding volume of all cells in grid for normalization purposes      
volume = zeros(Nx,Ny,Nr);
nv = zeros(Nx,Ny,Nr);
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nr
            
            volume(i,j,k) = ((RAC(i,j)*DRF(k)*hFacC(i,j,k)));
            if volume(i,j,k)>0
                nv(i,j,k)=log2(volume(i,j,k));
            end
        end
    end
    
end


%%%% index the FRIS area (-80 to -20) W, 74S -> !
new_vol=zeros(Nx,Ny,Nr);
Salt_JJA_new=zeros(Nx,Ny,Nr);
Temp_JJA_new=zeros(Nx,Ny,Nr);
  for i=1:Nx
    for j=1:Ny
        for k = 1:Nr
%              if (XC(i,1)>-80) && (XC(i,1))<-20
%                  if YC(i,j)<-75
%                     if SHELFICEtopo(i,j)<0
                
                        Temp_JJA_new(i,j,k) = Temp_JJA(i,j,k);
                        Salt_JJA_new(i,j,k) = Salt_JJA(i,j,k);
                        new_vol(i,j,k)=nv(i,j,k);
%                     end
%                  end
%              end
                  
            
        end
    end
  end


% Nx_new = size(Temp_JJA_new,1);
% Ny_new = size(Temp_JJA_new,2);
% Nz_new = size(Temp_JJA_new,3);

% theta_1d = reshape(Temp_JJA_new,[1],Nx_new*Ny_new*Nz_new);
% salt_1d = reshape(Salt_JJA_new,[1],Nx_new*Ny_new*Nz_new);
% new_vol = reshape(new_vol,[1],Nx_new*Ny_new*Nz_new);

% Salt_JJA_new(:,:,1:10)=[];
% Temp_JJA_new(:,:,1:10)=[];
% new_vol(:,:,1:10)=[];
new_vol(new_vol==0)=NaN;


%%%regrid by volume
salt_grid=33:.005:35;
temp_grid=-3:.005:1.5;
new_v_3445=zeros(size(salt_grid,2),size(temp_grid,2));

for i = 1:Nx
    for j = 1:Ny
       for k = 1:Nr-10
           if Salt_JJA_new(i,j,k)>0
            dist_s3445    = abs(Salt_JJA_new(i,j,k) - salt_grid);
            minDist_s3445 = min(dist_s3445);
            idx_s3445     = find(dist_s3445 == minDist_s3445);        
            
            dist_t_3445    = abs(Temp_JJA_new(i,j,k) - temp_grid);
            minDist_t3445 = min(dist_t_3445);
            idx_t3445     = find(dist_t_3445 == minDist_t3445);    
 
             new_v_3445((idx_s3445),(idx_t3445))=new_v_3445((idx_s3445),(idx_t3445))+(nv(i,j,k));              
             
           end
       end
    end
end

new_v_3445(new_v_3445==0)=NaN;
[s,t]=meshgrid(salt_grid,temp_grid);



%%%getting rid of zeros
% theta_1d2=theta_1d(~theta_1d==0);
% salt_1d2=salt_1d(~salt_1d==0);
% new_vol2=new_vol(~new_vol==0);
% new_vol2=log(new_vol2);
% new_vol2=([17.5:.2:25]);
cmap=jet;
figure(1)
clf
h=pcolor(s,t,new_v_3445');
shading interp
m = colorbar;
colormap(jet);
l = (min(min(new_v_3445)));
h=100;
caxis([10 1000]);
ylabel(m,'Log (Volume) m$^3$','interpreter','latex','fontsize',20)
xlim([33.5 35]);
xticks(33.5:.2:36);
% caxis([18 22]);
title('JJA Temperature/Salinity, 34.9 psu Restoring', 'interpreter','latex','FontSize',16);
xlabel('Salinity (psu)', 'interpreter','latex','FontSize',18);
ylabel('Temperature ($^{o}$C)', 'interpreter','latex','FontSize',18);
b = gca; legend(b,'off');



pt_max =1+.6;
pt_min = -3;
ss_max = 36;
ss_min = 33-1;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/100;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/100;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% Calculate potential density
pd = densmdjwf(SS_grid,PT_grid,0)-1000;
hold on
[C,f] = contour(SS_grid,PT_grid,pd,27.0:.1:33.0,'EdgeColor','k');
clabel(C,f)
% 
% 
ylim([-3 1.5]);
xlim([33.5 35]);
% 
% 
% 
