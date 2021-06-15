%%%%%%%%
%%% attempt to validate with KN section @ 17W




run ../setExpname;
run ../loadexp;
% deltaT=440;
% nIter0 = 424145;
Ny_orig = Ny;
yy_orig = yy;
%%% Ferequency of diagnostic output
dumpFreq = abs(diag_frequency(14));
nDumps = round(nTimeSteps*deltaT/dumpFreq);

dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);


Nyears = 9;

num = Nyears*12;

%%% now reading in 17W data
% load KN_section.mat
distance_shelf = ncread('KappNorvegiaCLM.nc','distance');
depth = ncread('KappNorvegiaCLM.nc','pressure');
salt_KN = ncread('KappNorvegiaCLM.nc','salt');
temp_KN = ncread('KappNorvegiaCLM.nc','ptemp');



%%%% finding the starting longitude in the model so that we can convert the
%%%% distance into longitude
start_lon = dsearchn(xx,-17);
for n = 1:Ny
    if (SHELFICEtopo(start_lon,n)==0 )&& (SHELFICEtopo(start_lon,n)-bathy(start_lon,n)) > 0 
        start = n+2;
        break
    end
end




%%%Now we turn the distances into longitudes using the starting lon
num_lats = 51;
lats_KN = zeros(1,num_lats);
for n = 1:num_lats
    lats_KN(1,n) = yy(start)+((distance_shelf(n)/1000)/110);
end

%%%% To use for our model data later
 min_lat = start;
 max_lat = dsearchn(yy',max(lats_KN));
 range = (max_lat-min_lat);
 yy((1:start)) = [];
 yy(range:end) = [];


nDepths = 101; 
nmonths = 12;
[depth2,lons_KN2] = meshgrid(depth,lats_KN);
    
%%% now we interpolate the data onto our grid...in depth/lat space;

    
 salt_KN2 = zeros(num_lats,length(zz),nmonths);
 for n=1:num_lats
    for k = 1:nmonths
        salt_temp = salt_KN(:,:,k);
        salt_KN2(n,:,k) = interp1(depth,salt_temp(:,n),squeeze(-zz),'linear');
   end
 end
 


%%%%% option for seasonal salt

start_yr = 9;
end_yr = 18;
Nyrs = end_yr-start_yr;
% 

Start_dump = start_yr*12;
end_dump = end_yr*12;

%%%%Chose an arbitrary time to output 
salt_all = NaN(Nx,Ny,Nr,end_dump-Start_dump);
temp_all = NaN(Nx,Ny,Nr,end_dump-Start_dump);

%%% n = 1:nDumps (but experiment hasn't finished yet %%%
for n=1:(end_dump-Start_dump+1)
  
  Salt = rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n+Start_dump));
  salt_all(:,:,:,n) = Salt;   
  
  Temp = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n+Start_dump));
  temp_all(:,:,:,n) = Temp;   
end

Seasonal_Salt = zeros(size(salt_all,1),size(salt_all,2),size(salt_all,3),Nyrs,4);
Seasonal_Temp = zeros(size(temp_all,1),size(temp_all,2),size(temp_all,3),Nyrs,4);


for n=1:Nyrs
  Seasonal_Salt(:,:,:,n,1) = mean(salt_all(:,:,:,(n-1)*12+[12 1 2]),4);
  Seasonal_Salt(:,:,:,n,2) = mean(salt_all(:,:,:,(n-1)*12+(3:5)),4);
  Seasonal_Salt(:,:,:,n,3) = mean(salt_all(:,:,:,(n-1)*12+(6:8)),4);
  Seasonal_Salt(:,:,:,n,4) = mean(salt_all(:,:,:,(n-1)*12+(9:11)),4);
  
  
  Seasonal_Temp(:,:,:,n,1) = mean(temp_all(:,:,:,(n-1)*12+[12 1 2]),4);
  Seasonal_Temp(:,:,:,n,2) = mean(temp_all(:,:,:,(n-1)*12+(3:5)),4);
  Seasonal_Temp(:,:,:,n,3) = mean(temp_all(:,:,:,(n-1)*12+(6:8)),4);
  Seasonal_Temp(:,:,:,n,4) = mean(temp_all(:,:,:,(n-1)*12+(9:11)),4);
end

%%% turn into seasonal mean, avg out year
Salt_S = squeeze(mean(Seasonal_Salt,4));
Salt_S = squeeze(Salt_S(start_lon,:,:,:));
Seastemp = squeeze(mean(Seasonal_Temp,4));
Seastemp = squeeze(Seastemp(start_lon,:,:,:));


Salt_S(1:start,:,:) = [];
Salt_S(range:end,:,:)=[];

Seastemp(1:start,:,:) = [];
Seastemp(range:end,:,:)=[];

Ny = size(Salt_S,1);

Oursalt_s = zeros(num_lats,length(zz),4);
for n=1:Ny
   for v = 1:length(zz)
        for k = 1:4
            s_temp = Salt_S(:,:,k);
              Oursalt_s(:,v,k) = interp1(yy,s_temp(:,v)',lats_KN,'linear');
        end
   
    
   end
 end





%%%%%%%% Back to the total mean.....


 Msalt_MAM = Oursalt_s(:,:,2);
 Msalt_SON = Oursalt_s(:,:,4);

 
 salt_KN2(salt_KN2==0)= NaN;
    
 KNs_JJA = nanmean(salt_KN2(:,:,(3:5)),3);
 KNs_DJF = nanmean(salt_KN2(:,:,9:11),3);

    
 salt_KN2avg = mean(salt_KN2,3);

    

%%%%%% moving on to temperature

temp_KN2 = zeros(num_lats,length(zz),nmonths);
 for n=1:num_lats
   for v = 1:nDepths
    for k = 1:nmonths
        temp = temp_KN(:,:,k);
%         [idx] = find(~isnan(temp(n,:)));
%         [b,i,j]=unique(depth2(n,idx));
        temp_KN2(n,:,k) = interp1(depth,temp(:,n),squeeze(-zz)','linear');
    end
   end
 end


 

 
 
Ourtemp = zeros(num_lats,length(zz),4);
for n = 1:Ny
    for v = 1:length(zz)
        for k = 1:4
            temp2 = Seastemp(:,:,k);

             Ourtemp(:,v,k) = interp1(yy,temp2(:,v)',lats_KN,'linear');
        end
    
    end
end

 
Ourtemp(Ourtemp==0)=NaN;


temp_KN2(temp_KN2==0) = NaN;
temp_KN2avg = nanmean(temp_KN2,3);
Mtemp_MAM = Ourtemp(:,:,2);
Mtemp_SON = Ourtemp(:,:,4);

KNt_DJF = nanmean(temp_KN2(:,:,(3:5)),3);
KNt_JJA = nanmean(temp_KN2(:,:,(9:11)),3);
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make a plot of temperature


%%%%%% Choose a depth ~ Surface and 500m






%%%%%% New meshgrid for plotting
[yg,gz] = meshgrid(lats_KN,zz);



%%%% plotting our model data on original model grid

[yz,zy] = meshgrid(yy_orig,zz);


% MTemp_DJF(KNt_DJF==0)=NaN;
% MTemp_JJA(KNt_JJA==0)=NaN;
% Msalt_DJF(KNs_DJF==0)=NaN;
figure(1)

clf
ww=subplot(2,2,1)
pcolor(yg,gz,Mtemp_MAM'), shading interp
s = set(gca,'position',[0.08 0.53 .42 .40]);
e =colorbar
colormap(e,jet(20))
colormap(ww,jet(20))
caxis([-2 1])
ylim([-1000 0]);
title('Mean MAM Model Temperature ($^\circ$C)','interpreter','latex','FontSize',12');


w = subplot(2,2,2)
%%%%figure plotting model output on A12 track
pcolor(yg,gz,(KNt_DJF-Mtemp_MAM)'), shading interp
set(gca,'position',[0.55 0.53 .42 .40]);
e1 =colorbar
colormap(e1,redblue(20))
colormap(w,redblue(20))
caxis([-2 2])
ylim([-1000 0]);
title('MAM Theta Anomaly (Observed minus Model) ($^\circ$C)','interpreter','latex','FontSize',12');

w2=subplot(2,2,3)
%%%figure plotting model output on A12 track
pcolor(yg,gz,Mtemp_SON'), shading interp
set(gca,'position',[0.08 0.06 .42 .40]);
e2 =colorbar
colormap(e2,jet(20))
colormap(w2,jet(20))
caxis([-2 1])
xlabel('Latitude','interpreter','latex','FontSize',12');
ylim([-1000 0])
title('Mean SON Model Temperature ($^\circ$C)','interpreter','latex','FontSize',12');


w3=subplot(2,2,4)
%%%%figure plotting model output on A12 track
pcolor(yg,gz,(KNt_JJA-Mtemp_SON)'), shading interp
set(gca,'position',[0.55 0.06 .42 .40]);
e3 =colorbar
colormap(e3,redblue(20))
colormap(w3,redblue(20))
caxis([-2 2])
xlabel('Latitude','interpreter','latex','FontSize',12');
ylim([-1000 0])
title('SON Theta Anomaly (Observed minus Model) ($^\circ$C)','interpreter','latex','FontSize',12'); 

figure(2)

clf
ww=subplot(2,2,1)
pcolor(yg,gz,Msalt_MAM'), shading interp
s = set(gca,'position',[0.08 0.53 .42 .40]);
e =colorbar
colormap(e,pmkmp(20))
colormap(ww,pmkmp(20))
caxis([34.2 34.9])
ylim([-1000 0]);
title('Mean MAM Model Salinity at Kappa Novegia (psu)','interpreter','latex','FontSize',12');


w = subplot(2,2,2)
%%%%figure plotting model output on A12 track
pcolor(yg,gz,(KNs_DJF-Msalt_MAM)'), shading interp
set(gca,'position',[0.55 0.53 .42 .40]);
e1 =colorbar
colormap(e1,redblue(20))
colormap(w,redblue(20))
caxis([-.3 .3])
ylim([-1000 0]);
title('MAM Salinity Anomaly (Observed minus Model) (psu)','interpreter','latex','FontSize',12');

w2=subplot(2,2,3)
%%%figure plotting model output on A12 track
pcolor(yg,gz,Msalt_SON'), shading interp
set(gca,'position',[0.08 0.06 .42 .40]);
e2 =colorbar
colormap(e2,pmkmp(20))
colormap(w2,pmkmp(20))
caxis([34.2 34.9])
xlabel('Latitude','interpreter','latex','FontSize',12');
ylim([-1000 0])
title('Mean SON Model Salinity (psu)','interpreter','latex','FontSize',12');


w3=subplot(2,2,4)
%%%%figure plotting model output on A12 track
pcolor(yg,gz,(KNs_JJA-Msalt_SON)'), shading interp
set(gca,'position',[0.55 0.06 .42 .40]);
e3 =colorbar
colormap(e3,redblue(20))
colormap(w3,redblue(20))
caxis([-.3 .3])
xlabel('Latitude','interpreter','latex','FontSize',12');
ylim([-1000 0])
title('SON Salinity Anomaly (Observed minus Model) (psu)','interpreter','latex','FontSize',12'); 
