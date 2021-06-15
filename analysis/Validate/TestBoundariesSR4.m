%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Validating Boundary Conditions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing Boundary Conditions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% First have to load experiment data
%%%%% loading SR4 data
run ../../newexp/defineGrid.m
run ../setExpname;
run ../loadexp;
%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(14));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);


tmin = 16*86400*365;
tmax = 25*86400*365;
Theta = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);
Salt = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);




datadir ='/data1/MITgcm_WS/data/SR4/';

%%%% SR04 latitude and longitude

SR4_latitude_89 = ncread(fullfile(datadir,'ANT_1989_VIII2_Fis.nc'),'Latitude');
SR4_longitude_89 = ncread(fullfile(datadir,'ANT_1989_VIII2_Fis.nc'),'Longitude');

%%%%%% extracting SR4 temperature


N_Stations = 64;

%%%%%% FIRST PUTTING SR4 DATA ON OUR VERTICAL GRID
saltsr4_89 = ncread(fullfile(datadir,'ANT_1989_VIII2_Fis.nc'),'VAR_3');
tempsr4_89 = ncread(fullfile(datadir,'ANT_1989_VIII2_Fis.nc'),'VAR_2');
depthsr4_89 = ncread(fullfile(datadir,'ANT_1989_VIII2_Fis.nc'),'VAR_1');

temp_newsr489 = zeros(length(zz),size(depthsr4_89,2));
for n=1:size(depthsr4_89,2)
  idx = find(~isnan(depthsr4_89(:,n)));
  temp_newsr489(:,n) = interp1(depthsr4_89(idx,n),tempsr4_89(idx,n),squeeze(-zz),'linear');
end

salt_newsr489 = zeros(length(zz),size(depthsr4_89,2));
for n=1:size(depthsr4_89,2)
  idx = find(~isnan(depthsr4_89(:,n)));
  salt_newsr489(:,n) = interp1(depthsr4_89(idx,n),saltsr4_89(idx,n),squeeze(-zz),'linear');
end


%%%%%%%%%%% Other cruises
%%%% 1993

SR4Cruise_93= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1993_X7_Fis.nc');
depthsr4_93 = ncread(fullfile(datadir,'ANT_1993_X7_Fis.nc'),'VAR_1');
tempsr4_93 = ncread(fullfile(datadir,'ANT_1993_X7_Fis.nc'),'VAR_2');
saltsr4_93 = ncread(fullfile(datadir,'ANT_1993_X7_Fis.nc'),'VAR_3');
SR4_2_Lon = ncread(fullfile(datadir,'ANT_1993_X7_Fis.nc'),'Longitude');

temp_newsr4_93 = zeros(length(zz),size(depthsr4_93,2));
for n=1:size(depthsr4_93,2)
  idx = find(~isnan(depthsr4_93(:,n)));
  temp_newsr4_93(:,n) = interp1(depthsr4_93(idx,n),tempsr4_93(idx,n),squeeze(-zz),'linear');
end

Theta_93_new_lon = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_temp_93 = temp_newsr4_93(n,:);
   Theta_93_new_lon(n,:) = interp1(SR4_2_Lon,sr4_temp_93,SR4_longitude_89,'linear');
end


salt_newsr4_93 = zeros(length(zz),size(depthsr4_93,2));
for n=1:size(depthsr4_93,2)
  idx = find(~isnan(depthsr4_93(:,n)));
  salt_newsr4_93(:,n) = interp1(depthsr4_93(idx,n),saltsr4_93(idx,n),squeeze(-zz),'linear');
end

Salt_93_new_lon = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_salt_93 = salt_newsr4_93(n,:);
   Salt_93_new_lon(n,:) = interp1(SR4_2_Lon,sr4_salt_93,SR4_longitude_89,'linear');
end

%%%%%%%%% 1996 cruise



SR4_96_Lon = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'Longitude');
SR4_96_Lat = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'Latitude');
saltsr4_96 = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'VAR_3');

depthsr4_96 = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'VAR_1');
tempsr4_96 = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'VAR_2');



temp_96 = zeros(length(zz),size(depthsr4_96,2));
for n=1:size(depthsr4_96,2)
  idx = find(~isnan(depthsr4_96(:,n)));
  temp_96(:,n) = interp1(depthsr4_96(idx,n),tempsr4_96(idx,n),squeeze(-zz),'linear');
end

Theta_96LO = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_temp_96 = temp_96(n,:);
   Theta_96LO(n,:) = interp1(SR4_96_Lon,sr4_temp_96,SR4_longitude_89,'linear');
end


salt_96 = zeros(length(zz),size(depthsr4_96,2));
for n=1:size(depthsr4_96,2)
  idx = find(~isnan(depthsr4_96(:,n)));
  salt_96(:,n) = interp1(depthsr4_96(idx,n),saltsr4_96(idx,n),squeeze(-zz),'linear');
end

Salt_96LO = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_salt_96 = salt_96(n,:);
   Salt_96LO(n,:) = interp1(SR4_96_Lon,sr4_salt_96,SR4_longitude_89,'linear');
end



%%%%%%%%%%%%%%%%%%%%%% 2008 cruise


depthsr4_08 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_1');
tempsr4_08 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_2');
saltsr4_08 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_3');

SR4_08_Lon = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'Longitude');
SR4_08_Lat = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'Latitude');


temp_08 = zeros(length(zz),size(depthsr4_08,2));
for n=1:size(depthsr4_08,2)
  idx = find(~isnan(depthsr4_08(:,n)));
  temp_08(:,n) = interp1(depthsr4_08(idx,n),tempsr4_08(idx,n),squeeze(-zz),'linear');
end

Theta_08LO = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_temp_08 = temp_08(n,:);
   Theta_08LO(n,:) = interp1(SR4_08_Lon,sr4_temp_08,SR4_longitude_89,'linear');
end


salt_08 = zeros(length(zz),size(depthsr4_08,2));
for n=1:size(depthsr4_08,2)
  idx = find(~isnan(depthsr4_08(:,n)));
  salt_08(:,n) = interp1(depthsr4_08(idx,n),saltsr4_08(idx,n),squeeze(-zz),'linear');
end

Salt_08LO = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_salt_08 = salt_08(n,:);
   Salt_08LO(n,:) = interp1(SR4_08_Lon,sr4_salt_08,SR4_longitude_89,'linear');
end




%%%%%%%%%%%% 1990

SR4Cruise_90= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1990_IX2_Fis.nc');   
depthsr4_90 = ncread(fullfile(datadir,'ANT_1990_IX2_Fis.nc'),'VAR_1');
tempsr4_90 = ncread(fullfile(datadir,'ANT_1990_IX2_Fis.nc'),'VAR_2');
SR4_90_Lon = ncread(fullfile(datadir,'ANT_1990_IX2_Fis.nc'),'Longitude');
SR4_90_Lat = ncread(fullfile(datadir,'ANT_1990_IX2_Fis.nc'),'Latitude');
saltsr4_90 = ncread(fullfile(datadir,'ANT_1990_IX2_Fis.nc'),'VAR_3');


temp_90 = zeros(length(zz),size(depthsr4_90,2));
for n=1:size(depthsr4_90,2)
  idx = find(~isnan(depthsr4_90(:,n)));
  [b,i,j]=unique(depthsr4_90(idx,n));
  temp_90(:,n) = interp1(b,tempsr4_90(i,n),squeeze(-zz),'linear');
end

Theta_90LO = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_temp_90 = temp_90(n,:);
   [b,i,j]=unique(SR4_90_Lon);
   Theta_90LO(n,:) = interp1(b,sr4_temp_90(i),SR4_longitude_89,'linear');
end


salt_90 = zeros(length(zz),size(depthsr4_90,2));
for n=1:size(depthsr4_90,2)
  idx = find(~isnan(depthsr4_90(:,n)));
  [b,i,j]=unique(depthsr4_90(idx,n));
  salt_90(:,n) = interp1(b,saltsr4_90(i,n),squeeze(-zz),'linear');
end

Salt_90LO = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_salt_90 = salt_90(n,:);
   [b,i,j]=unique(SR4_90_Lon);
   Salt_90LO(n,:) = interp1(b,sr4_salt_90(i),SR4_longitude_89,'linear');
end


%%%%%% 1998

SR4Cruise_98= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1998_XV4_Fis.nc'); 
SR4_98_Lon = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'Longitude');
SR4_98_Lat = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'Latitude');
depthsr4_98 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_1');
tempsr4_98 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_2');
saltsr4_98 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_3');




temp_98 = zeros(length(zz),size(depthsr4_98,2));
for n=1:size(depthsr4_98,2)
  idx = find(~isnan(depthsr4_98(:,n)));
  temp_98(:,n) = interp1(depthsr4_98(idx,n),tempsr4_98(idx,n),squeeze(-zz),'linear');
end

Theta_98LO = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_temp_98 = temp_98(n,:);
   Theta_98LO(n,:) = interp1(SR4_98_Lon,sr4_temp_98,SR4_longitude_89,'linear');
end

salt_98 = zeros(length(zz),size(depthsr4_98,2));
for n=1:size(depthsr4_98,2)
  idx = find(~isnan(depthsr4_98(:,n)));
  salt_98(:,n) = interp1(depthsr4_98(idx,n),saltsr4_98(idx,n),squeeze(-zz),'linear');
end

Salt_98LO = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_salt_98 = salt_98(n,:);
   Salt_98LO(n,:) = interp1(SR4_98_Lon,sr4_salt_98,SR4_longitude_89,'linear');
end





%%%%%%% 2011

SR4Cruise9= fullfile('/data1/MITgcm_WS/data/SR4/ANT_2011_XXVII2_Fis.nc');
depthsr4_11 = ncread(fullfile(datadir,'ANT_2011_XXVII2_Fis.nc'),'VAR_1');
SR4_11_Lon = ncread(fullfile(datadir,'ANT_2011_XXVII2_Fis.nc'),'Longitude');
SR4_11_Lat = ncread(fullfile(datadir,'ANT_2011_XXVII2_Fis.nc'),'Latitude');
tempsr4_11 = ncread(fullfile(datadir,'ANT_2011_XXVII2_Fis.nc'),'VAR_2');
saltsr4_11 = ncread(fullfile(datadir,'ANT_2011_XXVII2_Fis.nc'),'VAR_3');




temp_11 = zeros(length(zz),size(depthsr4_11,2));
for n=1:size(depthsr4_11,2)
  idx = find(~isnan(depthsr4_11(:,n)));
  temp_11(:,n) = interp1(depthsr4_11(idx,n),tempsr4_11(idx,n),squeeze(-zz),'linear');
end

Theta_11LO = zeros(length(zz),N_Stations);
for n=1:length(zz)
   sr4_temp_11 = temp_11(n,:);
   Theta_11LO(n,:) = interp1(SR4_11_Lon,sr4_temp_11,SR4_longitude_89,'linear');
end


% salt_11 = zeros(length(zz),size(depthsr4_11,2));
% for n=1:size(depthsr4_11,2)
%   idx = find(~isnan(depthsr4_11(:,n)));
%   salt_11(:,n) = interp1(saltsr4_11(idx,n),saltsr4_11(idx,n),squeeze(-zz),'linear');
% end
% 
% Salt_11LO = zeros(length(zz),N_Stations);
% for n=1:length(zz)
%    sr4_salt_11 = salt_11(n,:);
%    Salt_11LO(n,:) = interp1(SR4_11_Lon,sr4_salt_11,SR4_longitude_89,'linear');
% end


%%%%%% Find the mean of all of the cruise data

z = cat(3,Theta_11LO,Theta_96LO,Theta_98LO,Theta_93_new_lon,temp_newsr489,Theta_90LO);
Mean_CTD_Temp = nanmean(z,3);


s = cat(3,Salt_96LO,Salt_98LO,Salt_93_new_lon,salt_newsr489,Salt_90LO);
Mean_CTD_Salt = nanmean(s,3);


%%%%%%% trying to put OUR temperature results on SR4
%%%%%%% Using 1989 cruise as defined SR4 lat/lon points
SR4_new=zeros(size(SR4_latitude_89,1),1);
for i = 1:size(SR4_latitude_89,1)
    if SR4_latitude_89(i)>-65
        SR4_new(i) = SR4_latitude_89(i)-2;
    else
        SR4_new(i)= SR4_latitude_89(i);
    end
end
    

Theta_new_lon = zeros(N_Stations,length(ymc),length(zz));
for k = 1:Ny
    for n=1:length(zz)
        sr4_temp = Theta(:,k,n);
        Theta_new_lon(:,k,n) = interp1(xmc,sr4_temp,SR4_longitude_89,'linear');
    end
end

Theta_MOD = zeros(N_Stations,N_Stations,length(zz));
for k = 1:N_Stations
    for n=1:length(zz)
        sr4_temp_lat = Theta_new_lon(k,:,n);
        Theta_MOD(k,:,n) = interp1(ymc,sr4_temp_lat,SR4_new,'linear');
    end
end




Salt_new_lon = zeros(N_Stations,length(ymc),length(zz));
for k = 1:Ny
    for n=1:length(zz)
        sr4_salt = Salt(:,k,n);
        Salt_new_lon(:,k,n) = interp1(xmc,sr4_salt,SR4_longitude_89,'linear');
    end
end

Salt_MOD = zeros(N_Stations,N_Stations,length(zz));
for k = 1:N_Stations
    for n=1:length(zz)
        sr4_salt_lat = Salt_new_lon(k,:,n);
        Salt_MOD(k,:,n) = interp1(ymc,sr4_salt_lat,SR4_new,'linear');
    end
end

%%%%%% Interpolate hFacC?
hFacC_new_lon = zeros(N_Stations,length(ymc),length(zz));
for k = 1:Ny
    for n=1:length(zz)
        sr4_hFacC = hFacC(:,k,n);
        hFacC_new_lon(:,k,n) = interp1(xmc,sr4_hFacC,SR4_longitude_89,'linear');
    end
end

hFacC_MOD = zeros(N_Stations,N_Stations,length(zz));
for k = 1:N_Stations
    for n=1:length(zz)
        sr4_hFacC_lat = hFacC_new_lon(k,:,n);
        hFacC_MOD(k,:,n) = interp1(ymc,sr4_hFacC_lat,SR4_new,'linear');
    end
end








%%% make a meshgrid of 1989 Lon and our vertical grid

[SR4_LO,zzm]=meshgrid(SR4_longitude_89,zz);

Mean_CTD_Temp = Mean_CTD_Temp';
Theta_MOD = squeeze(Theta_MOD(:,18,:));

Mean_CTD_Salt = Mean_CTD_Salt';
Salt_MOD = squeeze(Salt_MOD(:,18,:));

hFacC_MOD = squeeze(hFacC_MOD(:,18,:));

Salt_MOD(hFacC_MOD<.6)=NaN;
Theta_MOD(hFacC_MOD<.6)=NaN;
Mean_CTD_Salt(Mean_CTD_Salt==0)=NaN;
Mean_CTD_Temp(Mean_CTD_Temp==0)=NaN;


figure(1)
clf
w1=subplot(2,2,1)

pcolor(SR4_LO',zzm',Theta_MOD), shading interp
s = set(gca,'position',[0.055 0.53 .42 .37]);
colormap(w1,jet(30))
h = colorbar
colormap(h,jet(30))
caxis([-2 1]);
ylim([-1000 0]);
% xlabel('Longitude','interpreter','latex','FontSize',12');
title('Mean Model Temperature ($^\circ$C), SR4','interpreter','latex','FontSize',13');
set(gca,'fontsize',15);

hold on

%%%%figure plotting model output on SR4 track
w=subplot(2,2,2);
pcolor(SR4_LO',zzm',(Mean_CTD_Temp-Theta_MOD)), shading interp
set(gca,'position',[0.55 0.53 .42 .37]);
colormap(w,redblue(20))

b=colorbar
colormap(b,redblue(20))
caxis([-2 2])
ylim([-1000 0]);
% xlabel('Longitude','interpreter','latex','FontSize',12');
title('SR4 Theta Anomaly ($^\circ$C) (Observed minus Model)','interpreter','latex','FontSize',13');
set(gca,'fontsize',15);

hold on
%%%%figure plotting model output on SR4 track
w2=subplot(2,2,3)
pcolor(SR4_LO',zzm',Salt_MOD), shading interp
colormap(w2,redblue(20))
set(gca,'position',[0.055 0.08 .42 .37]);
g = colorbar
colormap(g,pmkmp(20))
colormap(w2,pmkmp(20))
caxis([34.3 34.8])
ylim([-1000 0]);
xlabel('Longitude','interpreter','latex','FontSize',12');
title('Mean Model Salinity (psu), SR4','interpreter','latex','FontSize',13');
set(gca,'fontsize',15);


hold on
%%%%figure plotting model output on SR4 track
w4=subplot(2,2,4)
pcolor(SR4_LO',zzm',(Mean_CTD_Salt-Salt_MOD)), shading interp
set(gca,'position',[0.55 0.08 .42 .37]);
z=colorbar
colormap(z,redblue(20))
colormap(w4,redblue(20))
caxis([-.2 .2])
ylim([-1000 0]);
xlabel('Longitude','interpreter','latex','FontSize',12');
title('SR4 Salt Anomaly (psu) (Observed minus Model)','interpreter','latex','FontSize',13');
set(gca,'fontsize',15);









