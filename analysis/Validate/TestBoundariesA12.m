%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Boundaries A12
%%%%%%%%%%%%%%%%%%%%%%% % Similar to SR4 script

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
%%% Ferequency of diagnostic output
dumpFreq = abs(diag_frequency(14));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);



%%%%%% Loading our Temp/Salinity data from respective experiment

tmin = 9*86400*365;
tmax = 18*86400*365;
Theta = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
Salt = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);




datadir ='/data1/MITgcm_WS/data/A12/';

%%%% A12 latitude and longitude

Lon_limit=1;
A12_latitude_02 = ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'Latitude');
A12_longitude_02= ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'Longitude');
data_remove = (A12_longitude_02(:,:)>1);
A12_longitude_02(data_remove)=[];
A12_latitude_02(data_remove)=[];



%%%%%% extracting SR4 temperature


N_Stations = 37;

%%%%%% FIRST PUTTING SR4 DATA ON OUR VERTICAL GRID
salt_12 = ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'VAR_3');
temp_12 = ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'VAR_2');
depth_12 = ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'VAR_1');

salt_12(:,data_remove)=[];
temp_12(:,data_remove)=[];
depth_12(:,data_remove)=[];


temp_a12 = zeros(length(zz),size(depth_12,2));
for n=1:size(depth_12,2)
  idx = find(~isnan(depth_12(:,n)));
  [b,i,j]=unique(depth_12(idx,n));
  temp_a12(:,n) = interp1((b),temp_12(i,n),squeeze(-zz),'linear');
end

salt_a12 = zeros(length(zz),size(depth_12,2));
for n=1:size(depth_12,2)
  idx = find(~isnan(depth_12(:,n)));
  [b,i,j]=unique(depth_12(idx,n));
  salt_a12(:,n) = interp1(b,salt_12(i,n),squeeze(-zz),'linear');
end


%%%%%%%%%% Other cruises
%%% 1992

depth_92 = ncread(fullfile(datadir,'ANT_1992_X4_Fis.nc'),'VAR_1');  
temp_92 = ncread(fullfile(datadir,'ANT_1992_X4_Fis.nc'),'VAR_2');
salt_92 = ncread(fullfile(datadir,'ANT_1992_X4_Fis.nc'),'VAR_3');
A12_2_Lat = ncread(fullfile(datadir,'ANT_1992_X4_Fis.nc'),'Latitude');

salt_92(:,data_remove)=[];
temp_92(:,data_remove)=[];
depth_92(:,data_remove)=[];
A12_2_Lat(data_remove)=[];

temp_newa12_92 = zeros(length(zz),size(depth_92,2));
for n=1:size(depth_92,2)
  idx = find(~isnan(depth_92(:,n)));
  [b,i,j]=unique(depth_92(idx,n));
  temp_newa12_92(:,n) = interp1(b,temp_92(i,n),squeeze(-zz),'linear');
end

Theta_92 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_temp_92 = temp_newa12_92(n,:);
   [b,i,j]=unique(A12_2_Lat);
   Theta_92(n,:) = interp1(b,a12_temp_92(i),A12_latitude_02,'linear');
end


salt_newa12_92 = zeros(length(zz),size(depth_92,2));
for n=1:size(depth_92,2)
  idx = find(~isnan(depth_92(:,n)));
  [b,i,j]=unique(depth_92(idx,n));
  salt_newa12_92(:,n) = interp1(b,salt_92(i,n),squeeze(-zz),'linear');
end

Salt_92 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_salt_92 = salt_newa12_92(n,:);
   [b,i,j]=unique(A12_2_Lat);
   Salt_92(n,:) = interp1(b,a12_salt_92(i),A12_latitude_02,'linear');
end


% %%%%%%%%%%%%%%%%%%%
% %%%% 1996 %%%%%%%%%

depth_96 = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'VAR_1');  
temp_96 = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'VAR_2');
salt_96 = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'VAR_3');
A12_6_Lat = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'Latitude');

salt_96(:,data_remove)=[];
temp_96(:,data_remove)=[];
depth_96(:,data_remove)=[];
A12_6_Lat(data_remove)=[];


temp_newa12_96 = zeros(length(zz),size(depth_96,2));
for n=1:size(depth_96,2)
  idx = find(~isnan(depth_96(:,n)));
  [b,i,j]=unique(depth_96(idx,n));
  temp_newa12_96(:,n) = interp1(b,temp_96(i,n),squeeze(-zz),'linear');
end

Theta_96 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_temp_96 = temp_newa12_96(n,:);
   [b,i,j]=unique(A12_6_Lat);
   Theta_96(n,:) = interp1(b,a12_temp_96(i),A12_latitude_02,'linear');
end


salt_newa12_96 = zeros(length(zz),size(depth_96,2));
for n=1:size(depth_96,2)
  idx = find(~isnan(depth_96(:,n)));
  [b,i,j]=unique(depth_96(idx,n));
  salt_newa12_96(:,n) = interp1(b,salt_96(i,n),squeeze(-zz),'linear');
end

Salt_96 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_salt_96 = salt_newa12_96(n,:);
   [b,i,j]=unique(A12_6_Lat);
   Salt_96(n,:) = interp1(b,a12_salt_96(i),A12_latitude_02,'linear');
end

% %%%%%%%%%%%%%%%%%%%%%% 1998 cruise
% 

depth_98 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_1');  
temp_98 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_2');
salt_98 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_3');
A12_8_Lat = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'Latitude');
 
salt_98(:,data_remove)=[];
temp_98(:,data_remove)=[];
depth_98(:,data_remove)=[];
A12_8_Lat(data_remove)=[];


temp_newa12_98 = zeros(length(zz),size(depth_98,2));
for n=1:size(depth_98,2)
  idx = find(~isnan(depth_98(:,n)));
  [b,i,j]=unique(depth_98(idx,n));
  temp_newa12_98(:,n) = interp1(b,temp_98(i,n),squeeze(-zz),'linear');
end

Theta_98 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_temp_98 = temp_newa12_98(n,:);
   [b,i,j]=unique(A12_8_Lat);
   Theta_98(n,:) = interp1(b,a12_temp_98(i),A12_latitude_02,'linear');
end


salt_newa12_98 = zeros(length(zz),size(depth_98,2));
for n=1:size(depth_98,2)
  idx = find(~isnan(depth_98(:,n)));
  [b,i,j]=unique(depth_98(idx,n));
  salt_newa12_98(:,n) = interp1(b,salt_98(i,n),squeeze(-zz),'linear');
end

Salt_98 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_salt_98 = salt_newa12_98(n,:);
   [b,i,j]=unique(A12_6_Lat);
   Salt_98(n,:) = interp1(b,a12_salt_98(i),A12_latitude_02,'linear');
end



%%%%%%%%%%%%%%%%%%%%%% 2000 cruise


depth_00 = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'VAR_1');  
temp_00 = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'VAR_2');
salt_00 = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'VAR_3');

A12_0_Lat = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'Latitude');
A12_0_Lon = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'Longitude');
% 
% salt_00(:,data_remove)=[];
% temp_00(:,data_remove)=[];
% depth_00(:,data_remove)=[];
% A12_0_Lat(data_remove)=[];

temp_newa12_00 = zeros(length(zz),size(depth_00,2));
for n=1:size(depth_00,2)
  idx = find(~isnan(depth_00(:,n)));
  [b,i,j]=unique(depth_00(idx,n));
  temp_newa12_00(:,n) = interp1(depth_00(idx,n),temp_00(idx,n),squeeze(-zz),'linear');
end

Theta_00 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_temp_00 = temp_newa12_00(n,:);
   [b,i,j]=unique(A12_0_Lat);
   Theta_00(n,:) = interp1(b,a12_temp_00(i),A12_latitude_02,'linear');
end


salt_newa12_00 = zeros(length(zz),size(depth_00,2));
for n=1:size(depth_00,2)
  idx = find(~isnan(depth_00(:,n)));
  [b,i,j]=unique(depth_00(idx,n));
  salt_newa12_00(:,n) = interp1(depth_00(idx,n),salt_00(idx,n),squeeze(-zz),'linear');
end

Salt_00 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_salt_00 = salt_newa12_00(n,:);
   [b,i,j]=unique(A12_0_Lat);
   Salt_00(n,:) = interp1(b,a12_salt_00(i),A12_latitude_02,'linear');
end


%%%%%%%%%%%%%%%%%%%%%% 2005 cruise


depth_05 = ncread(fullfile(datadir,'ANT_2005_XXII3_Fis.nc'),'VAR_1');  
temp_05 = ncread(fullfile(datadir,'ANT_2005_XXII3_Fis.nc'),'VAR_2');
salt_05 = ncread(fullfile(datadir,'ANT_2005_XXII3_Fis.nc'),'VAR_3');
A12_05_Lat = ncread(fullfile(datadir,'ANT_2005_XXII3_Fis.nc'),'Latitude');

salt_05(:,data_remove)=[];
temp_05(:,data_remove)=[];
depth_05(:,data_remove)=[];
A12_05_Lat(data_remove)=[];

temp_newa12_05 = zeros(length(zz),size(depth_05,2));
for n=1:size(depth_05,2)
  idx = find(~isnan(depth_05(:,n)));
  [b,i,j]=unique(depth_05(idx,n));
  temp_newa12_05(:,n) = interp1(b,temp_05(i,n),squeeze(-zz),'linear');
end

Theta_05 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_temp_05 = temp_newa12_05(n,:);
   [b,i,j]=unique(A12_05_Lat);
   Theta_05(n,:) = interp1(b,a12_temp_05(i),A12_latitude_02,'linear');
end


salt_newa12_05 = zeros(length(zz),size(depth_05,2));
for n=1:size(depth_05,2)
  idx = find(~isnan(depth_05(:,n)));
  [b,i,j]=unique(depth_05(idx,n));
  salt_newa12_05(:,n) = interp1(b,salt_05(i,n),squeeze(-zz),'linear');
end

Salt_05 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_salt_05 = salt_newa12_05(n,:);
   [b,i,j]=unique(A12_05_Lat);
   Salt_05(n,:) = interp1(b,a12_salt_05(i),A12_latitude_02,'linear');
end




%%%%%%%%%%%%%%%%%%%%%% 2008 cruise


depth_08 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_1');  
temp_08 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_2');
salt_08 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_3');
A12_08_Lat = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'Latitude');

salt_08(:,data_remove)=[];
temp_08(:,data_remove)=[];
depth_08(:,data_remove)=[];
A12_08_Lat(data_remove)=[];

temp_newa12_08 = zeros(length(zz),size(depth_08,2));
for n=1:size(depth_08,2)
  idx = find(~isnan(depth_08(:,n)));
  [b,i,j]=unique(depth_08(idx,n));
  temp_newa12_08(:,n) = interp1(b,temp_08(i,n),squeeze(-zz),'linear');
end

Theta_08 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_temp_08 = temp_newa12_08(n,:);
   [b,i,j]=unique(A12_08_Lat);
   Theta_08(n,:) = interp1(b,a12_temp_08(i),A12_latitude_02,'linear');
end

salt_newa12_08 = zeros(length(zz),size(depth_08,2));
for n=1:size(depth_08,2)
  idx = find(~isnan(depth_08(:,n)));
  [b,i,j]=unique(depth_08(idx,n));
  salt_newa12_08(:,n) = interp1(b,salt_08(i,n),squeeze(-zz),'linear');
end

Salt_08 = zeros(length(zz),N_Stations);
for n=1:length(zz)
   a12_salt_08 = salt_newa12_08(n,:);
   [b,i,j]=unique(A12_08_Lat);
   Salt_08(n,:) = interp1(b,a12_salt_08(i),A12_latitude_02,'linear');
end



%%%%%% Find the mean of all of the cruise data

z = cat(3,Theta_92,Theta_00,Theta_05,Theta_08,Theta_96,Theta_98,temp_a12);
Mean_CTD_Temp = nanmean(z,3);


s = cat(3,Salt_92,salt_a12,Salt_00,Salt_05,Salt_08,Salt_96,Salt_98);
Mean_CTD_Salt = nanmean(s,3);


%%%%%%% trying to put OUR temperature results on SR4
%%%%%%% Using 1989 cruise as defined SR4 lat/lon points

N_Stations = 37;
Theta_new_lon = zeros(N_Stations,length(ymc),length(zz));
for k = 1:Ny
    for n=1:length(zz)
        a12_temp = Theta(:,k,n);
        Theta_new_lon(:,k,n) = interp1(xmc,a12_temp,A12_longitude_02,'linear');
    end
end

Theta_MOD = zeros(N_Stations,N_Stations,length(zz));
for k = 1:N_Stations
    for n=1:length(zz)
        a12_temp_lat = Theta_new_lon(k,:,n);
        Theta_MOD(k,:,n) = interp1(ymc,a12_temp_lat,A12_latitude_02,'linear');
    end
end


Theta_MOD2 = zeros(N_Stations,N_Stations,length(zz));
for k = 1:N_Stations
    for n=1:length(zz)
        a12_temp_lat = Theta_new_lon(k,:,n);
        Theta_MOD(k,:,n) = interp1(ymc,a12_temp_lat,A12_latitude_02,'linear');
    end
end




Salt_new_lon = zeros(N_Stations,length(ymc),length(zz));
for k = 1:Ny
    for n=1:length(zz)
        a12_salt = Salt(:,k,n);
        Salt_new_lon(:,k,n) = interp1(xmc,a12_salt,A12_longitude_02,'linear');
    end
end

Salt_MOD = zeros(N_Stations,N_Stations,length(zz));
for k = 1:N_Stations
    for n=1:length(zz)
        a12_salt_lat = Salt_new_lon(k,:,n);
        Salt_MOD(k,:,n) = interp1(ymc,a12_salt_lat,A12_latitude_02,'linear');
    end
end

%%%%%% Interpolate hFacC?
hFacC_new_lon = zeros(N_Stations,length(ymc),length(zz));
for k = 1:Ny
    for n=1:length(zz)
        a12_hFacC = hFacC(:,k,n);
        hFacC_new_lon(:,k,n) = interp1(xmc,a12_hFacC,A12_longitude_02,'linear');
    end
end

hFacC_MOD = zeros(N_Stations,N_Stations,length(zz));
for k = 1:N_Stations
    for n=1:length(zz)
        a12_hFacC_lat = hFacC_new_lon(k,:,n);
        hFacC_MOD(k,:,n) = interp1(ymc,a12_hFacC_lat,A12_latitude_02,'linear');
    end
end








%%% make a meshgrid of lat and our vertical grid

[A12_LA,zzm]=meshgrid(A12_latitude_02,zz);

Mean_CTD_Temp = Mean_CTD_Temp';
Theta_MOD = squeeze(Theta_MOD(end,:,:));

Mean_CTD_Salt = Mean_CTD_Salt';
Salt_MOD = squeeze(Salt_MOD(end,:,:));

hFacC_MOD = squeeze(hFacC_MOD(end,:,:));

Salt_MOD(hFacC_MOD<1)=NaN;
Theta_MOD(hFacC_MOD<1)=NaN;
Mean_CTD_Salt(Mean_CTD_Salt==0)=NaN;
Mean_CTD_Temp(Mean_CTD_Temp==0)=NaN;



%%%%figure plotting model output on A12 track
figure(1)
clf
ww=subplot(2,2,1)
pcolor(A12_LA',zzm',Theta_MOD), shading interp
s = set(gca,'position',[0.08 0.53 .42 .40]);
e =colorbar
colormap(e,jet(20))
colormap(ww,jet(20))
caxis([-2 1])
xlim([-70 -64]);
ylim([-1000 0]);
title('Mean Model Temperature ($^\circ$C), A12','interpreter','latex','FontSize',16');
set(gca,'fontsize',15);


w = subplot(2,2,2)
%%%%figure plotting model output on A12 track
pcolor(A12_LA',zzm',(Mean_CTD_Temp-Theta_MOD)), shading interp
set(gca,'position',[0.55 0.53 .42 .40]);
e1 =colorbar
colormap(e1,redblue(20))
colormap(w,redblue(20))
caxis([-2 2])
xlim([-70 -64]);
ylim([-1000 0]);
title('A12 Theta Anomaly ($^\circ$C) (Observed minus Model)','interpreter','latex','FontSize',16');
set(gca,'fontsize',15);

w2=subplot(2,2,3)
%%%figure plotting model output on A12 track
pcolor(A12_LA',zzm',Salt_MOD), shading interp
set(gca,'position',[0.08 0.06 .42 .40]);
e2 =colorbar
colormap(e2,pmkmp(20))
colormap(w2,pmkmp(20))
caxis([34.2 34.9])
xlabel('Latitude','interpreter','latex','FontSize',12');
xlim([-70 -64]);
ylim([-1000 0])
title('Mean Model Salinity (psu), A12','interpreter','latex','FontSize',16');
set(gca,'fontsize',15);


w3=subplot(2,2,4)
%%%%figure plotting model output on A12 track
pcolor(A12_LA',zzm',(Mean_CTD_Salt-Salt_MOD)), shading interp
set(gca,'position',[0.55 0.06 .42 .40]);
e3 =colorbar
colormap(e3,redblue(20))
colormap(w3,redblue(20))
caxis([-.2 .2])
xlabel('Latitude','interpreter','latex','FontSize',12');
xlim([-70 -64 ]);
ylim([-1000 0])
title('A12 Salt Anomaly (psu) (Observed minus Model)','interpreter','latex','FontSize',13'); 


set(gca,'fontsize',15);

%%%%% finding longitude locations
load KNlocs.mat

% %%%for appendix of WOCE sections
% Thetaappend=Theta(:,:,1);
% Thetaappend(Thetaappend==0)=NaN;
% 
% %17
% 
% clf
% figure(1)
% 
% pcolor(XC,YC,Thetaappend);
% shading interp
% colormap jet
% hold on
% 
% k=dsearchn(xx,0); %%find lon coordinate corresponding to 0 lon
% plot(XC(243,188:end),YC(243,188:end),'k','linewidth',3);
% hold on
% plot(XC(193,start:start+17),YC(193,start:start+17),'k','linewidth',3);
% hold on
% m = dsearchn(xx,SR4_longitude_89);
% o =  dsearchn(yy',SR4_latitude_89);
% hold on
% plot(xx(m),yy(o),'k','linewidth',3)
% 




