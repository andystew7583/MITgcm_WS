%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Validating Boundary Conditions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing Boundary Conditions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% First have to load experiment data

run ../setExpname;
run ../loadexp;
%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(7));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);


tmin = 8*86400*365;
tmax = 9*86400*365;
Theta = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);
Salt = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);



%%%%% loading SR4 data
run ../../newexp/defineGrid.m
datadir ='/data1/MITgcm_WS/data/SR4/';

%%%% SR04 latitude and longitude

SR4_latitude_89 = ncread(fullfile(datadir,'ANT_1989_VIII2_Fis.nc'),'Latitude');
SR4_longitude_89 = ncread(fullfile(datadir,'ANT_1989_VIII2_Fis.nc'),'Longitude');

%%%%%% extracting SR4 temperature


N_Stations = 64;

%%%%%% FIRST PUTTING SR4 DATA ON OUR VERTICAL GRID


tempsr4_89 = ncread(fullfile(datadir,'ANT_1989_VIII2_Fis.nc'),'VAR_2');
depthsr4_89 = ncread(fullfile(datadir,'ANT_1989_VIII2_Fis.nc'),'VAR_1');

temp_newsr489 = zeros(length(zz),size(depthsr4_89,2));
for n=1:size(depthsr4_89,2)
  idx = find(~isnan(depthsr4_89(:,n)));
  temp_newsr489(:,n) = interp1(depthsr4_89(idx,n),tempsr4_89(idx,n),squeeze(-zz),'linear');
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




%%%%%%%%% 1996 cruise



SR4_96_Lon = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'Longitude');
SR4_96_Lat = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'Latitude');

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





%%%%%%%%%%%%%%%%%%%%%% 2008 cruise


depthsr4_08 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_1');
tempsr4_08 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_2');

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






%%%%%%%%%%%%% 1990
% 
% SR4Cruise_90= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1990_IX2_Fis.nc');   
% depthsr4_90 = ncread(fullfile(datadir,'ANT_1990_IX2_Fis.nc'),'VAR_1');
% tempsr4_90 = ncread(fullfile(datadir,'ANT_1990_IX2_Fis.nc'),'VAR_2');
% SR4_90_Lon = ncread(fullfile(datadir,'ANT_1990_IX2_Fis.nc'),'Longitude');
% SR4_90_Lat = ncread(fullfile(datadir,'ANT_1990_IX2_Fis.nc'),'Latitude');
% 
% 
% temp_90 = zeros(length(zz),size(depthsr4_90,2));
% for n=1:size(depthsr4_90,2)
%   idx = find(~isnan(depthsr4_90(:,n)));
%   temp_90(:,n) = interp1(depthsr4_90(idx,n),tempsr4_90(idx,n),squeeze(-zz),'linear');
% end
% 
% Theta_90LO = zeros(length(zz),N_Stations);
% for n=1:length(zz)
%    sr4_temp_90 = temp_90(n,:);
%    Theta_90LO(n,:) = interp1(SR4_90_Lon,sr4_temp_90,SR4_longitude_89,'linear');
% end
% 

%%%%%% 1998

SR4Cruise_98= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1998_XV4_Fis.nc'); 
SR4_98_Lon = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'Longitude');
SR4_98_Lat = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'Latitude');
depthsr4_98 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_1');
tempsr4_98 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_2');




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







%%%%%%% 2011

SR4Cruise9= fullfile('/data1/MITgcm_WS/data/SR4/ANT_2011_XXVII2_Fis.nc');

SR4_11_depthsr49 = ncread(fullfile(datadir,'ANT_2011_XXVII2_Fis.nc'),'VAR_1');
SR4_11_Lon = ncread(fullfile(datadir,'ANT_2011_XXVII2_Fis.nc'),'Longitude');
SR4_11_Lat = ncread(fullfile(datadir,'ANT_2011_XXVII2_Fis.nc'),'Latitude');



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


















%%%%%%% trying to put OUR temperature results on SR4
%%%%%%% Using 1989 cruise as defined SR4 lat/lon points


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
        Theta_MOD(k,:,n) = interp1(ymc,sr4_temp_lat,SR4_latitude_89,'linear');
    end
end

%%% make a meshgrid

[SR4_LO,zzm]=meshgrid(SR4_longitude_89,zz);

Theta_MOD = squeeze(Theta_MOD(:,18,:));

%%%%figure plotting model output on SR4 track
pcolor(SR4_LO',zzm',Theta_MOD), shading interp
colormap jet(100)
colorbar
caxis([-3 1])















% Weddell_LO = XMC';
% Weddell_LA = YMC';
% 
% temp = scattercasts(temp_newsr4,SR4_LO,SR4_LA,Weddell_LO,Weddell_LA);


%%%%%%% Using delaunayn to see closest longitude/latitude to the sr4 track

% x = dsearchn(xmc,SR4_longitude);
% 
% y = dsearchn(ymc',SR4_latitude);
% 
% 
% n_xmc = xmc(x);
% n_ymc = ymc(y);
% 
% N_XMC = XMC(x,y);
% N_YMC = YMC(x,y);
% 
% 
% [n_x,n_z]=meshgrid(x,zz);
% 
% depthsr4 = inpaint_nans(depthsr4,4);
% d = depthsr4(:,1);
% 
% NStations = 64;
% CastTemp = NaN(Nz,NStations);
% Cast_temp = zeros(size(tempsr4,2));
% for i = 1:size(tempsr4,1)
%     Cast_temp = tempsr4(i,:);
%     CastTemp(:,i) = interp1(depthsr4,Cast_temp,n_z); 
% end
% 
% %%%%% Now I need to interpolate in the horizontal and vertical to my grid
% 
% 






