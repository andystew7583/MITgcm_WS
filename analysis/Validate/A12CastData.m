%%%%%%%%%%v%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Testing Boundaries with Cast Data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Linearly interpolate Cast data to model grid in vertical to our grid
%%%%%%%%%%%%%



%%%%%% A12 -> Eastern Boundary
%%%%%% load data
load Cruisedata.mat



%%%%%%% Fitting A12 Lon/Lat to eastern boundary

addpath 'A12castplots';

datadir ='/data1/MITgcm_WS/data/A12/';


%%%%%%%%%%% depths


depth1 = ncread(fullfile(datadir,'ANT_1992_X4_Fis.nc'),'VAR_1');  
depth2 = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'VAR_1');
depth3 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_1');
% depth4 = ncread(fullfile(datadir,'ANT_1999_XVI2_Fis.nc'),'VAR_1');
depth4 = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'VAR_1');
depth6 = ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'VAR_1');
depth7 = ncread(fullfile(datadir,'ANT_2005_XXII3_Fis.nc'),'VAR_1');
depth8 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_1');



%%%%%% Calculate Temp of Sections

Temp1= ncread(fullfile(datadir,'ANT_1992_X4_Fis.nc'),'VAR_2');
Temp2= ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'VAR_2');
Temp3 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_2');
% Temp4 = ncread(fullfile(datadir,'ANT_1999_XVI2_Fis.nc'),'VAR_2');
Temp4 = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'VAR_2');
Temp6 = ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'VAR_2');
Temp7 = ncread(fullfile(datadir,'ANT_2005_XXII3_Fis.nc'),'VAR_2');
Temp8 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_2');

%%%%%% Calculate Temp of Sections

Salt1= ncread(fullfile(datadir,'ANT_1992_X4_Fis.nc'),'VAR_3');
Salt2= ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'VAR_3');
Salt3 = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'VAR_3');
% Salt4 = ncread(fullfile(datadir,'ANT_1999_XVI2_Fis.nc'),'VAR_3');
Salt4 = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'VAR_3');
Salt6 = ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'VAR_3');
Salt7 = ncread(fullfile(datadir,'ANT_2005_XXII3_Fis.nc'),'VAR_3');
Salt8 = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'VAR_3');



%%%%%%%%%%% Cruise 1

A12Lon1 = Cruise_A121.Longitude;
A12Lon1(1:34) = [];
A12Lon1(43:end) = [];


A12Lat1 = Cruise_A121.Latitude;
A12Lat1(1:34) = [];
A12Lat1(43:end) = [];


Temp1(:,1:34) = [];
Temp1(:,43:end) = [];

Salt1(:,1:34) = [];
Salt1(:,43:end) = [];

depth1 = depth1(:,1);
[lat1,d1]=meshgrid(A12Lat1,depth1);


%%%%%%% plot
figure(1)
pcolor(lat1,d1,Temp1),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(50)
caxis([-3 4]);
title('A12 Temperature 1992')
print(figure(1),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A121992Temp.png']);

figure(2)
pcolor(lat1,d1,Salt1),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(30)
title('A12 Salinity 1992')
print(figure(2),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A121992Salt.png']);

%%%%%%%%%%%% Section 2 

A12_2_Lon = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'Longitude');
A12_2_Lat = ncread(fullfile(datadir,'ANT_1996_XIII4_Fis.nc'),'Latitude');

A12_2_Lon(1:45) = [];
A12_2_Lon(43:end) = [];

A12_2_Lat(1:45) = [];
A12_2_Lat(43:end) = [];

Temp2(:,1:45) = [];
Temp2(:,43:end) = [];

Salt2(:,1:45) = [];
Salt2(:,43:end) = [];

depth2 = depth2(:,1);

[lat2,d2]=meshgrid(A12_2_Lat,depth2);



%%%%%%% plot
figure(3)
pcolor(lat2,d2,Temp2),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(50)
caxis([-2 4]);
title('A12 Temperature 1996')
print(figure(3),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A121996Temp.png']);

figure(4)
pcolor(lat2,d2,Salt2),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(30)
title('A12 Salinity 1996')
print(figure(3),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A121996Salt.png']);



%%%%%%%%%%%%%%%%%%%%%%% Section 3 

A12_3_Lon = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'Longitude');
A12_3_Lat = ncread(fullfile(datadir,'ANT_1998_XV4_Fis.nc'),'Latitude');

A12_3_Lon(1:86) = [];
A12_3_Lon(58:end) = [];

A12_3_Lat(1:86) = [];
A12_3_Lat(58:end) = [];

Temp3(:,1:86) = [];
Temp3(:,58:end) = [];

Salt3(:,1:86) = [];
Salt3(:,58:end) = [];

depth3 = depth3(:,1);

[lat3,d3]=meshgrid(A12_3_Lat,depth3);


%%%%%%% plot
figure(5)
pcolor(lat3,d3,Temp3),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(50)
caxis([-3 4]);
title('A12 Temperature 1998')
print(figure(5),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A121998Temp.png']);

figure(6)
pcolor(lat3,d3,Salt3),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(30)
title('A12 Salinity 1998')
print(figure(6),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A121998Salt.png']);



%%%%%%%%%%%%%%%%%%%%%%% Section 4

A12_4_Lon = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'Longitude');
A12_4_Lat = ncread(fullfile(datadir,'ANT_2000_2001_XVIII3_Fis.nc'),'Latitude');

A12_4_Lon(50:end) = [];

A12_4_Lat(50:end) = [];

Temp4(:,50:end) = [];
Salt4(:,50:end) = [];

depth4 = depth4(:,1);

[lat4,d4]=meshgrid(A12_4_Lat,depth4);


%%%%%%% plot
figure(7)
pcolor(lat4,d4,Temp4),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(50)
caxis([-3 4]);
title('A12 Temperature 2000')
print(figure(7),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A122000Temp.png']);

figure(8)
pcolor(lat4,d4,Salt4),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(30)
title('A12 Salinity 2000')
print(figure(8),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A122000Salt.png']);




%%%%%%%%%%%%%%%%%%%%%%% Section 5

A12_6_Lon = ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'Longitude');
A12_6_Lat = ncread(fullfile(datadir,'ANT_2002_2003_XX2_Fis.nc'),'Latitude');

A12_6_Lon(1:6) = [];
A12_6_Lon(39:end) = [];

A12_6_Lat(1:6) = [];
A12_6_Lat(39:end) = [];

Temp6(:,1:6) = [];
Temp6(:,39:end) = [];

Salt6(:,1:6) = [];
Salt6(:,39:end) = [];


depth6 = depth6(:,1);

[lat6,d6]=meshgrid(A12_6_Lat,depth6);


%%%%%%% plot
figure(9)
pcolor(lat6,d6,Temp6),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(50)
caxis([-3 4]);
title('A12 Temperature 2002')
print(figure(9),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A122002Temp.png']);

figure(10)
pcolor(lat6,d6,Salt6),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(30)
caxis([33.5 34.5])
title('A12 Salinity 2002')
print(figure(10),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A122002Salt.png']);



%%%%%%%%%%%%%%%%%%%%%%% Section 6

A12_7_Lon = ncread(fullfile(datadir,'ANT_2005_XXII3_Fis.nc'),'Longitude');
A12_7_Lat = ncread(fullfile(datadir,'ANT_2005_XXII3_Fis.nc'),'Latitude');

% A12_7_Lon(1:6) = [];
A12_7_Lon(44:end) = [];

% A12_7_Lat(1:6) = [];
A12_7_Lat(44:end) = [];

% Temp7(1:6) = [];
Temp7(:,44:end) = [];
Salt7(:,44:end) = [];

depth7 = depth7(:,1);

[lat7,d7]=meshgrid(A12_7_Lat,depth7);


%%%%%%% plot
figure(11)
pcolor(lat7,d7,Temp7),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(50)
caxis([-3 4]);
title('A12 Temperature 2005')
print(figure(11),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A122005Temp.png']);

figure(12)
pcolor(lat7,d7,Salt7),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(30)
title('A12 Salinity 2005')
print(figure(12),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A122005Salt.png']);




%%%%%%%%%%%%%%%%%%%%%%% Section 7

A12_8_Lon = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'Longitude');
A12_8_Lat = ncread(fullfile(datadir,'ANT_2008_XXIV3_Fis.nc'),'Latitude');

A12_8_Lon(1:14) = [];
A12_8_Lon(76:end) = [];

A12_8_Lat(1:14) = [];
A12_8_Lat(76:end) = [];

Temp8(:,1:14) = [];
Temp8(:,76:end) = [];

Salt8(:,1:14) = [];
Salt8(:,76:end) = [];


depth8 = depth8(:,1);

[lat8,d8]=meshgrid(A12_8_Lat,depth8);


%%%%%%% plot
figure(13)
pcolor(lat8,d8,Temp8),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(50)
caxis([-3 4]);
title('A12 Temperature 2008')
print(figure(13),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A122008Temp.png']);

figure(14)
pcolor(lat8,d8,Salt8),shading interp;
set (gca,'YDir','reverse');
colorbar
colormap jet(30)
caxis([33 35])
title('A12 Salinity 2008')
print(figure(14),'-dpng',['/data1/MITgcm_WS/analysis/ValidateSOSE/A12castplots/A122008Salt.png']);


