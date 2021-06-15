
%%%%%%%%%% Reading in depth, temp, salinity





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reading A12 data%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Cruise1 = fullfile('/data1/MITgcm_WS/data/A12/ANT_1992_X4_Fis.nc');


% create a cell array of variables to load
CRUISE1_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(CRUISE1_load)
    % extract the jth variable (type = string)
    var1 = CRUISE1_load{j};

    % use dynamic field name to add this to the structure
    Cruise_A121.(var1) = ncread(Cruise1,var1);
    
    % convert from single to double
    if isa(Cruise_A121.(var1),'single')
        Cruise_A121.(var1) = double(Cruise_A121.(var1));
    end
end





Cruise2 = fullfile('/data1/MITgcm_WS/data/A12/ANT_1996_XIII4_Fis.nc');


% create a cell array of variables to load
CRUISE2_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(CRUISE2_load)
    % extract the jth variable (type = string)
    var2 = CRUISE2_load{j};

    % use dynamic field name to add this to the structure
    Cruise_A122.(var2) = ncread(Cruise2,var2);

    % convert from single to double
    if isa(Cruise_A122.(var2),'single')
        Cruise_A122.(var2) = double(Cruise_A122.(var2));
    end
end



Cruise3 = fullfile('/data1/MITgcm_WS/data/A12/ANT_1998_XV4_Fis.nc');

% create a cell array of variables to load
CRUISE3_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(CRUISE3_load)
    % extract the jth variable (type = string)
    var3 = CRUISE3_load{j};

    % use dynamic field name to add this to the structure
    Cruise_A123.(var3) = ncread(Cruise3,var3);

    % convert from single to double
    if isa(Cruise_A123.(var3),'single')
        Cruise_A123.(var3) = double(Cruise_A123.(var3));
    end
end


Cruise4 = fullfile('/data1/MITgcm_WS/data/A12/ANT_1999_XVI2_Fis.nc');

% create a cell array of variables to load
CRUISE4_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(CRUISE4_load)
    % extract the jth variable (type = string)
    var4 = CRUISE4_load{j};

    % use dynamic field name to add this to the structure
    Cruise_A124.(var4) = ncread(Cruise4,var4);

    % convert from single to double
    if isa(Cruise_A124.(var4),'single')
        Cruise_A124.(var4) = double(Cruise_A124.(var4));
    end
end



Cruise5 = fullfile('/data1/MITgcm_WS/data/A12/ANT_2000_2001_XVIII3_Fis.nc');

% create a cell array of variables to load
CRUISE5_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(CRUISE5_load)
    % extract the jth variable (type = string)
    var5 = CRUISE5_load{j};

    % use dynamic field name to add this to the structure
    Cruise_A125.(var5) = ncread(Cruise5,var5);

    % convert from single to double
    if isa(Cruise_A125.(var5),'single')
        Cruise_A125.(var5) = double(Cruise_A125.(var5));
    end
end


Cruise6 = fullfile('/data1/MITgcm_WS/data/A12/ANT_2002_2003_XX2_Fis.nc');

% create a cell array of variables to load
CRUISE6_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(CRUISE6_load)
    % extract the jth variable (type = string)
    var6 = CRUISE6_load{j};

    % use dynamic field name to add this to the structure
    Cruise_A126.(var6) = ncread(Cruise6,var6);

    % convert from single to double
    if isa(Cruise_A126.(var6),'single')
        Cruise_A126.(var6) = double(Cruise_A126.(var6));
    end
end





Cruise7 = fullfile('/data1/MITgcm_WS/data/A12/ANT_2005_XXII3_Fis.nc');

% create a cell array of variables to load
CRUISE7_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(CRUISE7_load)
    % extract the jth variable (type = string)
    var7 = CRUISE7_load{j};

    % use dynamic field name to add this to the structure
    Cruise_A127.(var7) = ncread(Cruise7,var7);

    % convert from single to double
    if isa(Cruise_A127.(var7),'single')
        Cruise_A127.(var7) = double(Cruise_A127.(var7));
    end
end




Cruise8 = fullfile('/data1/MITgcm_WS/data/A12/ANT_2008_XXIV3_Fis.nc');
 
% create a cell array of variables to load
CRUISE8_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(CRUISE8_load)
    % extract the jth variable (type = string)
    var8 = CRUISE8_load{j};

    % use dynamic field name to add this to the structure
    Cruise_A128.(var8) = ncread(Cruise8,var8);

    % convert from single to double
    if isa(Cruise_A128.(var8),'single')
        Cruise_A128.(var8) = double(Cruise_A128.(var8));
    end
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% SR4 $$$$$$$%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SR4Cruise1= fullfile('/data1/MITgcm_WS/data/SR4/AJAX_Leg2.nc');




% create a cell array of variables to load
SR4CRUISE1_load = {'longitude','latitude','var2','var3','var1'};

% loop over the variables
for j=1:numel(SR4CRUISE1_load)
    % extract the jth variable (type = string)
    SR4var1 = SR4CRUISE1_load{j};

    % use dynamic field name to add this to the structure
    Cruise_SR41.(SR4var1) = ncread(SR4Cruise1,SR4var1);

    % convert from single to double
    if isa(Cruise_SR41.(SR4var1),'single')
        Cruise_SR41.(SR4var1) = double(Cruise_SR41.(SR4var1));
    end
end


SR4Cruise2= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1993_X7_Fis.nc');

% create a cell array of variables to load
SR4CRUISE2_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(SR4CRUISE2_load)
    % extract the jth variable (type = string)
    SR4var2 = SR4CRUISE2_load{j};

    % use dynamic field name to add this to the structure
    Cruise_SR42.(SR4var2) = ncread(SR4Cruise2,SR4var2);

    % convert from single to double
    if isa(Cruise_SR42.(SR4var2),'single')
        Cruise_SR42.(SR4var2) = double(Cruise_SR42.(SR4var2));
    end
end



SR4Cruise3= fullfile('/data1/MITgcm_WS/data/SR4/ANT_2005_XXII3_Fis.nc');
% create a cell array of variables to load
SR4CRUISE3_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(SR4CRUISE3_load)
    % extract the jth variable (type = string)
    SR4var3 = SR4CRUISE3_load{j};

    % use dynamic field name to add this to the structure
    Cruise_SR43.(SR4var3) = ncread(SR4Cruise3,SR4var3);

    % convert from single to double
    if isa(Cruise_SR43.(SR4var3),'single')
        Cruise_SR43.(SR4var3) = double(Cruise_SR43.(SR4var3));
    end
end



SR4Cruise4= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1989_VIII2_Fis.nc');  
% create a cell array of variables to load
SR4CRUISE4_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(SR4CRUISE4_load)
    % extract the jth variable (type = string)
    SR4var4 = SR4CRUISE4_load{j};

    % use dynamic field name to add this to the structure
    Cruise_SR44.(SR4var4) = ncread(SR4Cruise4,SR4var4);

    % convert from single to double
    if isa(Cruise_SR44.(SR4var4),'single')
        Cruise_SR44.(SR4var4) = double(Cruise_SR44.(SR4var4));
    end
end



SR4Cruise5= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1996_XIII4_Fis.nc');
% create a cell array of variables to load
SR4CRUISE5_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(SR4CRUISE5_load)
    % extract the jth variable (type = string)
    SR4var5 = SR4CRUISE5_load{j};

    % use dynamic field name to add this to the structure
    Cruise_SR45.(SR4var5) = ncread(SR4Cruise5,SR4var5);

    % convert from single to double
    if isa(Cruise_SR45.(SR4var5),'single')
        Cruise_SR45.(SR4var5) = double(Cruise_SR45.(SR4var5));
    end
end



SR4Cruise6= fullfile('/data1/MITgcm_WS/data/SR4/ANT_2008_XXIV3_Fis.nc');
% create a cell array of variables to load
SR4CRUISE6_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(SR4CRUISE6_load)
    % extract the jth variable (type = string)
    SR4var6 = SR4CRUISE6_load{j};

    % use dynamic field name to add this to the structure
    Cruise_SR46.(SR4var6) = ncread(SR4Cruise6,SR4var6);

    % convert from single to double
    if isa(Cruise_SR46.(SR4var6),'single')
        Cruise_SR46.(SR4var6) = double(Cruise_SR46.(SR4var6));
    end
end



SR4Cruise7= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1990_IX2_Fis.nc');    
% create a cell array of variables to load
SR4CRUISE7_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(SR4CRUISE7_load)
    % extract the jth variable (type = string)
    SR4var7 = SR4CRUISE7_load{j};

    % use dynamic field name to add this to the structure
    Cruise_SR47.(SR4var7) = ncread(SR4Cruise7,SR4var7);

    % convert from single to double
    if isa(Cruise_SR47.(SR4var7),'single')
        Cruise_SR47.(SR4var7) = double(Cruise_SR47.(SR4var7));
    end
end


SR4Cruise8= fullfile('/data1/MITgcm_WS/data/SR4/ANT_1998_XV4_Fis.nc'); 
% create a cell array of variables to load
SR4CRUISE8_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(SR4CRUISE8_load)
    % extract the jth variable (type = string)
    SR4var8 = SR4CRUISE8_load{j};

    % use dynamic field name to add this to the structure
    Cruise_SR48.(SR4var8) = ncread(SR4Cruise8,SR4var8);

    % convert from single to double
    if isa(Cruise_SR48.(SR4var8),'single')
        Cruise_SR48.(SR4var8) = double(Cruise_SR48.(SR4var8));
    end
end


SR4Cruise9= fullfile('/data1/MITgcm_WS/data/SR4/ANT_2011_XXVII2_Fis.nc');
% create a cell array of variables to load
SR4CRUISE9_load = {'Longitude','Latitude','VAR_2','VAR_3','VAR_1'};

% loop over the variables
for j=1:numel(SR4CRUISE9_load)
    % extract the jth variable (type = string)
    SR4var9 = SR4CRUISE9_load{j};

    % use dynamic field name to add this to the structure
    Cruise_SR49.(SR4var9) = ncread(SR4Cruise9,SR4var9);

    % convert from single to double
    if isa(Cruise_SR49.(SR4var9),'single')
        Cruise_SR49.(SR4var9) = double(Cruise_SR49.(SR4var9));
    end
end



save Cruisedata.mat Cruise_A121 Cruise_A122 Cruise_A123 Cruise_A124 Cruise_A125 Cruise_A126 Cruise_A127 Cruise_A128 ...
    Cruise_SR41 Cruise_SR42 Cruise_SR43 Cruise_SR44 Cruise_SR45 Cruise_SR46 ...
    Cruise_SR47 Cruise_SR48








%%%%%%%%%%%%%%%%%%%%%%%% Trying to extract relevant casts
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cruise 1    


% LONCruise1 = ncread(Cruise1,'Longitude');
% LONCruise1(1:58)=[];
% LONCruise1(19:end)=[];
% 
% 
% LATCruise1 = ncread(Cruise1,'Latitude');
% LATCruise1(1:58)=[];
% LATCruise1(19:end)=[];
% 
% 
% TempCruise1 = ncread(Cruise1,'VAR_2');
% TempCruise1(:,1:58)=[];
% TempCruise1(:,19:end)=[];
% TempCruise1=TempCruise1(1,:);
% 
% TempC1 =meshgrid(TempCruise1,LONCruise1,LATCruise1);
% 
% 
% 
% SalCruise1 = ncread(Cruise1,'VAR_3');
% SalCruise1(:,1:58)=[];
% SalCruise1(:,19:end)=[];
% SalCruise1=SalCruise1(1,:);
% 
% SaltC1 =meshgrid(SalCruise1,LONCruise1,LATCruise1);
% 
% 
% 
% DepthCruise1 = ncread(Cruise1,'VAR_1');
% DepthCruise1(:,1:58)=[];
% DepthCruise1(:,19:end)=[];
% DepthCruise1=DepthCruise1(1,:);
% 
% DepthC1 =meshgrid(DepthCruise1,LONCruise1,LATCruise1);
% 
% 
% 
% 
% %%%% Cruise 2
% TempCruise2 = ncread(Cruise2,'VAR_2');
% TempCruise2(1:69)=[];
% TempCruise2(19:end)=[];
% 
% SalCruise2 = ncread(Cruise2,'VAR_3');
% SalCruise2(1:69)=[];
% SalCruise2(19:end)=[];
% 
% DepthCruise2 = ncread(Cruise2,'VAR_1');
% DepthCruise2(1:69)=[];
% DepthCruise2(19:end)=[];
% 
% LONCruise2 = ncread(Cruise2,'Longitude');
% LONCruise2(1:69) = [];
% LONCruise2(19:end) = [];
% 
% LATCruise2 = ncread(Cruise2,'Latitude');
% LATCruise2(1:69)=[];
% LATCruise2(19:end)=[];
% 
% 
% %%%% Cruise 3
% TempCruise3 = ncread(Cruise3,'VAR_2');
% TempCruise3(1:86)=[];
% TempCruise3(14:end)=[];
% 
% SalCruise3 = ncread(Cruise3,'VAR_3');
% SalCruise3(1:86)=[];
% SalCruise3(14:end)=[];
% 
% 
% DepthCruise3 = ncread(Cruise3,'VAR_1');
% DepthCruise3(1:86)=[];
% DepthCruise3(14:end)=[];
% 
% 
% 
% LONCruise3 = ncread(Cruise3,'Longitude');
% LONCruise3(1:86)=[];
% LONCruise3(14:end)=[];
% 
% LATCruise3 = ncread(Cruise3,'Latitude');
% LATCruise3(1:86)=[];
% LATCruise3(14:end)=[];
% 
% 
% %%%% Cruise 4
% TempCruise4 = ncread(Cruise4,'VAR_2');
% SalCruise4 = ncread(Cruise4,'VAR_3');
% DepthCruise4 = ncread(Cruise4,'VAR_1');
% LONCruise4 = ncread(Cruise4,'Longitude');
% LATCruise4 = ncread(Cruise4,'Latitude');
% 
% %%%%% Cruise 5    
% TempCruise5 = ncread(Cruise5,'VAR_2');
% SalCruise5 = ncread(Cruise5,'VAR_3');
% DepthCruise5 = ncread(Cruise5,'VAR_1');
% LONCruise5 = ncread(Cruise5,'Longitude');
% LATCruise5 = ncread(Cruise5,'Latitude');
% 
%     
% %%%% Cruise 6
% TempCruise6 = ncread(Cruise6,'VAR_2');
% SalCruise6 = ncread(Cruise6,'VAR_3');
% DepthCruise6 = ncread(Cruise6,'VAR_1');
% LONCruise6 = ncread(Cruise6,'Longitude');
% LATCruise6 = ncread(Cruise6,'Latitude');
% 
% %%%% Cruise 7
% TempCruise7 = ncread(Cruise7,'VAR_2');
% SalCruise7 = ncread(Cruise7,'VAR_3');
% DepthCruise7 = ncread(Cruise7,'VAR_1');
% LONCruise7 = ncread(Cruise7,'Longitude');
% LATCruise7 = ncread(Cruise7,'Latitude');
% 
% %%%% Cruise 8
% TempCruise8 = ncread(Cruise8,'VAR_2');
% SalCruise8 = ncread(Cruise8,'VAR_3');
% DepthCruise8 = ncread(Cruise8,'VAR_1');
% LONCruise8 = ncread(Cruise8,'Longitude');
% LATCruise8 = ncread(Cruise8,'Latitude');






%%%%%%%%%%  TO LOAD A12 FROM ASII FILE %%%%%%%%%%%%
% 
% A12 = dlmread('./ewoce/data/whp/ctd/cfgAtlanticSections/A12.sec','',1,0);
% 
% S04A = dlmread('./ewoce/data/whp/ctd/cfgSouthOceanSections/S04A.sec','',1,0);