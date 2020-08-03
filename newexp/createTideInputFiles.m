%%%
%%% createTideInputFiles.m
%%%
%%% Creates tide OBCS files.
%%%

%%% Load model grid 
defineGrid;


%%%%%%% Tide BC code needs model grid files to prescribe BCs
addpath ../newexp_utils

data = XMC';
writeDataset(data,fullfile(inputconfigdir,'XMC.bin'),ieee,prec);
clear data;

data = YMC';
writeDataset(data,fullfile(inputconfigdir,'YMC.bin'),ieee,prec);
clear data2;

data = XUMG';
writeDataset(data,fullfile(inputconfigdir,'XUMG.bin'),ieee,prec);
clear data;

data = YUMG';
writeDataset(data,fullfile(inputconfigdir,'YUMG.bin'),ieee,prec);
clear data;

data = zzf;
writeDataset(data,fullfile(inputconfigdir,'zzf.bin'),ieee,prec);
clear data;

data = AngleCS;
writeDataset(data,fullfile(inputconfigdir,'AngleCS.bin'),ieee,prec);
clear data;

data = AngleSN;
writeDataset(data,fullfile(inputconfigdir,'AngleSN.bin'),ieee,prec);
clear data;



%%% Generate tidal BCs
current_dir = pwd();
cd('tides/tmd_toolbox');
creatingtides;
cd(current_dir);



%%% Copy newly-generated MITgcm input files to the relevant input file
%%% configuration directory
copyfile(fullfile('tides/tmd_toolbox','OBEamFile'),fullfile(inputconfigdir,'OBEamFile.bin'));
copyfile(fullfile('tides/tmd_toolbox','OBEphFile'),fullfile(inputconfigdir,'OBEphFile.bin'));
copyfile(fullfile('tides/tmd_toolbox','OBNamFile'),fullfile(inputconfigdir,'OBNamFile.bin'));
copyfile(fullfile('tides/tmd_toolbox','OBNphFile'),fullfile(inputconfigdir,'OBNphFile.bin'));
