%%%%%%%%%% Using NISDC data instead of SOSE for sea ice area
run defineGrid.m
load ../analysis/Validate/Ice_a.mat
addpath ../newexp_utils/


%%% find the monthly mean of the NODC data

Ice_monthly = nanmean(Ice_Seasonal,4);

Ice_monthly = Ice_monthly/100;

%%% finding the ice concentration on the EB

EB_Ice = squeeze(Ice_monthly(end,:,:));
NB_Ice = squeeze(Ice_monthly(:,end,:));



data = EB_Ice;
writeDataset(data,fullfile(inputfolder,OBEaFile),ieee,prec);
clear data

data = NB_Ice;
writeDataset(data,fullfile(inputfolder,OBNaFile),ieee,prec);
clear data

