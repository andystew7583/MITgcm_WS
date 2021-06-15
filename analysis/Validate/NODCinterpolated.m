%%%
%%% NODCinterpolated.m
%%%
%%% Validates model experiments against NODC sea ice concentrations.
%%%



%%%%%%%%%%%%%%%%%%
%%%% NODC data

%%% NODC grid dimensions
Nxi = 790;
Neta = 830;

%%% transform our lat/lon meshgrids


%%%%% load datadirectory
NODCiceconc = NaN(Nxi,Neta,12,30);
n=0;
for years = 1978:2007
  
  datadir = fullfile(['/data3/MITgcm_WS/analysis/Validate/NODCseaice/0-data/PolarStereographic/' num2str(years)]);

  for months = 1:12

    NODCiceconc_tmp = NaN(Nxi,Neta,31);

    for days = 1:31

      %%% Print a message to keep track of progress
      ['year = ',num2str(years),', month = ',num2str(months),', day = ',num2str(days)]

      %%% Velocity data file
      filename = fullfile(datadir,['ice_conc_sh_polstere-100_reproc_' num2str(years) num2str(sprintf('%02d',months)) num2str(sprintf('%02d', days)) num2str(1200),'.nc']);
      if (exist(filename)==2)
          lon = ncread(filename, 'lon');
          lat = ncread(filename, 'lat');
          NODCiceconc_tmp(:,:,days) = ncread(filename,'ice_conc');            
      end

    end
    
    %%% Average over available days in each month
    NODCiceconc(:,:,months,years-1977) = nanmean(NODCiceconc_tmp,3);
    
  end
  
end

%%% Save to output file
save(['NODCiceconc.mat'],'NODCiceconc','lon','lat');

