function avgcastdata = castavg (years,years2,roman,depth,ncname)

    %%% Calculate cast average variable
    avgcastdata = 0*depth;    
    navg = 0;
    for n=1:9
      castdata = fullfile(datadir,['ANT_' num2str(years) '_' num2str(years2) roman '_Fis.nc']);                
      if (exist(castdata))
        avgcastdata = avgcastdata + ncread(castdata,ncname);
        navg = navg + 1;
      end
    end
    
    if (navg > 0)
      avgcastdata = avgcastdata / navg;
    end
end

