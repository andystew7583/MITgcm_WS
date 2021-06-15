%%%v
%%% readIters.m
%%%
%%% Reads ans sums all iterations of a specified MITgcm output field
%%% between specified times, and calculates the time average.
%%%
function avg = readIters (exppath,field,dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr)
 
  avg = zeros(Nx,Ny,Nr);
  navg = 0;
  
  %%% Loop through output iterations
  for n=1:length(dumpIters)
    tseconds =  dumpIters(n)*deltaT;
  
    if ((tseconds >= tmin) && (tseconds <= tmax))
    
      fname = fullfile(exppath,'results',field);      
      A = rdmdsWrapper(fname,dumpIters(n));      
      if (isempty(A))
        error(['Could not find ',fname,' iteration ',num2str(dumpIters(n))]);
      end
      avg = avg + A;
      navg = navg + 1;

    end
    
  end
  
  %%% Calculate average
  if (navg > 0)
    avg = avg / navg;
  else
    error('No output files found');
  end

end

