%%%
%%% interpCast
%%%
%%% Convenience function to interpolate vertical cast data onto a specified
%%% pressure level. Currently configured to perform simple linear
%%% interpolation with nearest-neighbor extrapolation.
%%%
function ff_interp = interpCast (ff_cast,pp_cast,kmin,kmax,pp_interp)

 if (pp_interp<=pp_cast(kmin))
    ff_interp = ff_cast(kmin);    
  else
    if (pp_interp>=pp_cast(kmax))
      ff_interp = ff_cast(kmax);      
    else
      knext = find(pp_cast>pp_interp,1,'first');
      kprev = knext-1;
      wnext = (pp_interp-pp_cast(kprev)) ./ (pp_cast(knext)-pp_cast(kprev));
      wprev = 1-wnext;
      ff_interp = wprev.*ff_cast(kprev) + wnext.*ff_cast(knext);      
    end
 end
  
end