%%%
%%% densDiff
%%%
%%% Computes the function E defined by Jackett and McDougall (1997), which
%%% measures the density difference between a parcel and a point on a
%%% reference cast if they were moved to their midpoint pressure. Used to
%%% find the point on a cast that is neutral to a given parcel of water.
%%%
function E = densDiff (ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,kmin,kmax,pp_ref)

  %%% Parcel and reference point on cast will be referred to their midpoint
  %%% pressure
  pp_mean = 0.5 * (pp_ref+pp_parcel);
  
  %%% Determine reference potential temperature and salinity at the point
  %%% pp_ref on the cast via linear interpolation
  ss_ref = interpCast(ss_cast,pp_cast,kmin,kmax,pp_ref);
  pt_ref = interpCast(pt_cast,pp_cast,kmin,kmax,pp_ref);
  
  %%% Compute density difference between parcel and referred point on cast,
  %%% if moved to pressure mid-way between the two
  E = densjmd95(ss_parcel,pt_parcel,pp_mean) - densjmd95(ss_ref,pt_ref,pp_mean);
  
  %%% Alternative definition from Stewart and Thompson (2015)
%   rho0 = 1000;  
%   ss_mean = 0.5 * (ss_ref+ss_parcel);
%   pt_mean = 0.5 * (pt_ref+pt_parcel);
%   [alpha_mean,beta_mean] = calcAlphaBeta(ss_mean,pt_mean,pp_mean);
%   E = rho0 * ( beta_mean*(ss_parcel-ss_ref) - alpha_mean*(pt_parcel-pt_ref) );

end