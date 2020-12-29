%%%
%%% defineMOCgrid.m
%%%
%%% Creates an alternative coordinate system in which to compute the
%%% overturning circulation.
%%%
function ETA = defineMOCgrid (XC,YC,icedraft,bathy,deform_cavity)

  %%% Coordinates marking ends of FRIS front
  lat1 = -74.8;
  lon1 = -61; 
  lat2 = -78.35;
  lon2 = -36.7;

  %%% Rotation angle for coordinates
  phic = atan((lat2-lat1)/(lon2-lon1));

  %%% Rotated coordinate grid
  ETA0 = (YC-lat1)*cos(phic) - (XC-lon1)*sin(phic);
  sin(phic)/cos(phic)
  sin(phic)
  cos(phic)

  %%% If needed, modify the coordinates within the FRIS cavity
  if (deform_cavity)

    %%% Transition value of eta: coordinates look like ETA0 to the ''north''
    %%% of this value
    eta_t = 4; 

    %%% Cavity coordinate system: distance from NW FRIS
%     ETA1 = sqrt(((XC-lon1).*cosd(YC)).^2 + .25*((YC-lat1)).^2);    
    ETA1 = sqrt(((XC-lon1).*cosd(YC)).^2 + 0.25*((YC-lat1)).^2);

    %%% Normalize to ensure eta remains smaller than the transition eta
    %%% throughout the cavity
%     etafrac = 1;
    etafrac = .75;
    ETA1 = ETA1 * etafrac*eta_t / max(ETA1((icedraft-bathy>0) & (ETA0<eta_t)));

    %%% Sinusoidal weights for transition between cavity and open ocean grids
    wt0 = sin(pi/2*ETA0/eta_t);
    wt0(ETA0<0) = 0;
    wt0(ETA0>eta_t) = 1;
    wt1 = cos(pi/2*ETA0/eta_t);
    wt1(ETA0<0) = 1;
    wt1(ETA0>eta_t) = 0;

    %%% Compute the combined eta grid
    ETA = wt1 .* ETA1 + wt0 .* ETA0;

  %%% If no cavity coordinates are needed then just return the rotated grid
  else

    ETA = ETA0;

  end
  
end
