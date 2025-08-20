%%%
%%% defineMOCgrid.m
%%%
%%% Creates an alternative coordinate system in which to compute the
%%% overturning circulation.
%%%
function ETA = defineMOCgrid (XC,YC,icedraft,bathy,deform_cavity,gl_coord)

  %%% Coordinates marking ends of FRIS front
  lat1 = -74.8;
  lon1 = -61; 
  lat2 = -78.35;
  lon2 = -36.7;

  %%% Rotation angle for coordinates
  phic = atan((lat2-lat1)/(lon2-lon1));

  %%% Rotated coordinate grid
  ETA0 = (YC-lat1)*cos(phic) - (XC-lon1)*sin(phic);  

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

  %%% Hyper-engineered coordinate system that places the grounding lines
  %%% (as close as possible) to an isoline of the coordinate system
  elseif (gl_coord)

    %%% Transition value of eta: coordinates look like ETA0 to the ''north''
    %%% of this value
    eta_t = 1;    

    %%% Value of eta at "mid-cavity" - used for successive interpolation
    %%% toward complex grounding lines
    eta_mc = -5;

    eta_off = 2.5;

    %%% Value of eta at the grounding lines
    eta_gl = -8;     

    %%% Value of eta at the ice front
    eta_if = 0;

    %%% Value of eta at "deep cavity" location - used for successive interpolation
    %%% toward complex grounding lines
    eta_dc = -7.75;


    
    %%% Grid inside the cavity, notionally varies from eta_if to eta_gl
    ETA2 = ETA0;
    ETA2(icedraft == 0) = eta_if;
    ETA2((icedraft<=bathy)) = eta_gl+eta_off;
    ETA2((icedraft > bathy) & (icedraft < 0)) = NaN;
    ETA2((XC > -67) & (YC>-81) & (XC<-42) & (YC<-77.5) & (icedraft<=bathy)) = NaN; %%% Manually remove islands in FRISs
    ETA2((XC<=-67) & (YC>-80.5) & (XC>-72) & (YC<-78.1) & (icedraft<=bathy)) = NaN;
    ETA2 = inpaint_nans(ETA2);
    % idx_known = find(~isnan(ETA2));
    % idx_unknown = find(isnan(ETA2));
    % F = scatteredInterpolant(XC(idx_known),YC(idx_known),ETA2(idx_known));
    % ETA2(idx_unknown) = F(XC(idx_unknown),YC(idx_unknown));

    %%% Grid outside cavity - approaches standard ETA grid at eta_t
    ETA1 = ETA2;
    ETA1(icedraft == 0 | (ETA0>-1)) = NaN;
    % ETA1((ETA2 > eta_if) & (ETA0 < eta_t)) = NaN;    
    ETA1(ETA0>=eta_t) = ETA0(ETA0>=eta_t);
    ETA1(icedraft <= bathy) = NaN;
    ETA1((icedraft<0) & (icedraft > bathy) & (ETA0 < 0)) = eta_if;
    ETA1 = inpaint_nans(ETA1);
    % idx_known = find(~isnan(ETA1));
    % idx_unknown = find(isnan(ETA1));
    % F = scatteredInterpolant(XC(idx_known),YC(idx_known),ETA1(idx_known));
    % ETA1(idx_unknown) = F(XC(idx_unknown),YC(idx_unknown));

    %%% Grid deeper inside cavity - interpolates better toward complex GL
    ETA3 = ETA2;
    ETA3((icedraft<=bathy)) = eta_gl;
    ETA3((icedraft > bathy) & (icedraft < 0)) = NaN;
    ETA3((XC > -67) & (YC>-81) & (XC<-42) & (YC<-77.5) & (icedraft<=bathy)) = NaN;
    ETA3((XC<=-67) & (YC>-80.5) & (XC>-72) & (YC<-78.1) & (icedraft<=bathy)) = NaN;
    ETA3(ETA2>=eta_mc) = ETA2(ETA2>=eta_mc);
    ETA3 = inpaint_nans(ETA3);

    %%% Grid very close to grounding line: it's very difficult to make it
    %%% exactly eta_gl, so we allow for some flexibility by defining an
    %%% external box around the cavity
    ETA4 = ETA3;
    ETA4((icedraft<=bathy)) = NaN;
    ETA4((icedraft > bathy) & (icedraft < 0)) = NaN;
    ETA4((XC > -67) & (YC>-81) & (XC<-42) & (YC<-77.5) & (icedraft<=bathy)) = NaN;
    ETA4((XC<=-67) & (YC>-80.5) & (XC>-72) & (YC<-78.1) & (icedraft<=bathy)) = NaN;
    ETA4(ETA3>=eta_dc) = ETA3(ETA3>=eta_dc);
    y75 = find(YC(1,:)>-75,1);
    x75 = find(icedraft(:,y75)>bathy(:,y75),1);
    y78 = find(YC(1,:)>-78,1);
    x78 = find(icedraft(:,y78)>bathy(:,y78),1,'last');
    ymin = find(sum(icedraft-bathy,1)>0,1);
    xmin = find(sum(icedraft-bathy,2)>0,1);
    xmax = find(sum(icedraft(:,1:y78)-bathy(:,1:y78),2)>0,1,'last');
    ETA4(1:x75-1,y75:end) = eta_gl-0.5;
    ETA4(x78+1:end,y78:end) = eta_gl-0.5;
    ETA4(:,1:ymin-1) = eta_gl-0.5;
    ETA4(1:xmin-1,1:y75) = eta_gl-0.5;
    ETA4(xmax+1:end,1:y78) = eta_gl-0.5;
    ETA4 = inpaint_nans(ETA4);


    ETA = ETA1;    
    ETA(ETA1 <= eta_if) = ETA2(ETA1 <= eta_if);
    ETA(ETA<eta_mc) = ETA3(ETA<eta_mc);
    ETA(ETA<eta_dc) = ETA4(ETA<eta_dc);

  else
    

    ETA = ETA0;

  end
  
end

%%%
%%% Helper function to find minimum distance between a point and a contour
%%% returned by the contour function.
%%%
function dist = minDist (lat,lon,cntr)
  dist = 360;
  for m=1:size(cntr,2)
    dist_tmp = distance(lat,lon,cntr(2,m),cntr(1,m));
    if (dist_tmp < dist)
      dist = dist_tmp;
    end
  end
end


%%%
%%% Helper function to find the longest contour in a group of contours
%%% returned by the contour function.
%%%
function cntr = findLongestContour (C)
  idx = 1;
  maxlen = 0;
  maxidx = 2;
  while (idx < size(C,2))
    len = C(2,idx);
    if (maxlen<len)
      maxlen = len;
      maxidx = idx + 1;
    end
    idx = idx + len + 1;
  end
  cntr = C(:,maxidx:maxidx+maxlen-1);
end