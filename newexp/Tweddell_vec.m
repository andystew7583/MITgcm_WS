%%%v
%%% Tweddell_vec.m
%%%
%%% "Twists" lat/lon gridded data longitudinally to create a modified
%%% Weddell Sea domain. This version of the function applies to vector
%%% data, and additionally performs rotation of the vector as appropriate
%%% for the distortion created by the "twist".
%%%
%%% N.B. Assumes longitude ranges from -180 to 180, and that u and v
%%% gridpoints are collocated (u and v zonal and meridional components of
%%% the vector).
%%%
function [u_twist,v_twist] = Tweddell_vec (xc,yc,u,v)

  %%% Input grid sizes
  Ny = length(yc);
  Nz = size(u,3);
  
  %%% First do the twist on both u and v
  u_twist = Tweddell(xc,yc,u);
  v_twist = Tweddell(xc,yc,v);
  
  %%% Adjust the twisted vector according to the distortion caused by the
  %%% twisting
  for j=1:Ny      
    
    %%% Calculate gradient of twisting function - this defines the
    %%% distortion angle, phi
    dy = 0.01;
    dtwist_dy = (twistfunc(yc(j)+dy/2)-twistfunc(yc(j)-dy/2))/dy;
    phi = atan(dtwist_dy*cosd(yc(j))); %%% Distortion angle
    
    %%% If there is no distortion here then don't bother
    if (phi==0)
      continue;
    end
    
    %%% Rotate the vector and renormalize
    for k=1:Nz
      
      u_vec = u_twist(:,j,k);
      v_vec = v_twist(:,j,k);
      u_abs = sqrt(u_vec.^2+v_vec.^2+2*u_vec.*v_vec*sin(phi));
      if (u_abs > 0)
        u_twist(:,j,k) = (u_vec + v_vec*sin(phi)) / u_abs;
        v_twist(:,j,k) = v_vec*cos(phi) / u_abs;
      end
      
    end
  end

end