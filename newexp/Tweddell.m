%%%v
%%% Tweddell.m
%%%
%%% "Twists" lat/lon gridded data longitudinally to create a modified
%%% Weddell Sea domain.
%%%
%%% N.B. Assumes longitude ranges from -180 to 180.
%%%
function twisted_data = Tweddell (xc,yc,mata)

  %%% Load grid configuration
  defineGrid;
  
  %%% Load the twist function
  Ny = length(yc);
  Nz = size(mata,3);
  twist = twistfunc(yc);
  
  %%% Do the twist
  twisted_data = mata;
  for j=1:Ny      
    if (twist(j)==0)
      continue;
    end
    xc_twist = xc+twist(j);
    outrange = xc_twist>180; 
    inrange = xc_twist<=180; 
    xc_twist(outrange) = xc_twist(outrange) - 360;
    xc_twist = [xc_twist(outrange)' xc_twist(inrange)']';    
    for k=1:Nz
      data_vec = [mata(outrange,j,k)' mata(inrange,j,k)']';
      twisted_data(:,j,k) = interp1(xc_twist,data_vec,xc,'linear');
    end
  end

end