%%%
%%% twistfunc.m
%%%
%%% Provides the degree of twisting (in degrees longitude) as a function of
%%% latitude (yc).
%%%
function twist = twistfunc (yc)

  %%% Load grid parameters
  defineGrid;

  %%% Twisting parameters
  ytmin = -79;
  ytmax = ymax;
%   twistmax = 15;
%   twistmin = 0;  
  twistmax = 0;
  twistmin = 0;
  
  %%% Construct twist function    
  twistrange = yc>ytmin & yc<ytmax;
  twist = zeros(1,Ny);
    
  %%% Linear twisting function
%   twist(twistrange) = twistmin + twistmax*(yc(twistrange)-ytmax)/(ytmin-ytmax);

  %%% Twisting function that creates a uniform distortion of the grid in
  %%% Cartesian space
  C1 = (twistmin-twistmax) / (log(tand(ytmax)+secd(ytmax)) - log(tand(ytmin)+secd(ytmin)));
  C2 = twistmin - C1*log(tand(ytmax)+secd(ytmax));
  twist(twistrange) = C1*log(tand(yc(twistrange))+secd(yc(twistrange))) + C2;
  
  
  twist(yc>=ytmax) = 0;
  twist(yc<=ytmin) = twistmax;
  
end

