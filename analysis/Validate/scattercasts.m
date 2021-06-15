%%%
%%% scattercasts.m
%%%
%%% Function to take CTD casts and interpolates them
%%% onto our grid
%%%


function [casts] = scattercasts(datafile,lon,lat,W_LO,W_LA)

  % Interpolate U Grd 
  xsize = size(lat,1);
  ysize = size(lat,2);
                   
  lat = reshape(lat,xsize*ysize,1);
  lon = reshape(lon,xsize*ysize,1);
  V = datafile;
  V = reshape(V,xsize*ysize,1); 
  

  F = scatteredInterpolant(double(lon),double(lat),double(V),'natural');

                 
  %%% Store the velocity data in the storage matrix for Ugrd
  casts(:,:) = F(double(W_LO),double(W_LA));
end
      