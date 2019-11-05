%%%
%%% Surfaceforcing.m
%%%
%%% Function to take polar wrf data and create 12 month climatology, interpolates for
%%% each data variable
%%%
%%% Input will be PolarWrf data, size of grize lat/lon of designated grid
%%% (highest res of polar wrf ) xxxxxx

function [sfc_forcing,F] = force(datafile,lon,lat,W_LO,W_LA,F)

  % Interpolate U Grd 
  xsize = size(lat,1);
  ysize = size(lat,2);
                   
  lat = reshape(lat,xsize*ysize,1);
  lon = reshape(lon,xsize*ysize,1);
  V = datafile;
  V = reshape(V,xsize*ysize,1); 
  
  if (isempty(F))
    F = scatteredInterpolant(double(lon),double(lat),double(V),'natural');
  else
    F.Values = double(V); 
  end
                 
  %%% Store the velocity data in the storage matrix for Ugrd
  sfc_forcing(:,:) = F(double(W_LO),double(W_LA));
end
                
                

  