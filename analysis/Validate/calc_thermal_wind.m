function [velocity transport] = calc_thermal_wind (dens,lats,depths)
%%%
%%% USAGE: [velocity, transport] = calc_thermal_wind (dens, lats, depths)
%%%
%%% Calculates the thermal wind velocity from the potential density
%%%
%%% Arguments:
%%% dens - matrix of in situ density values 
%%% lats - vector of latitudes
%%% depths - vector of depths
%%%    
%%% dens must be an M x N matrix. lats must be a vector of length M. 
%%% depths must be a vector of length N.
%%%

%%% Error checking
if (length(lats) ~= size(dens,1))
  error('make_TS_plot: Length of latitude data must match number of rows of density data');
end
if (length(depths) ~= size(dens,2))
  error('make_TS_plot: Length of depth data must match number of columns of density data');
end

%%% Easier variable names
zz = depths;
dd = dens;
Ny = size(dd,1);
Nz = size(dd,2);

%%% Fixed parameters
Rp = 6371000;
rho0 = 1027;
g = 9.81;
Omega = 2*pi*366/365/86400;
f = 2*Omega*sin(2*pi*lats/360);

%%% Convert latitudes into distances (approximately)
yy = 0*lats;
dy = (lats(2:Ny)-lats(1:Ny-1))*(2*pi/360)*Rp;
yy(2:Ny) = cumsum(dy);

%%% Crudely approximate thermal wind
uu = NaN*dd;
T = 0;
for j=2:Ny-1;
  for k=Nz:-1:1
    
    %%% If we're in the bottom topography, ignore this point
    if (isnan(dd(j,k)+dd(j-1,k)+dd(j+1,k)))
      continue;
    end
    
    %%% Set velocity to zero just above bottom topography
    if ((k==Nz) || isnan(uu(j,k+1)))
      uu(j,k) = 0;
      continue;
    end
    
    %%% Otherwise calculate thermal wind shear
    dz = zz(k)-zz(k+1);
    uu(j,k) = uu(j,k+1) + dz*(g/rho0/f(j))*(dd(j+1,k)-dd(j-1,k))/(dy(j)+dy(j-1));
    T = T + uu(j,k)*dz*(dy(j)+dy(j-1))/2;
    
  end
end

%%% Return the calculated velocity
velocity = uu;
transport = T;