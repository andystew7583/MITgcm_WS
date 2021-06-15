function contour_casts(castfiles,contours)
%%%
%%% USAGE: contour_cast_velocities (castfiles,contours)
%%%
%%% Creates a contour plot of the geostrophic velocities between neighboring
%%% casts. The plot will be created with castfiles{1} on the left and castfiles{end} on
%%% the right. Positive velocity indicates velocity into the page, and
%%% negative velocity indicates flow out of the page.
%%%
%%% Arguments:
%%% castfiles, - Cell array of .mat file names containing castdata objects. It is
%%%              assumed that the casts have fields called 'temperatures' 
%%%              and 'depths', which are equal-length arrays. The 'depths'
%%%              arrays must use the same depth intervals, though one may
%%%              be longer than the other (covering a deeper range of
%%%              depths). The cast objects must also have 'latitude' and
%%%              'longitude' fields with their locations.
%%% contours - Simply passed to the 'contours' argument of contourf.
%%%
%%%%%%%%% Now we have to put cast data together to plot and then plot
%%%%%%%%% against SOSE data


%%% Error-checking
if (length(castfiles) <= 2)
  error('castfiles must be a cell array containing at least three file names');
end

%%% Extract necessary fields from cast objects
N = length(castfiles)-1;
lat = zeros(1,N+1);
lon = zeros(1,N+1);
tt = zeros(1,N+1);
len = zeros(1,N);
zz = cell(1,N);
vv = cell(1,N);
for n=1:N+1
  load(castfiles{n});
  lat(n) = castdata.latitude;
  lon(n) = castdata.longitude;
  tt(n) = castdata.temperatures;
  if (n <= N)
    [vv{n},zz{n},unused] = calc_thermal_wind(castfiles{n},castfiles{n+1});
    len(n) = length(zz{n});
  end
end
latoff = lat(1);
lonoff = lon(1);
lat = 0.5 * (lat(1:end-1) + lat(2:end));
lon = 0.5 * (lon(1:end-1) + lon(2:end));


%%% Calculate distance between sections
xx = zeros(size(vv));
Rp = 6371000;
dlat = lat(1)-latoff;
dlon = lon(1)-lonoff;
mlat = 0.5*(lat(1)+latoff);
dx = dlon*(2*pi/360)*Rp*cos(mlat*2*pi/360);
dy = dlat*(2*pi/360)*Rp;
dl = sqrt(dx^2+dy^2);
xx(1) = dl;
for n=2:length(vv)
  dlat = lat(n)-lat(n-1);
  dlon = lon(n)-lon(n-1);
  mlat = 0.5*(lat(n)+lat(n-1));
  dx = dlon*(2*pi/360)*Rp*cos(mlat*2*pi/360);
  dy = dlat*(2*pi/360)*Rp;
  dl = sqrt(dx^2+dy^2);
  xx(n) = xx(n-1) + dl;
end

%%% Extend casts to length of deepest cast and combine into a matrix
[maxlen,idx] = max(len);
% depths = zz{idx(1)};
% [ZZ XX] = meshgrid(depths,xx);
XX = repmat(xx',[1 maxlen]);
ZZ = zeros(N,maxlen);
VV = zeros(N,maxlen);
TT = zeros(N,maxlen);
for n=1:N
  if (len(n) < maxlen)    
%     vv_old = vv{n};
%     vv_new = vv{idx(1)};
%     vv_new(1:len(n)) = vv_old;
%     vv_new(len(n)+1:end) = NaN;  
%     vv{n} = vv_new;

    zz_old = zz{n};
    vv_old = vv{n};
    tt_old = tt{n};
    
    dz_new = (zz_old(end)-zz_old(1))/(maxlen-1);
    zz_new = zz_old(1):dz_new:zz_old(end);
    
    tt_new = interp1(zz_old,tt_old,tt_new,'linear','extrap');
    vv_new = interp1(zz_old,vv_old,zz_new,'linear','extrap');
    
    tt{n} = tt_new;
    vv{n} = vv_new;
    zz{n} = zz_new;
  end
%   ZZ(n,:) = depths;
  ZZ(n,:) = zz{n};
  VV(n,:) = vv{n};
  TT(n,:) = tt{n};
end

%%% Create the plot
figure(1);
contourf(XX/1000,ZZ,VV,contours);
set(gca,'YDir','Reverse');
colormap cool;
colorbar;
caxis([-0.3 0.3]);
set(gca,'FontSize',16);
xlabel('Distance (km)');
ylabel('Depth (m)');

figure(2);
contourf(XX/100,ZZ,TT,contours);
set(gca,'YDir','Reverse');
colormap cool
colorbar
set(gca,'Fontsize',16);
xlabel('Distance km');
ylabel('Depth (m)');


end

