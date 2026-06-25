fname = '/Volumes/Stewart-RAID1-A/UCLA/Data/OCEAN_ICE_climatology/OI_Climatology.nc';
info = ncinfo(fname);
longitude = ncread(fname,'longitude');
latitude = ncread(fname,'latitude');
pressure = ncread(fname,'pressure');
ct = ncread(fname,'ct');
sa = ncread(fname,'sa');

clim = [-2.6 -1];
cmap = cmocean('balance',32);
cmap = cmap(8:23,:);
fontsize = 18;

figure(300);
axesm('stereo',...
  'fontsize',18,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[-90 -45], ...
  'MapLonLimit',[-180 180], ...   
  'PLineLocation', 10, ...
  'MLineLocation', 60,...
  'MeridianLabel', 'off', ...
  'ParallelLabel', 'off');  
axis off;
[XC,YC] = meshgrid(longitude,latitude);
pcolorm(YC,XC,ct(:,:,25));
shading interp
colormap(gca,cmocean('thermal',70));
set(gca,'FontSize',18);
caxis([-2 2]);
h = colorbar;
title(h,'$^\circ$C','Fontsize',18,'interpreter','latex');


%%% Define grid to extract data sections
startLat = -82.3;
endLat = -72;
startLon = -62;
endLon = -40;
Nsec = 2001;
dLat = (endLat-startLat)/(Nsec-1);
dLon = (endLon-startLon)/(Nsec-1);
secLats = startLat:dLat:endLat;
secLons = startLon:dLon:endLon;
Nr = length(pressure);

%%% Extract data along defined sections
secS = zeros(Nsec,Nr);
secT = zeros(Nsec,Nr);
for n=1:Nsec
  jm = find(secLats(n)>=latitude,1,'last');
  im = find(secLons(n)>=longitude,1,'last');
  jp = jm+1;
  ip = im+1;
  wp_y = (secLats(n)-latitude(jm))/(latitude(jp)-latitude(jm));
  wm_y = 1-wp_y;
  wp_x = (secLons(n)-longitude(im))/(longitude(ip)-longitude(im));
  wm_x = 1-wp_x;
  if (wp_x > wm_x)
    xidx = ip;
  else
    xidx = im;
  end
  if (wp_y > wm_y)
    yidx = jp;
  else
    yidx = jm;
  end
  secS(n,:) = squeeze(sa(yidx,xidx,:));
  secT(n,:) = squeeze(ct(yidx,xidx,:));     
end


[ZZ_sec,LA_sec] = meshgrid(pressure,secLats);
figure(400);
pcolor(LA_sec,ZZ_sec,secT);
set(gca,'YDir','reverse');
shading interp;
hold on
[C,h] = contour(LA_sec,ZZ_sec,secS,[34.1:.1:35],'EdgeColor','k');
clabel(C,h,'FontSize',fontsize-4,'LabelSpacing',300);
hold off;
set(gca,'YLim',[0 1400]);
set(gca,'XLim',[startLat (-73)]);
caxis(clim);
colormap(gca,cmap);
cbhandle = colorbar;
% set(cbhandle,'Position',cb_pos);
title(cbhandle,'$^\circ$C','Fontsize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);
xlabel('Latitude')
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);
% title(['Instantaneous potential temperature, salinity, ',datestr(datenum('01-Jan-2008')+iter*deltaT/86400)]);
