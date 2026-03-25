
%%%%%making a histogram of winds

run ../newexp/defineGrid.m

days =3287;

inputpath = fullfile(gendir,'/MITgcm_WS/newexp/differentResolutions/a_128_orig');
inputpath2 = fullfile(gendir,'/MITgcm_WS/newexp/differentResolutions/a_128_orig');

zonal = zeros(Nx,Ny,days);
merid = zeros(Nx,Ny,days);
fid = fopen(fullfile(inputpath,zwind),'r','b');
fid2 = fopen(fullfile(inputpath,mwind),'r','b');
for k=1:days
  zonal(:,:,k) = fread(fid,[Nx Ny],'real*8');
  merid(:,:,k) = fread(fid2,[Nx Ny],'real*8');
end
fclose(fid);
fclose(fid2);

fid = fopen(fullfile(inputpath2,'bathyFile.bin'),'r','b');
bathy = fread(fid,[Nx Ny],'real*8');
fclose(fid);

fid = fopen(fullfile(inputpath2,SHELFICEtopoFile),'r','b');
SHELFICEtopo = fread(fid,[Nx Ny],'real*8');
fclose(fid);

figure(1)
zonal(zonal==0)=NaN;
z = histogram(zonal,'Normalization','pdf');
title('Pdf of Zonal Wind, Control','Interpreter','Latex');


% zonal_timemean = mean(zonal,1);
zonal_timemean = squeeze(nanmean(zonal,3));
zonal_timemean(SHELFICEtopo-bathy==0)=NaN;
zonal_timemean(zonal_timemean==0)=NaN;
figure(2)
pcolor(XMC,YMC,zonal_timemean'),shading interp;
title('Time-mean Zonal Wind, Control','Interpreter','Latex');
colormap jet(30)
caxis([-15 5]);
colorbar

figure(3)
merid(merid==0)=NaN;
m = histogram(merid,'Normalization','pdf');
title('Pdf of Meridional Wind, Control','Interpreter','Latex');


merid_timemean = squeeze(nanmean(merid,3));
merid_timemean(SHELFICEtopo-bathy==0)=NaN;


figure(4)
pcolor(XMC,YMC,merid_timemean'),shading interp;
title('Time-mean Meridional Wind, Control','Interpreter','Latex');
colormap jet(30)
colorbar
caxis([-15 10]);