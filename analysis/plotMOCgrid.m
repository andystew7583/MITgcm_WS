%%% 
%%% plotMOCgrid.m
%%%
%%% Plots the coordinate transformation used for our MOC calculations.
%%%

%%% Load experiment
% loadexp;

ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,false,true);

bathy_plot = bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
ETA_plot = ETA;
ETA_plot(SHELFICEtopo-bathy<=0) = NaN;


figure(99);
clf;
set(gcf,'Position',[243         485        1077         438]);
pcolor(XC,YC,SHELFICEtopo-bathy_plot);
shading interp;
colormap(flip(haxby,1))
colorbar;
hold on
% [C,h]=contour(XC,YC,ETA_plot,[-86:0.2:-70],'EdgeColor','k');
[C,h]=contour(XC,YC,ETA_plot,[-9:.1:-.1 -0.09:0.01:0.09 0.1:0.1:20],'EdgeColor','k');
% [C,h]=contour(XC,YC,ETA_plot,[0:.5:2 3:1:20],'EdgeColor','k');
hold off;
clabel(C,h);
caxis([0 5000]);
ylabel('Latitude');
xlabel('Longitude');
title('MOC coordinate \eta');
set(gca,'Position',[0.0743    0.1096    0.8292    0.8196]);

