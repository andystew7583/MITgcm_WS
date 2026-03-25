% %%%% plotting Validation Locatiovns
% 
% [im1,map1] = imread('KN_pic.jpg');
% [im2,map2] = imread('woce_pick.jpg');
% 
% figure(1)
% clf
% w = subplot(2,1,1);
% imshow(im1,map1);
% set(w,'Position',[.2 .45 .6 .55]);
% title('KN section, from A. Behrendt, W. Dierking, and H. Witte (2015)','interpreter','latex','fontsize',14)
% hold on
% w1=subplot(2,1,2)
% imshow(im2,map2)
% set(w1,'Position',[.2 .0 .6 .45]);
% title('WOCE A12 and SR4 Sections, from R. Kerr et al. (2018)','interpreter','latex','fontsize',14)





loadexp
%%% Set up the figure
figure(2)
clf
scrsz = get(0,'ScreenSize');
fontsize = 18;





latMin = min(min(YC));
latMax = -60;
lonMin = min(min(XC));
lonMax = max(max(XC));

%'FLatLimit', [latMin latMax], ...
     %frrrrrrrrrrrrrrrrrrrr     'FLonLimit', [lonMin lonMax], ...


axesm('eqaconicstd',...
  'fontsize',13,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[-84 -64], ...
  'MapLonLimit',[-80 20], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
% setm(gca,'MLabelParallel',-20)

bathy(SHELFICEtopo-bathy==0)=NaN;
% bathy(SHELFICEtopo<0)=NaN;
% set(gca,'color',[.5 .5 .5]);



hold on
h = NaN(Nx,Ny);
cc = NaN(Nx,Ny);
mm = NaN(Nx,Ny);

for i = 1:Nx
     for j = 1:Ny
         if  XC(i,j) <-30 && XC(i,j)>-80
             
             if YC(i,j) <-75
                 if SHELFICEtopo(i,j)<0
                   if bathy(i,j)<0  
                   cc(i,j)=XC(i,j);
                   mm(i,j)=YC(i,j);
                 
                   h(i,j) = bathy(i,j);
                   end
                 end
             end
         end
     end
end



pcolorm(YC,XC,bathy),shading interp
set(gca,'Position',[0.07 0.03 .85 .92]);
colormap haxby(100);  
hold on
% contourm(YC,XC,(h),50,'linestyle',':','linewidth',1,'linecolor','k');
l=fillm(mm,cc,'k','markerfacecolor','k');

colorbar

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left+.05 bottom ax_width-.08 ax_height];
% hp=findobj(h,'type','patch');
% hatchfill(hp);
% ax=pcolor(XC,YC,h),shading flat;
set(gca,'color',[.5 .5 .5]);
title('WSRM FRIS Analysis Region, Overlaid ontop of Bathymetry','fontsize',15,'interpreter','latex')
xlabel('longitude','fontsize',15,'interpreter','latex');
ylabel('latitude','fontsize',15,'interpreter','latex');
