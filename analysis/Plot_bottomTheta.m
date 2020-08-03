%%%%%%%

%%%plotting bottom theta values for FRESH and REF simulations
load SavedFiles/thetavals_343445



        
figure(2)
clf
scrsz = get(0,'ScreenSize');
% fontsize = 18;
hold on
p1=subplot(2,1,1);

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
setm(gca,'MLabelParallel',-20)




contourfm(YG,XG,Tc2,40,'EdgeColor','none');
set(gca,'Position',[0.2 0.32 .6 .8]);
colormap jet(40);  
h = colorbar;
set(h,'Position',[0.85 0.1 0.02 .75])
set(gca,'FontSize',10);
title(h,'$^\circ$C','Fontsize',20,'interpreter','latex');

hold on

         
        
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
% h = title('Mean Bottom Temperature ($^\circ$C), 34.45psu FRIS Cavity Restoring','interpreter','latex','Fontsize',15);
h = title('Mean Bottom Temperature ($^\circ$C), Initial 34.45psu Cavity','interpreter','latex','Fontsize',15);

set(h,'Position',[0 -.95 0])

hold on
% ax2 = gca;
% outerpos = ax2.OuterPosition;
% ti = ax2.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax2_width = outerpos(3) - ti(1) - ti(3);
% ax2_height = outerpos(4) - ti(2) - ti(4);
% ax2.Position = [left+.03 bottom ax2_width-.03 ax2_height];

text(.3,-1.3,'a','fontsize',14,'interpreter','latex')

hold on

%%%%%%%%%%%%%%%%
p2=subplot(2,1,2);


% 
% latMin = min(min(YC));
% latMax = -60;
% lonMin = min(min(XC));
% lonMax = max(max(XC));

% 'FLatLimit', [latMin latMax], ...
%           'FLonLimit', [lonMin lonMax], ...


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
setm(gca,'MLabelParallel',-20)




contourfm(YG,XG,Tw2,40,'EdgeColor','none');
set(gca,'Position',[0.2 0.0 .6 .45]);
colormap jet(100);  

% h = colorbar;
% set(gca,'FontSize',10);
% set(h,'Position',[0.8 0.1 0.02 .75])
% title(h,'$^\circ$C','Fontsize',20,'interpreter','latex');
         
hold on
        
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
        

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
% h = title('Mean Bottom Temperature ($^\circ$C), 34psu FRIS Cavity Restoring','interpreter','latex','Fontsize',15);
h = title('Mean Bottom Temperature ($^\circ$C), Initial 34psu Cavity','interpreter','latex','Fontsize',15);

set(h,'Position',[0 -.95 0])

% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% text(.3,-1.3,'b','fontsize',14,'interpreter','latex')

