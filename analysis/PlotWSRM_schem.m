
%%%
%%%  make schematic of WSRM.
%%%load data
setExpname
%%% Set true if plotting on a Mac
mac_plots=0;

%%%load shelf ice elevation data
load ../newexp/ELEV.mat
%%% Read experiment data
loadexp;

%%% Select diagnostic variable to animate
diagnum = 14;
outfname =diag_fileNames{1,diagnum};

dumpFreq = abs(diag_frequency(diagnum));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((0:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
% Mice_elev=Mice_elev';
for i = 1:Nx
    for j = 1:Ny
        if (bathy(i,j)-SHELFICEtopo(i,j)>=0)
            bathy(i,j)=NaN;
            SHELFICEtopo(i,j)=NaN;
        end  
    end
end
topog_msk = ones(Nx,Ny);

for i=1:Nx
    for j=1:Ny       
           
           if (((SHELFICEtopo(i,j)) - (bathy(i,j)) < 0  ||  (bathy(i,j) == SHELFICEtopo(i,j)) ))
               topog_msk(i,j) = 0;
           end
    end
end
bathy(bathy==0)=NaN;
% SHELFICEtopo(SHELFICEtopo==0) = NaN;
% Mice_elev(Mice_elev==0) = NaN;



%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 18;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.15 0.15 0.7 0.76];
  framepos = [100    500   800  800];
end

handle = figure(20);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);





%%%%tried to get 3D image working....
% set(gcf,'color','w');
% clf;
% p = surf(XC,YC,bathy); 
% p.FaceColor = 'blue';
% hold on
% %         colorbar
% p1 = surf(XC,YC,SHELFICEtopo);shading flat
% hold on
% p1.FaceColor = [.5 .5 .5];
% hold on
% p2 = surf(XC,YC,Mice_elev),shading interp;
% hold on
% p2.FaceColor = [48 129 238]/256;
%         
%        
% hold off;        
% axis([-80 20 -84 -64])
% view(170,40);
% lighting gouraud;
% camlight('headlight');        
% xlabel('x (Longitude)','interpreter','latex','fontsize',15);
% ylabel('y (Latitude)','interpreter','latex','fontsize',15);
% zlabel('z (m)','interpreter','latex','fontsize',15);
% set(gca,'FontSize',16);
% %         set(gcf,'Color','w');
        
a=subplot(2,2,[1 2])
axesm('eqaconicstd',...
  'fontsize',8,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[-84 -64], ...
  'MapLonLimit',[-80 20], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 2, ...
  'MLineLocation', 5,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)

contourfm(YC,XC,bathy,20);
shading interp
colormap(a,haxby(20)); 
handle(1)=colorbar;
set(handle,'Position',[0.85 0.55 0.01 .4]);
hold on
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1);
text(.18,-1.4,'\textbf{a}','interpreter','latex','fontsize',18)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left-.05 bottom-.05 ax_width+.15 ax_height+.05];
h=title('WSRM Bathymetry (m)','fontsize',14,'interpreter','latex');
set(h,'Position',[0 -.95 0]);

b=subplot(2,2,3)

% rem = YC(1,:)>-75;
% 
% YG(:,rem)=[];
% XG(:,rem)=[];
% XC(:,rem)=[];
% YC(:,rem)=[];  
 
for i=1:Nx
    for j = 1:Ny
        if YC(i,j)>-75
            SHELFICEtopo(i,j)=0;
            bathy(i,j)=0;
        end
    end
end



Nx = size(XG,1);
Ny = size(YG,2);


axesm('eqaconicstd',...
  'fontsize',10,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[-84 -75], ...
  'MapLonLimit',[-80 -25], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 2, ...
  'MLineLocation', 5,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)


contourfm(YC,XC,SHELFICEtopo,20)
shading interp 
colormap(b,flipud(brewermap(12,'BuPu')));
caxis([-1400 0])
mandle= colorbar;
set(mandle,'Position',[.48,.13,.01,.3]);
text(.06,-1.35,'\textbf{b}','interpreter','latex','fontsize',18)

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left-.025 bottom ax_width+.05 ax_height];
m=title('FRIS Cavity Ice Draft Depth (m)','fontsize',14,'interpreter','latex');
set(m,'Position',[0,-1.158,0]);






c=subplot(2,2,4);

axesm('eqaconicstd',...
  'fontsize',10,...
  'Grid','on', ...    
  'Frame','off', ...
  'MapLatLimit',[-84 -75], ...
  'MapLonLimit',[-80 -25], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 2, ...
  'MLineLocation', 5,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)


contourfm(YC,XC,(SHELFICEtopo-bathy),15)
shading interp 
colormap(c,brewermap(20,'Blues'));
mandle= colorbar;
set(mandle,'Position',[.94,.13,.01,.3]);

text(.06,-1.35,'\textbf{c}','interpreter','latex','fontsize',18)

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width+.025 ax_height];
s=title('FRIS Cavity Water Column Thickness (m)','fontsize',14,'interpreter','latex');
set(s,'Position',[0,-1.158,0]);