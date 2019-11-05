
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
Mice_elev=Mice_elev';
for i = 1:Nx
    for j = 1:Ny
        if (bathy(i,j)-SHELFICEtopo(i,j)>=0)
            bathy(i,j)=NaN;
            SHELFICEtopo(i,j)=NaN;
%             Mice_elev(i,j)=NaN;
        end  
    end
end

% bathy(bathy==0)=NaN;
SHELFICEtopo(SHELFICEtopo==0) = NaN;
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
set(gcf,'color','w');
clf;
p = surface(XC,YC,bathy); 
p.EdgeColor = 'None';
p.FaceColor = 'TextureMap';
hold on
%         colorbar
p1 = surface(XC,YC,SHELFICEtopo);
% shading flat
hold on
p1.EdgeColor = 'None';
p1.FaceColor = [.5 .5 .5];
hold on
% p2 = surface(XC,YC,Mice_elev);
% % shading interp;
% hold on
% p2.EdgeColor = 'None';
% p2.FaceColor = [48 129 238]/256;
        
       
hold off;        
axis([-80 20 -84 -64])
view(170,40);
lighting gouraud;
camlight('headlight');        
xlabel('x (Longitude)','interpreter','latex','fontsize',15);
ylabel('y (Latitude)','interpreter','latex','fontsize',15);
zlabel('z (m)','interpreter','latex','fontsize',15);
set(gca,'FontSize',16);
%         set(gcf,'Color','w');
        

        
  
  
