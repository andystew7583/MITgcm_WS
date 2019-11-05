%%%v
%%% Calculate and plot time-mean BSF
%%%



%%% Load velocity

loadexp;
deltaT=440;
% nIter0 = 0;
%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Calculate time-averaged velocity
tmin = 17*86400*360;
tmax = 26*86400*360;
uu = readIters(exppath,'UVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);



%%% Grid spacing matrices
DX = repmat(DXG,[1 1 Nr]);
DY = repmat(DYG,[1 1 Nr]);
DZ = repmat(reshape(DRF,[1 1 Nr]),[Nx Ny 1]);

%%% Calculate depth-averaged zonal velocity
UU = sum(uu.*DZ.*hFacW,3);

%%% Calculate barotropic streamfunction
Psi = zeros(Nx+1,Ny+1);
Psi(2:Nx+1,2:Ny+1) = -cumsum(UU.*DYG,2);
Psi = Psi(1:Nx,1:Ny);



%%% if wanting to look at ice shelf cavity in particular -->
rem = YC(1,:)>-75;





YG(:,rem)=[];
XG(:,rem)=[];
Psi(:,rem)=[];
bathy(:,rem)=[];
SHELFICEtopo(:,rem)=[];
XC(:,rem)=[];
YC(:,rem)=[];  
 



Nx = size(XG,1);
Ny = size(YG,2);

topog_msk = ones(Nx,Ny);

for i=1:Nx
    for j=1:Ny       
           
           if (((SHELFICEtopo(i,j)) - (bathy(i,j)) < 0  ||  (bathy(i,j) == SHELFICEtopo(i,j)) ))
               topog_msk(i,j) = 0;
           end
    end
end


Psi(SHELFICEtopo==0)=NaN;

Psi_plot = Psi;
Psi_plot(bathy>=SHELFICEtopo)=NaN;
Psi_plot(bathy==0) = NaN;
Psi_plot = Psi_plot/1e6;
Psi_plot(Psi_plot==0)=NaN;


for i = 1:Nx
         for j = 1:Ny
                if (SHELFICEtopo(i,j)-bathy(i,j) <= 0)
                    bathy(i,j) = NaN;
                end
          end
end


% topog_msk(Ny,:) = 1;


%%% Set up the figure
figure(2)
clf
scrsz = get(0,'ScreenSize');
% fontsize = 18;




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
  'MapLatLimit',[-84 -75], ...
  'MapLonLimit',[-80 -30], ... 
  'MapParallels',[-85 -65], ...
  'PLineLocation', 5, ...
  'MLineLocation', 10,...
  'MeridianLabel', 'on', ...
  'ParallelLabel', 'on');    %, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
axis off;
setm(gca,'MLabelParallel',-20)




contourfm(YG,XG,Psi_plot,40,'EdgeColor','none');
set(gca,'Position',[0.07 0.03 .85 .92]);
colormap haxby(100);  
colormap(flipud(colormap));       

h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.92 0.1 0.02 .8])
title(h,'Sv','Fontsize',20,'interpreter','latex');

hold on
contourm(YG,XG,Psi_plot,1:.2:5,'EdgeColor',[.5,.5,.5],'LineWidth',.6)
hold on


[cs,C] = contourm(YC,XC,bathy,[-5000:1000:-1000 -500 -200 -100],'EdgeColor','black'); 
hh = clabelm(cs,C);
set (hh,'fontsize',10,'BackgroundColor','none','Edgecolor','none')
         
hold on
        
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
        
hh = clabelm(cs,C);
set (hh,'fontsize',10,'BackgroundColor','none')

%     hh = clabelm(cs,C);
%     set (hh,'fontsize',7,'BackgroundColor','none')

caxis([0 5])

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
h = title('Barotropic Stream Function (Sv) (FRIS Cavity)','interpreter','latex','Fontsize',18);
set(h,'Position',[0 .9 0])

hold on



%%%%%%%%%%%%%%%%

% contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)




