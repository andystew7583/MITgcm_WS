%%%Barotropic Stream Function

%%% Load velocity
setExpname
loadexp

load SavedFiles/bsf_2cases.mat
load SavedFiles/TS_bot.mat
load SavedFiles/melt343445.mat
loadexp
%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
deltaT = 440;
% nIter0 = 0;
nDumps = round(nTimeSteps*(deltaT/dumpFreq));
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);
 
deltaT_4 = 300;
nIter0_4 = 1;
nDumps_4 = round(nTimeSteps*(deltaT_4/dumpFreq));
dumpIters_4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters_4 = dumpIters_4(dumpIters_4 >= nIter0_4);



deltaT_5 = 200;
nIter0_5 = 2877120;
nDumps_5 = round(nTimeSteps*10*(deltaT_5/dumpFreq));
dumpIters_5 = round((1:nDumps_5)*dumpFreq/deltaT_5); 
dumpIters_5 = dumpIters_5(dumpIters_5 >= nIter0_5);

deltaT2 = 440;
nDumps2 = round(nTimeSteps*10*(deltaT2/dumpFreq));
dumpIters_2 = round((1:nDumps2)*dumpFreq/deltaT2); 
dumpIters2 = dumpIters_2(dumpIters_2 >= nIter0);
nDumps2 = length(dumpIters2);

%%% Calculate time-averaged velocity
% tmin = 9*86400*360;
% tmax = 18*86400*360;
% 
% exppath1 = '/data3/MITgcm_WS/experiments/n_34452';
% exppath2 = '/data3/MITgcm_WS/experiments/a_34_20boundary';
% uuc = readIters(exppath1,'UVEL',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
% uuw = readIters(exppath2,'UVEL',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,Nr);
% 
% 
% %%% Grid spacing matrices
% DX = repmat(DXG,[1 1 Nr]);
% DY = repmat(DYG,[1 1 Nr]);
% DZ = repmat(reshape(DRF,[1 1 Nr]),[Nx Ny 1]);
% 
% %%% Calculate depth-averaged zonal velocity
% UUc= sum(uuc.*DZ.*hFacW,3);
% UUw = sum(uuw.*DZ.*hFacW,3);
% %%% Calculate barotropic streamfunction
% Psic = zeros(Nx+1,Ny+1);
% Psic(2:Nx+1,2:Ny+1) = -cumsum(UUc.*DYG,2);
% Psic = Psic(1:Nx,1:Ny);
% 
% Psiw = zeros(Nx+1,Ny+1);
% Psiw(2:Nx+1,2:Ny+1) = -cumsum(UUw.*DYG,2);
% Psiw = Psiw(1:Nx,1:Ny);
% 
% %%% if wanting to look at ice shelf in particular
rem = YC(1,:)>-75;
YG(:,rem)=[];
XG(:,rem)=[];


bathy(:,rem)=[];
SHELFICEtopo(:,rem)=[];
XC(:,rem)=[];
YC(:,rem)=[];  
Nx = size(XG,1);
Ny = size(YG,2);
Tc2(:,rem)=[];
Tw2(:,rem)=[];

Tc2(SHELFICEtopo==0)=NaN;
Tw2(SHELFICEtopo==0)=NaN;

FreshWaterw(:,rem)=[];
FreshWater(:,rem)=[];

% FreshWaterw(SHELFICEtopo==0)=NaN;
% FreshWater(SHELFICEtopo==0)=NaN;

topog_msk = ones(Nx,Ny);

for i=1:Nx
    for j=1:Ny       
           
           if (((SHELFICEtopo(i,j)) - (bathy(i,j)) < 0  ||  (bathy(i,j) == SHELFICEtopo(i,j)) ))
               bathy(i,j)=NaN;
               topog_msk(i,j) = 0;
           end
    end
end

% 
% 
% Psi_plotc = Psic;
% Psi_plotc(bathy>=SHELFICEtopo)=NaN;
% Psi_plotc(bathy==0) = NaN;
% Psi_plotc = Psi_plotc/1e6;
% Psi_plotc(Psi_plotc==0)=NaN;
% 
% Psi_plot = Psiw;
% Psi_plot(bathy>=SHELFICEtopo)=NaN;
% Psi_plot(bathy==0) = NaN;
% Psi_plot = Psi_plot/1e6;
% Psi_plot(Psi_plot==0)=NaN;
% 
% Psi_plot(SHELFICEtopo==0) = NaN;
% Psi_plotc(SHELFICEtopo==0)=NaN;
% for i = 1:Nx
%          for j = 1:Ny
%                 if (SHELFICEtopo(i,j)-bathy(i,j) <= 0)
%                     bathy(i,j) = NaN;
%                 end
%           end
% end


% topog_msk(Ny,:) = 1;


%%% Set up the figure
figure(2)
clf
set(gcf,'Position',[ 983    43   824   904]);

% fontsize = 18;
a=subplot(3,2,1);


axesm('eqaconicstd',...
  'fontsize',10,...
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




contourfm(YG,XG,Psi_plotc,11,'EdgeColor','none');

colormap(flipud(pmkmp(11)))

hold on
contourm(YG,XG,Psi_plotc,0:.5:5,'EdgeColor',[.5,.5,.5],'LineWidth',.6)
hold on


% [cs,C] = contourm(YC,XC,bathy,[-5000:1000:-1000 -500 -200 -100],'EdgeColor','black'); 
% hh = clabelm(cs,C);
% set (hh,'fontsize',10,'BackgroundColor','none','Edgecolor','none')
%          
% hold on
        
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
hold on   

% set(gca,'Position',[0.07 0.66 .475 .26]);


caxis([0 5])
ax2 = gca;
outerpos = ax2.OuterPosition;
ti = ax2.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax2_width = outerpos(3) - ti(1) - ti(3);
ax2_height = outerpos(4) - ti(2) - ti(4);
ax2.Position = [left-.005 bottom-.06 ax2_width+.05 ax2_height+.08];

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
% h = title('BSF (Sv), Initial 34.45psu FRIS Cavity','interpreter','latex','Fontsize',18);
% set(h,'Position',[0 -.95 0])
text(.115,-1.3,'\textbf{a}','fontsize',14,'interpreter','latex')

arrow([-0.025 -1.199],[-.0364 -1.213],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)
arrow([-0.01394 -1.224],[0 -1.24],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)
arrow([0.03815 -1.264],[0.0484 -1.247],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)
arrow([0.05334 -1.265],[0.04511 -1.286],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)

hold on

handle = title('REF','FontSize',16);
set(handle,'Position',[0.0000   -1.15    0.0000])




%%%%%%%%%%%%%%%%
b=subplot(3,2,2);



axesm('eqaconicstd',...
  'fontsize',10,...
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




contourfm(YG,XG,Psi_plot,11,'EdgeColor','none');

h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.9 0.7 0.01 .25])
title(h,'Sv','Fontsize',15,'interpreter','latex');
colormap(b,flipud(pmkmp(11)));  
colormap(h,flipud(pmkmp(11)));  


hold on
contourm(YG,XG,Psi_plot,0:.5:5,'EdgeColor',[.5,.5,.5],'LineWidth',.6)
hold on


% [cs,C] = contourm(YC,XC,bathy,[-5000:1000:-1000 -500 -200 -100],'EdgeColor','black'); 
% hh = clabelm(cs,C);
% set (hh,'fontsize',10,'BackgroundColor','none','Edgecolor','none')
%          
hold on
        
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
        
% hh = clabelm(cs,C);
% set (hh,'fontsize',10,'BackgroundColor','none')
% set(gca,'Position',[0.5 0.66 .475 .26]);

caxis([0 5])

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
% h2 = title('BSF (Sv), Initial 34psu FRIS Cavity','interpreter','latex','Fontsize',18);

ax3 = gca;
outerpos = ax3.OuterPosition;
ti = ax3.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax3_width = outerpos(3) - ti(1) - ti(3);
ax3_height = outerpos(4) - ti(2) - ti(4);
ax3.Position = [left-.07 bottom-.06 ax3_width+.05 ax3_height+.08];
text(.115,-1.3,'\textbf{b}','fontsize',14,'interpreter','latex')

handle = title('FRESH','FontSize',16);
set(handle,'Position',[0.0000   -1.15    0.0000])

arrow([-0.01942 -1.294],[-.03864 -1.281],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)
arrow([-0.05059 -1.218],[-0.04256 -1.206],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)
arrow([-0.0066 -1.278],[-0.01181 -1.26],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)
arrow([-0.01416 -1.24],[0.00699 -1.239],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)
arrow([0.03815 -1.264],[0.0484 -1.247],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)
arrow([0.05334 -1.265],[0.04511 -1.286],'Color',[.5 .5 .5],'Length',10,'Width',1,'TipAngle',25)



c=subplot(3,2,3)

axesm('eqaconicstd',...
  'fontsize',10,...
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


hold on
contourfm(YG,XG,Tc2,40,'EdgeColor','none');
% set(c,'Position',[0.07 0.35 .475 .26]);
colormap(c,jet(40));  
h = colorbar;
set(h,'Position',[0.9 0.375 0.01 .25])
set(gca,'FontSize',10);
title(h,'$^\circ$C','Fontsize',15,'interpreter','latex');
caxis([-2.5 -.5]);

hold on
    
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
% h = title('Mean Bottom Temperature ($^\circ$C), Initial 34.45psu Cavity','interpreter','latex','Fontsize',15);

% set(h,'Position',[0 -.95 0])

hold on
ax4 = gca;
outerpos = ax4.OuterPosition;
ti = ax4.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax4_width = outerpos(3) - ti(1) - ti(3);
ax4_height = outerpos(4) - ti(2) - ti(4);
ax4.Position = [left bottom-.07 ax4_width+.05 ax4_height+.08];
text(.115,-1.3,'\textbf{c}','fontsize',14,'interpreter','latex')

d=subplot(3,2,4)
axesm('eqaconicstd',...
  'fontsize',10,...
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




contourfm(YG,XG,Tw2,40,'EdgeColor','none');
colormap(d,jet(40));  
caxis([-2.5 -.5]);
         
hold on
        
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
% set(gca,'Position',[0.5 0.35 .475 .26]);
        

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
% h = title('Mean Bottom Temperature ($^\circ$C), Initial 34psu Cavity','interpreter','latex','Fontsize',15);

ax5 = gca;
outerpos = ax5.OuterPosition;
ti = ax5.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax5_width = outerpos(3) - ti(1) - ti(3);
ax5_height = outerpos(4) - ti(2) - ti(4);
ax5.Position = [left-.07 bottom-.07 ax5_width+.05 ax5_height+.08];
text(.115,-1.3,'\textbf{d}','fontsize',14,'interpreter','latex')

e=subplot(3,2,5)

axesm('eqaconicstd',...
  'fontsize',10,...
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


% 

contourfm(YG,XG,-FreshWater,100,'EdgeColor','none');
% set(gca,'Position',[0.07 0.02 .475 .26]);
colormap(e,redblue(100));  
caxis([-5 5])

% h = colorbar;
% set(gca,'FontSize',10);
% set(h,'Position',[0.8 0.1 0.02 .75])
% title(h,'$^\circ$C','Fontsize',20,'interpreter','latex');
         
hold on
        
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
        

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
% h = title('Mean Bottom Temperature ($^\circ$C), Initial 34psu Cavity','interpreter','latex','Fontsize',15);


ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom-.08 ax_width+.05 ax_height+.08];
text(.115,-1.3,'\textbf{e}','fontsize',14,'interpreter','latex')


f=subplot(3,2,6)
axesm('eqaconicstd',...
  'fontsize',10,...
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


hold on
contourfm(YG,XG,-FreshWaterw,100,'EdgeColor','none');
% set(gca,'Position',[0.5 0.02 .475 .26]);

colormap(f,redblue(100));  
h = colorbar;
set(h,'Position',[0.9 0.05 0.01 .25])
title(h,'m/yr','Fontsize',15,'interpreter','latex');
caxis([-5 5]);
set(gca,'FontSize',10);

hold on
    
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
% h = title('Shelf Ice Melt (m/yr),','interpreter','latex','Fontsize',15);

% set(h,'Position',[0 -.95 0])

hold on
ax6 = gca;
outerpos = ax6.OuterPosition;
ti = ax6.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax6_width = outerpos(3) - ti(1) - ti(3);
ax6_height = outerpos(4) - ti(2) - ti(4);
ax6.Position = [left-.07 bottom-.08 ax6_width+.05 ax6_height+.08];

text(.115,-1.3,'\textbf{f}','fontsize',14,'interpreter','latex')