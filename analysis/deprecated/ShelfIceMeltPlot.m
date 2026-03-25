%%%%%  Plotting Melt Rate from Ice Shelves
%%%%%% Plotting fresh water flux/convert to ice melt rate GT/year



%%%%%%%has to be halfshelfice experiment

%%% Read experiment data
setExpname
loadexp;



%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(15));
deltaT = 440;
% nIter0 = 0;
nDumps = round(nTimeSteps*(deltaT/dumpFreq));
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);
 
deltaT_4 = 200;
nIter0_4 = 1;
nDumps_4 = round(nTimeSteps*(deltaT_4/dumpFreq));
dumpIters_4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters_4 = dumpIters_4(dumpIters_4 >= nIter0_4);



deltaT_5 = 200;
nIter0_5 = 2877120;
nDumps_5 = round(nTimeSteps*10*(deltaT_5/dumpFreq));
dumpIters_5 = round((1:nDumps_5)*dumpFreq/deltaT_5); 
dumpIters_5 = dumpIters_5(dumpIters_5 >= nIter0_5);

deltaT2 = 400;
nDumps2 = round(nTimeSteps*10*(deltaT2/dumpFreq));
dumpIters_2 = round((1:nDumps2)*dumpFreq/deltaT2); 
dumpIters2 = dumpIters_2(dumpIters_2 >= nIter0);
nDumps2 = length(dumpIters2);

mac_plots = false;

%%% Select diagnostic variable 
diagnum = 45;
outfname = diag_fileNames{1,diagnum};

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 1;


rho_I = 917;
rho_W = 1027;

%%% Longitudinal positions

XC;


%%% Latitudinal positions

YC;


% nDump_start = 15*12;
% nDump_Finish = 18*12;



topog_msk = ones(Nx,Ny);

for i=1:Nx
    for j=1:Ny        
           
           if (((SHELFICEtopo(i,j)) - (bathy(i,j)) < 0  ||  (bathy(i,j) == SHELFICEtopo(i,j)) ))
               topog_msk(i,j) = 0;
           end
    end
end

tmin = 18*86400*360;
tmax = 25*86400*360;

exppath1 = '/data3/MITgcm_WS/experiments/a_3445_20boundary';
exppath2 = '/data3/MITgcm_WS/experiments/a_34_20boundary';
freshwaterc = readIters(exppath1,'SHIfwFlx',dumpIters2,deltaT2,tmin,tmax,Nx,Ny,1);
freshwaterw = readIters(exppath2,'SHIfwFlx',dumpIters,deltaT,tmin,tmax,Nx,Ny,1);



% freshwater = NaN(Nx,Ny,nDump_Finish-nDump_start);
% for n=1:(nDump_Finish-nDump_start)
%   
%   SHImelt = rdmdsWrapper(fullfile(exppath,'/results/SHIfwFlx'),dumpIters(n+nDump_start-1));
%   freshwater(:,:,n) = (SHImelt);
%   
% end
  
% freshwater_season = zeros(size(freshwater,1),size(freshwater,2),3,4);
% 
% for n = 1:3
%   freshwater_season(:,:,n,1) = nanmean(freshwater(:,:,(n-1)*12+[12 1 2]),3);
%   freshwater_season(:,:,n,2) = nanmean(freshwater(:,:,(n-1)*12+(3:5)),3);
%   freshwater_season(:,:,n,3) = nanmean(freshwater(:,:,(n-1)*12+(6:8)),3);
%   freshwater_season(:,:,n,4) = nanmean(freshwater(:,:,(n-1)*12+(9:11)),3);
%   
% end
% freshwater_season = squeeze(mean(freshwater_season,3));
% freshwater_season= freshwater_season(:,:,2);
%  

FreshWater = (freshwaterc*86400*365)/(rho_I); %%% convert to GT, then normalize by densitoes
FreshWaterw = (freshwaterw*86400*365)/(rho_I); %%% convert to GT, then normalize by densitoes

for i = 1:Nx
    for j = 1:Ny
        FreshWater(i,j) = (FreshWater(i,j));  
        FreshWaterw(i,j) = (FreshWaterw(i,j));  
       
    end
end

  
FreshWater(bathy==0) = NaN;
FreshWater(FreshWater==0)=NaN;

FreshWaterw(bathy==0) = NaN;
FreshWaterw(FreshWaterw==0)=NaN;
%%% Set up the figure
figure(2);
clf;
scrsz = get(0,'ScreenSize');
subplot(2,1,1);
fontsize = 16;


  
  latMin = min(min(YC));
  latMax = max(max(YC));
  lonMin = min(min(XC));
  lonMax = max(max(XC));



        axesm('eqaconicstd',...
          'fontsize',12,...
          'Grid','on', ...    
          'Frame','off', ...
          'MapLatLimit',[-84 -63], ...
          'MapLonLimit',[-83 20], ... 
          'MapParallels',[-85 -65], ...
          'PLineLocation', 5, ...
          'MLineLocation', 10,...
          'MeridianLabel', 'on', ...
          'ParallelLabel', 'on');%, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
        axis off;
        setm(gca,'MLabelParallel',-20)

             
 hold on 
        
contourm(YC,XC,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)

        
hold on
        
contourfm(YC,XC,-FreshWater,60,'Edgecolor','none');
set(gca,'FontSize',14);

set(gca,'Position',[0.25 0.5 .52 .43]);
% handle(1) = colorbar;
colormap redblue(60);
% colormap(flipud(colormap));
caxis([-20 20])
% set(handle,'Position',[0.92 0.1 0.02 .8]);
        

xlabel('Longitude','interpreter','latex','FontSize',14);
ylabel('Latitude','interpreter','latex','FontSize',14);

      
h = title('Shelf Ice Melt (m/yr), Initial 34.45psu Cavity','interpreter','latex','Fontsize',15);
      set(h,'Position',[0 -.92 0]) 

hold on
ax2 = gca;
outerpos = ax2.OuterPosition;
ti = ax2.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax2_width = outerpos(3) - ti(1) - ti(3);
ax2_height = outerpos(4) - ti(2) - ti(4);
ax2.Position = [left bottom ax2_width ax2_height];
text(.3,-1.3,'a','fontsize',14,'interpreter','latex')


contourm(YC,XC,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)

        
hold on
subplot(2,1,2)  

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

contourfm(YC,XC,-FreshWaterw,60,'Edgecolor','none');
set(gca,'FontSize',14);

caxis([-20 20]);
set(gca,'Position',[0.25 0.01 .52 .43]);
handle(1) = colorbar;
colormap redblue(60);
% colormap(flipud(colormap));

caxis([-20 20])
set(handle,'Position',[0.85 0.1 0.03 .8]);
hold on
contourm(YC,XC,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
hold on
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom+.02 ax_width ax_height];       

xlabel('Longitude','interpreter','latex','FontSize',14);
ylabel('Latitude','interpreter','latex','FontSize',14);
text(.3,-1.3,'b','fontsize',14,'interpreter','latex')
     
h2 = title('Shelf Ice Melt (m/yr), Initial 34psu Cavity','interpreter','latex','Fontsize',15);
      set(h2,'Position',[0 -.93 0]) 
      