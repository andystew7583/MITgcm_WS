
%%%
%%% animatecon.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie(in conical form) of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% Read experiment data
loadexp;

%%% Set true if plotting on a Mac
mac_plots = false;

%%% Select diagnostic variable to animate
diagnum = 3;
outfname = diag_fileNames{1,diagnum};

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 15;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 0;

%%% Set true to plot the field in the topmost wet cell at each horizontal
%%% location
topplot = 0;

%%% Set true to plot the field in the middle of the water column at each
%%% horizontal location
midplot = 0;


%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics
dumpFreq = abs(diag_frequency(diagnum));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((0:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);



%%% Longitudinal positions

XC;


%%% Latitudinal positions

YC;

%%% Mesh grids for plotting
kmax = ones(Nx,Ny);
kmin = ones(Nx,Ny);
kmid = ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    idx = find(squeeze(hFacC(i,j,:))>0);
    if (~isempty(idx))
      kmin(i,j) = min(idx);
      kmax(i,j) = max(idx);
      kmid(i,j) = round((kmin(i,j) + kmax(i,j)) / 2);
    end
  end
end



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



%%% Set up the figure
handle = figure(1);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');



dy = dumpIters;

M = moviein(length(nDumps));
for n=1:1:length(dy) 
    
    t = dumpIters(n)*deltaT/86400/365;
    tyears(n) = t;
    A = rdmdsWrapper(fullfile(exppath,'results',outfname),dy(n));          
    if (isempty(A))
      error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
    end 
    
    
  [idx1 idx2] = find(isnan(A))
  
  if (~isempty(idx1))
    break;
  end
    

  FF = squeeze(A(:,:,xylayer,outfidx));        
  FF(hFacC(:,:,xylayer)==0) = NaN;
  

  
        
        figure(1);
        clf
        set(gcf,'Color','w');
%         quikpcolor(xc,yc,tmp')
  latMin = min(min(YC));
  latMax = max(max(YC));
  lonMin = min(min(XC));
  lonMax = max(max(XC));

% 'FLatLimit', [latMin latMax], ...
%           'FLonLimit', [lonMin lonMax], ...

        axesm('eqaconicstd',...
          'Grid','on', ...    
          'Frame','on', ...
          'MapLatLimit',[-84 -63], ...
          'MapLonLimit',[-83 -3], ... 
          'MapParallels',[-85 -65], ...
          'PLineLocation', 5, ...
          'MLineLocation', 10,...
          'MeridianLabel', 'on', ...
          'ParallelLabel', 'on');%, ...
%           'origin',[yc(round(size(yc,1)/2),round(size(yc,2)/2)),xc(round(size(xc,1)/2),round(size(xc,2)/2))])
        axis off;
        setm(gca,'MLabelParallel',-20)
%         pcolor(xc,yc,tmp);

        FF(FF==0) = NaN;
        
        
%%%%%%%% Get rid of BL lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        FF(bathy==0)=NaN;
        for i = 1:Nx
            for j = 1:Ny
                if (SHELFICEtopo(i,j)-bathy(i,j) == 0)
                    bathy(i,j) = NaN;
                end
            end
        end
        
       

        FF(end,:,:) = NaN;
        
             
%         setm(gca,'mlabellocation',.2,'plabellocation',.2,'mlinelocation',.2,'plinelocation',.2)
%         tightmap;


        shading interp;
%         caxis(cx)
        caxis([-1 1])
        set(gca,'FontSize',14);
        set(gca,'Position',[0.03 0 .91 .95]);
        handle = colorbar;
        colormap jet(160);
        
        set(handle,'FontSize',8);
        set(handle,'Position',[0.92 0.1 0.02 .85]);

        FF(FF==0) = NaN;

        pcolorm(YC,XC,FF);        
%         setm(gca,'mlabellocation',.2,'plabellocation',.2,'mlinelocation',.2,'plinelocation',.2)
%         tightmap;
        xlabel('Longitude','FontSize',18);
        ylabel('Latitude','FontSize',18);
        shading interp;
%         caxis(cx)
        set(gca,'FontSize',18);
        set(gca,'Position',[-0.05 0 1 1]);

%         set(gca,'XTick',[-80:10:0]);
%         set(gca,'YTick',[-80:5:-60]);
%         axesm ('eqaazim', 'Frame', 'on', 'Grid', 'on');
       
       title(['$t=',num2str(tyears(n)*365,'%.1f'),'$ days'],'interpreter','latex');

        pause(.01)
        M(n) = getframe(gcf);
 end

