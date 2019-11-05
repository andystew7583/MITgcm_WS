%%%%isolating FRIS cavity region

%%%% Doing an average of the salt flux in winter, then plotting.  Plotting
%%%% bottom salinity

setExpname
loadexp;


%%%%%%%has to be halfshelfice experiment

%%% Read experiment data


%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(38));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
deltaT_4 = 200;
nIter0_4 = 1;
nDumps_4 = round(nTimeSteps*deltaT_4/dumpFreq);
dumpIters4 = round((1:nDumps_4)*dumpFreq/deltaT_4); 
dumpIters4 = dumpIters4(dumpIters4 >= nIter0_4);






%%% Longitudinal positions

XC;


%%% Latitudinal positions

YC;


%%%% index the FRIS area (-80 to -20) W, 75S -> !
msk = NaN(Nx,Ny);
  for i=1:Nx
    for j=1:Ny
             if (XC(i,1))<-30
                 if YC(i,j)<-75
                    if SHELFICEtopo(i,j)<0 &&SHELFICEtopo(i,j)-bathy(i,j) >0
                        
                        msk(i,j)=1;

                    end
                 end
             end
                  
            
        
    end
  end

msk_shade =0*msk;
msk_alpha = 60*msk;




%%% finding minimum k level for grid

kmax = ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    idx = find(squeeze(hFacC(i,j,:))>0);
    if (~isempty(idx))
      kmax(i,j) = max(idx);
    end
  end
end


topog_msk = ones(Nx,Ny);
hold on
for i=1:Nx
    for j=1:Ny        
           
           if (((SHELFICEtopo(i,j)) - (bathy(i,j)) < 0  ||  (bathy(i,j) == SHELFICEtopo(i,j)) ))
               topog_msk(i,j) = 0;
               bathy(i,j)=NaN;
           end
    end
end

%%% Set up the figure
figure(2);
clf;
scrsz = get(0,'ScreenSize');
fontsize = 14;


  
  latMin = min(min(YC));
  latMax = max(max(YC));
  lonMin = min(min(XC));
  lonMax = max(max(XC));



        axesm('eqaconicstd',...
          'fontsize',16,...
          'Grid','on', ...    
          'Frame','on', ...
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

w=pcolorm(YC,XC,bathy),shading interp;
colormap haxby(50);

hold on

c=colorbar;
set(c,'Position',[0.93 0.1 0.02 .75])
title(c,'Depth(m)','interpreter','latex','fontsize',16);
hold on;
h = pcolorm(YC,XC,msk_shade);
set(h,'FaceAlpha','flat');
set(h,'AlphaDataMapping','direct');
set(h,'AlphaData',msk_alpha);
set(h,'facecolor','k');
% shading flat;
hold off;
 
        
        
% contourm(YC,XC,topog_msk,1,'EdgeColor','black','LineWidth',1)
% set(gca,'color',[.5 .5 .5]);


ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width-.05 ax_height];
h = title('FRIS Cavity Region','interpreter','latex','Fontsize',18);
set(h,'Position',[0 -.93 0])
% ylim([-84 -66])