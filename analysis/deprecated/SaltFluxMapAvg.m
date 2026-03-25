%%%% Doing an average of the salt flux in winter, then plotting.  Plotting
%%%% bottom salinity

setExpname
loadexp;

ronneregion 

%%%%%%%has to be halfshelfice experiment

%%% Read experiment data


%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(38));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);


nDumps = 120;


nDump_start = 9*12;
nDump_finish = 18*12;

%%% Longitudinal positions

XC;


%%% Latitudinal positions

YC;


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

for i=1:Nx
    for j=1:Ny        
           
           if (((SHELFICEtopo(i,j)) - (bathy(i,j)) < 0  ||  (bathy(i,j) == SHELFICEtopo(i,j)) ))
               topog_msk(i,j) = 0;
           end
    end
end



%%%%Chose an iteration number to output (n)
SFLUX_tot = NaN(Nx,Ny,nDump_finish-nDump_start);
% SALT_tot = NaN(Nx,Ny,Nr,nDump_finish-nDump_start);
%%% n = 1:nDumps (but experiment hasn't finished yet %%%
for n=1:nDump_finish-nDump_start
  
  SFLUX = rdmdsWrapper(fullfile(exppath,'/results/SFLUX'),dumpIters(n+nDump_start));
  SFLUX_tot(:,:,n) = (SFLUX)/1000; %%%convert to kg/m2/s
  

end

%%%%getting the seasonal Salt Flux
Nyears = (nDump_finish-nDump_start)/12;

SFLUX_season = zeros(size(SFLUX_tot,1),size(SFLUX_tot,2),Nyears,4);

for n = 1:Nyears
  SFLUX_season(:,:,n,1) = mean(SFLUX_tot(:,:,(n-1)*12+[12 1 2]),3);
  SFLUX_season(:,:,n,2) = mean(SFLUX_tot(:,:,(n-1)*12+(3:5)),3);
  SFLUX_season(:,:,n,3) = mean(SFLUX_tot(:,:,(n-1)*12+(6:8)),3);
  SFLUX_season(:,:,n,4) = mean(SFLUX_tot(:,:,(n-1)*12+(9:11)),3);
  
end


SFLUX_seasonavg = squeeze(nanmean(SFLUX_season,3));
% SALT_seasonavg = squeeze(nanmean(SALT_season,4));

salt_final = NaN(Nx,Ny,4);
for i = 1:Nx
    for j = 1:Ny
        for k = 1:4
          SFLUX_seasonavg(i,j,k) = (SFLUX_seasonavg(i,j,k)*RAC(i,j)); %%% convert to kg 
%           salt_final(i,j,k) = SALT_seasonavg(i,j,kmax(i,j),k);
        end
    end
end


SFLUX_winter = SFLUX_seasonavg(:,:,3);
% SALT_winter = salt_final(:,:,3);

SFLUX_winter(bathy==0) = NaN;
% SFLUX_winter(bathy>=SHELFICEtopo)=NaN;
% SFLUX_winter(bathy==0) = NaN;
SFLUX_winter(SFLUX_winter==0) = NaN;   
       


% SALT_winter(bathy==0) = NaN;
% SALT_winter(bathy>=SHELFICEtopo)=NaN;
% SALT_winter(bathy==0) = NaN;
% SALT_winter(SFLUX_winter==0) = NaN;   
         
%%%get rid of northern boundary
SFLUX_winter(:,220:end)=NaN;

%%%% index the ronne polynya region (defined in ronneregion.m)
msk = NaN(Nx,Ny);
 for i=1:Nx
    for j=1:Ny
             if start(i,j)==1
                 msk(i,j)=1;
             end
             
                  
            
        
    end
  end

msk_shade =0*msk;
msk_alpha = .15*msk;



%%% Set up the figure
figure(2);
clf;
set(gcf,'Position',[ 477         225        1075         724]);
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



w=pcolorm(YC,XC,SFLUX_winter),shading interp;
colormap redblue(50);

% set(gca,'FontSize',10)
caxis([-1000 1000]);
hold on

c=colorbar;
set(c,'Position',[0.94 0.1 0.02 .75])

h = pcolorm(YC,XC,msk_shade);
set(h,'FaceAlpha','flat');
set(h,'AlphaDataMapping','none');
set(h,'AlphaData',msk_alpha);
set(h,'facecolor','k');
        
        
contourm(YC,XC,topog_msk,1,'EdgeColor','black','LineWidth',1)
% set(gca,'color',[.5 .5 .5]);


bathy(start==1)=NaN;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width-.05 ax_height];
t = title('Salt Flux in Winter (JJA) (kg)/s','interpreter','latex','Fontsize',18);
set(t,'Position',[0 -.93 0])
% ylim([-84 -66])
