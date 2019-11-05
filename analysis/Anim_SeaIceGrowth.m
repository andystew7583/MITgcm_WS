% %%%%% Sea ice growth/melt
% So the net growth rate of ice concentration is SIdA (SIdA = SIdAbATO + SIdAbATC + SIdAbOCN), 
% which equals the (SIarea - SIareaPT)/deltaT. And the net growth rate of ice thickness is (SIdHbATO + SIdHbATC + SIdHbOCN + SIdHbFLO), 
% which equals the (SIheff - SIheffPT)/deltaT. 

%%%
%%% anim.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = false;

%%% Read experiment data
loadexp;


sp30days=  3456000;

%%% Call variables from SI area rate of change (PT)

 %%% area rate of change by open ocean atmospheric flux (m^2/m^2/s)



diagnum1 = 24;  %%% SIdAbATO
outfname1 =diag_fileNames{1,diagnum1};



 %%% area rate of change by atm flux over ice (m^2/m^2/s)

diagnum2 = 25;  %%% SIdAbATC
outfname2 =diag_fileNames{1,diagnum2};


 %%% area rate of change by ocean ice flux (m^2/m^2/s)

diagnum3 = 26;  %%% SIdAbOCN
outfname3 =diag_fileNames{1,diagnum3};


%%% Data index in the output data files
outfidx = 1;


%%% Specify color range
set_crange = 1;
crange = [-5 5]; %%%change in area


%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics
dumpFreq = abs(diag_frequency(diagnum));
nIter0 = 6480;
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((0:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);



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
handle = figure(20);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');
M = moviein(nDumps);

Amean = [];
tyears = [];

% for n = 29 
 for n = 1:length(dumpIters)
  dumpIters(n);
    
  t = dumpIters(n)*deltaT/86400/365;
  
  
  tyears(n) = t;
  %%% area rate of change by open ocean atmosphere flux (m^2/m^2/s)
  
  SIdAbATO = rdmdsWrapper(fullfile(exppath,'results',outfname1),dumpIters(n));
  if (isempty(SIdAbATO))
    error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
  end     
  
  SIdAbATO = SIdAbATO*sp30days;
  
  if (~isempty(find(isnan(SIdAbATO))))
    break
  end
  
  
   [idx1 idx2] = find(isnan(SIdAbATO));
   
   
  %%% area rate of change by atm flux over ice(m^2/m^2/s)

   
   
  SIdAbATC = rdmdsWrapper(fullfile(exppath,'results',outfname2),dumpIters(n));
  if (isempty(SIdAbATC))
    error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
  end
  
  SIdAbATC = SIdAbATC*sp30days;

  
 %%% area rate of change by ocean ice flux (m^2/m^2/s)

  
  SIdAbOCN = rdmdsWrapper(fullfile(exppath,'results',outfname3),dumpIters(n));
  if (isempty(SIdAbOCN))
    error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
  end  
  
  
  SIdAbOCN = SIdAbOCN*sp30days;

  
  
  %%% sum together to get total change in SIconc
  
  A = SIdAbATO + SIdAbATC + SIdAbOCN;
  
  
%   if (~isempty(idx1))
%     break;
%   end
  
%%% x/y plot

        FF = zeros(Nx,Ny);
        for i=1:Nx
          for j=1:Ny
            FF(i,j) = A(i,j);
          end
        end

    
%     FF(FF==0) = NaN;
FF(hFacC(:,:,1)==0)=NaN;
    
    max(max(FF))
    min(min(FF))
    
%     contourf(XC,YC,FF,100,'EdgeColor','None');  
    pcolor(XC,YC,FF);
    shading interp;        
    hold on;
%     contour(XC,YC,bathy,[-5000:1000:-1000 -500 -200],'EdgeColor','k');
    hold off;
    xlabel('x (km)');
    ylabel('y (km)');

    xlabel('Longitude','interpreter','latex');
    ylabel('Latitude','interpreter','latex');

    
  
    
  %%% Finish the plot
  handle=colorbar;
  colormap jet(200);
  set(handle,'FontSize',fontsize);
  title(['$t=',num2str(tyears(n)*365,'%.1f'),'$ days'],'interpreter','latex');
  if (set_crange)  
    caxis(crange);
  end
  
  set(gca,'Position',plotloc);
  set(gca,'FontSize',fontsize);
%   annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  M(n) = getframe(gcf);
  
end







