% Validate_Cryosat , polar-sterographic projection


%%%% Grid to interpolate to....
run ../../newexp/defineGrid.m
%%%%% get our grid

XMC = XMC';
YMC = YMC';


datafile = 'CS2.nc';
DOT = ncread(datafile,'DOT');
lat = ncread(datafile,'Latitude');
lon = ncread(datafile,'Longitude');
MDT = ncread(datafile,'MDT');


Nyears = 6;
Nmonths = 12;
C_x = size(lon,2);
C_y = size(lat,1);

time_total = Nyears*Nmonths;


%%% transform our lat/lon meshgrids


%%%%% load datadirectory
interpSSH = NaN(Nx,Ny,time_total);

n=0;
for n = 1:time_total
      
        %%% Print a message to keep track of progress
        ['year = ',num2str(n)]
      
        %%% Velocity data file
                Old_DOT = DOT(:,:,n);
                LA = reshape(lat,C_y*C_x,1);
                LO = reshape(lon,C_y*C_x,1); 
                V = Old_DOT;
                V = reshape(V,C_y*C_x,1);
                F = scatteredInterpolant(double(LO),double(LA),double(V),'natural');

                % Store the velocity data in the storage matrix
                
                interpSSH(:,:,n) = F(double(XMC),double(YMC));

            
            
        
    
end
% 
% 
save SSHvalidation interpSSH

%%%%% Now we validate
load SSHvalidation



%%% Read experiment data
run ../setExpname.m;
run ../loadexp;
%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(7));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);
g=9.81;

start_yr = 13;
end_yr = 19;
Nyrs = end_yr-start_yr;

Start_dump = start_yr*12;
end_dump = end_yr*12;



%%%%Chose an arbitrary time to output 
SSH_DOT = NaN(Nx,Ny,end_dump-Start_dump);

%%% n = 1:nDumps (but experiment hasn't finished yet %%%
for n=Start_dump:end_dump-1  
  
  SSH = rdmdsWrapper(fullfile(exppath,'/results/PHIHYD'),dumpIters(n));
  SSH_DOT(:,:,n-(Start_dump-1)) = (SSH(:,:,1)*100)./g;   %%%convert to cm
  
  tt(n) =  (dumpIters(n)*deltaT)/86400;
  %tt(n) is the day of the run
  
end

Seasonal_SSH = zeros(size(SSH_DOT,1),size(SSH_DOT,2),Nyrs,4);

Seasonal_cry = zeros(size(SSH_DOT,1),size(SSH_DOT,2),Nyrs,4);
Nyrs_2=6;
for n=1:Nyrs
  Seasonal_SSH(:,:,n,1) = mean(SSH_DOT(:,:,(n-1)*12+(2:4)),3);
  Seasonal_SSH(:,:,n,2) = mean(SSH_DOT(:,:,(n-1)*12+(5:7)),3);
  Seasonal_SSH(:,:,n,3) = mean(SSH_DOT(:,:,(n-1)*12+(8:10)),3);
  Seasonal_SSH(:,:,n,4) = mean(SSH_DOT(:,:,(n-1)*12+([11 12 1])),3);
end
for n = 1:Nyrs_2
  %%%seasonal avg of cryosat-2
  Seasonal_cry(:,:,n,1) = mean(interpSSH(:,:,(n-1)*12+([12 1 2])),3);
  Seasonal_cry(:,:,n,2) = mean(interpSSH(:,:,(n-1)*12+(3:5)),3);
 
  Seasonal_cry(:,:,n,3) = mean(interpSSH(:,:,(n-1)*12+(6:8)),3);
  Seasonal_cry(:,:,n,4) = mean(interpSSH(:,:,(n-1)*12+(9:11)),3);
  
  
end


SSH = mean(Seasonal_SSH,3);

Cry = mean(Seasonal_cry,3);


%%%%NaNing out shelf ice
SSH_DJF = SSH(:,:,1);

SSH_SON = SSH(:,:,4);


SSH_JJA = SSH(:,:,3);

Cry_DJF = Cry(:,:,1);



Cry_SON = Cry(:,:,4);


Cry_JJA = Cry(:,:,3);


%%%%%remove the mean difference between obs-model

SSH_DJF(SHELFICEtopo-bathy<=0)=NaN;
SSH_JJA(SHELFICEtopo-bathy<=0)=NaN;
SSH_JJA(SHELFICEtopo<0)=NaN;
SSH_DJF(SHELFICEtopo<0)=NaN;
SSH_DJF(:,224:end)=NaN;
SSH_JJA(:,224:end)=NaN;
Cry_DJF(:,224:end)=NaN;
Cry_JJA(:,224:end)=NaN;


diff_djf = (Cry_DJF-SSH_DJF);
diff_djf_full = Cry_DJF-SSH_DJF;
diff_djf_new = diff_djf_full-diff_djf;


figure(1)
clf
subplot(2,1,1);



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

pcolorm(YG,XG,SSH_DJF),shading interp
set(gca,'Position',[0.2 0.48 .6 .45]);
% h=colorbar;
% set(h,'Position',[0.8 0.47 0.02 .4])
caxis([-220 -180]);
colormap((brewermap(40,'Spectral')));
h=title('DJF Dynamic SSH Anomaly (cm)','Interpreter','latex','FontSize',17);
set(h,'Position',[0 -.94 0])
xlabel('Longitude','interpreter','latex','FontSize',12);
ylabel('Latitude','interpreter','latex','FontSize',12);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom+.02 ax_width+.02 ax_height-.1];



diff_jja =(Cry_JJA-SSH_JJA);
diff_jja_full = Cry_JJA-SSH_JJA;
diff_jja_new = diff_jja_full-diff_jja;
hold on
subplot(2,1,2);


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



pcolorm(YG,XG,SSH_JJA),shading interp
colormap((brewermap(40,'Spectral')));
o=colorbar;
caxis([-220 -180]);
set(o,'Position',[0.83 0.02 0.02 .86])
h=title('JJA Dynamic SSH Anomaly (cm)','Interpreter','latex','FontSize',17);
set(h,'Position',[0 -.94 0])
xlabel('Longitude','interpreter','latex','FontSize',12);
ylabel('Latitude','interpreter','latex','FontSize',12);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom-.05 ax_width ax_height];


