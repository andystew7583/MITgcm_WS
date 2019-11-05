%%%
%%% plotRd.m
%%%
%%% Calculates and plots the first Rossby radius of deformation as a 
%%% function of latitude, calculated from the time/zonal mean
%%% temperature and salinity output.
%%%

%%% Calculate time and zonal average
avg_t;
avg_xt;

%%% Plotting options
mac_plots = 0;
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
  fontsize = 26;
else
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.9 scrsz(4)/2];
  fontsize = 20;
end

%%% Temperature stratification   
TT = tt_avg;
DZ = repmat(delR,[Ny 1]);
Tz = zeros(size(TT));  
Tz(:,1) = (TT(:,1)-TT(:,2)) ./ (0.5*(DZ(:,1)+DZ(:,2)));
Tz(:,2:Nr-1) = (TT(:,1:Nr-2)-TT(:,3:Nr)) ./ (DZ(:,2:Nr-1) + 0.5*(DZ(:,1:Nr-2)+DZ(:,3:Nr)));
Tz(:,Nr) = (TT(:,Nr-1)-TT(:,Nr)) ./ (0.5*(DZ(:,Nr-1)+DZ(:,Nr)));  
N2 = tAlpha*gravity*Tz;

%%% Calculate Rd
Rd = zeros(Ny,1);
for j=1:Ny    
  Cig = 0;
  for k=1:Nr
    if (TT(j,k)~=0)
      Cig = Cig + sqrt(N2(j,k))*delR(k);      
    end    
  end  
  Rd(j) = Cig/(pi*abs(f0));
end

%%% Create a topography-following grid
[ZZ YY] = meshgrid(zz,yy);
for j=1:Ny
  kmax = sum(squeeze(hFacC(1,j,:))~=0);
  if (kmax > 0)
    ZZ(j,kmax) = bathy(1,j);
    N2(j,kmax) = N2(j,kmax-1); % + (ZZ(j,kmax)-ZZ(j,kmax-1))*(N2(j,kmax-1)-N2(j,kmax-2))/(ZZ(j,kmax-1)-ZZ(j,kmax-2));
  end
end

%%% Plot the result
figure(5);
clf;
axes('FontSize',16);
plot(yy,Rd,'k-');
xlabel('y (km)');
ylabel('Rd (m)');
  
%%% Plot stratification
handle = figure(6);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
contourf(YY/1000,ZZ/1000,log10(N2),[-9:0.1:-4],'EdgeColor','None');
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\mathrm{log}_{10}N^2$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'clim',[-9 -4]);