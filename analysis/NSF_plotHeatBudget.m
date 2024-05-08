%%%
%%% NSF_plotHeatBudget.m
%%% 
%%% Plots eddy heat fluxes for our NSF proposal.
%%%

%%% Load experiment data
expdir = '../experiments';
% expname = 'hires_seq_onethird_notides_RTOPO2';
% expname = 'hires_seq_onesixth_notides_RTOPO2';
% expname = 'hires_seq_onetwelfth_notides_RTOPO2';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
deform_cavity = false;
% loadexp;

%%% Iteration number of model output to plot
% iter = 1949760 - 316800;
% iter = 1949760 - 100*86400/60;
iter = 1949760 - 50*86400/60;
% iter = 1949760 - 250*86400/60;

%%% MOC grid
ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
eta = -9:.1:11;
Neta = length(eta);

%%% T/S flux storage file
infname = [expname,'_TSfluxes'];
if (deform_cavity)
  infname = [infname,'_deform'];
end
infname = [infname,'.mat'];

%%% Load pre-computed vertical flux data
load(fullfile('products',infname),'wvel_tavg','theta_tavg','salt_tavg','wvelth_tavg','wvelslt_tavg','tflux_tavg','sflux_tavg');
tflux_tavg(sum(hFacC,3)==0) = NaN;
sflux_tavg(sum(hFacC,3)==0) = NaN;

%%% Calculate vertical fluxes
theta_f = -1.85;
salt_ref = 34.45;
theta_w = 0.5*(theta_tavg(:,:,1:Nr)+theta_tavg(:,:,[2:Nr 1]));
wvelth_ref = theta_f.*wvel_tavg;
wvelth_mean = theta_w.*wvel_tavg;
wvelth_eddy = wvelth_tavg - wvelth_mean;
salt_w = 0.5*(salt_tavg(:,:,1:Nr)+salt_tavg(:,:,[2:Nr 1]));
wvelslt_ref = salt_ref.*wvel_tavg;
wvelslt_mean = salt_w.*wvel_tavg;
wvelslt_eddy = wvelslt_tavg - wvelslt_mean;
theta_tavg(hFacC==0) = NaN;
clear('theta_w','salt_w','theta_tavg','salt_tavg'); %%% Free up mmeory

%%% Define grid to extract data sections
startLat = -82.3;
endLat = -72;
startLon = -62;
endLon = -40;
Nsec = 2001;
dLat = (endLat-startLat)/(Nsec-1);
dLon = (endLon-startLon)/(Nsec-1);
secLats = startLat:dLat:endLat;
secLons = startLon:dLon:endLon;
% plotm(secLats,secLons,'k-');

%%% Extract data along defined sections
S = rdmdsWrapper(fullfile(exppath,'results','SALT_12hourly'),iter);
T = rdmdsWrapper(fullfile(exppath,'results','THETA_12hourly'),iter);
load(fullfile('products',[expname,'_EKE.mat']),'EKE');
S(hFacC==0) = NaN;
T(hFacC==0) = NaN;
EKE(hFacC==0) = NaN;
secS = zeros(Nsec,Nr);
secT = zeros(Nsec,Nr);
secE = zeros(Nsec,Nr);
secB = zeros(Nsec,1);
secI = zeros(Nsec,1);
secH = zeros(Nsec,Nr);
for n=1:Nsec
  jm = find(secLats(n)>=YC(1,:),1,'last');
  im = find(secLons(n)>=XC(:,1),1,'last');
  jp = jm+1;
  ip = im+1;
  wp_y = (secLats(n)-YC(1,jm))/(YC(1,jp)-YC(1,jm));
  wm_y = 1-wp_y;
  wp_x = (secLons(n)-XC(im,1))/(XC(ip,1)-XC(im,1));
  wm_x = 1-wp_x;
%   secS(n,:) = squeeze(wm_x*wm_y*S(im,jm,:) + wp_x*wm_y*S(ip,jm,:) + wm_x*wp_y*S(im,jp,:) + wp_x*wp_y*S(ip,jp,:));
%   secT(n,:) = squeeze(wm_x*wm_y*T(im,jm,:) + wp_x*wm_y*T(ip,jm,:) + wm_x*wp_y*T(im,jp,:) + wp_x*wp_y*T(ip,jp,:));
%   secH(n,:) = squeeze(wm_x*wm_y*hFacC(im,jm,:) + wp_x*wm_y*hFacC(ip,jm,:) + wm_x*wp_y*hFacC(im,jp,:) + wp_x*wp_y*hFacC(ip,jp,:)); 
  if (wp_x > wm_x)
    xidx = ip;
  else
    xidx = im;
  end
  if (wp_y > wm_y)
    yidx = jp;
  else
    yidx = jm;
  end
  secS(n,:) = squeeze(S(xidx,yidx,:));
  secT(n,:) = squeeze(T(xidx,yidx,:));
  secE(n,:) = squeeze(EKE(xidx,yidx,:));
  secH(n,:) = squeeze(hFacC(xidx,yidx,:));
  secB(n) = bathy(xidx,yidx);
  secI(n) = SHELFICEtopo(xidx,yidx);
end











                
%%% Plotting options                
fontsize = 14;
framepos = [417   100   986   880];
labelspacing = 200;
axpos = zeros(2,4);
axpos(1,:) = [0.06 0.78 0.87 0.18];
axpos(2,:) = [0.06 0.53 0.87 0.18];
axpos(3,:) = [0.06 0.05 0.44 0.4];
axpos(4,:) = [0.63 0.05 0.3 0.4];
cb1_pos = [0.95 0.78 0.015 0.18];
cb2_pos = [0.95 0.53 0.015 0.18];
cb3_pos = [0.52 0.05 0.015 0.4];
axlabels = {'(a)','(b)','(c)','(d)'};
icecolor = [186 242 239]/255;



figure(205);
clf;
set(gcf,'Position',framepos);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TEMP AND EKE PANELS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ZZ,LA] = meshgrid(RC,secLats);
for n=1:Nsec
  kmin = find(secH(n,:)>0,1,'first');
  kmax = find(secH(n,:)>0,1,'last');
  if (isempty(kmin) || isempty(kmax))
    continue;
  end
  if (kmin == kmax)
    secT(n,kmin) = NaN;
    secS(n,kmin) = NaN;
    continue;
  end
  ZZ(n,kmin) = RF(kmin);% - (1-sechFac(n,kmin))*DRF(kmin);
  ZZ(n,kmax) = RF(kmax+1);% - (1-sechFac(n,kmax))*DRF(kmax);
  
end


subplot('Position',axpos(1,:));
pcolor(LA,-ZZ,secT);
set(gca,'YDir','reverse');
shading interp;
hold on
[C,h] = contour(LA,-ZZ,secS,[34.1:.1:35],'EdgeColor','k');
clabel(C,h,'FontSize',12);
plot(secLats,-secB,'k-','LineWidth',3);
Iidx = find(secI==0,1,'first');
plot(secLats(1:Iidx),-secI(1:Iidx),'k-','LineWidth',3);
hold off;
caxis([-2.5 -1]);
set(gca,'YLim',[0 1400]);
set(gca,'XLim',[startLat (-73)]);
colormap(gca,cmocean('thermal',15))
cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
title(cbhandle,'$^\circ$C','Fontsize',14,'interpreter','latex');
set(gca,'FontSize',fontsize);
% xlabel('Latitude')
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);
title(['Instantaneous potential temperature, salinity, ',datestr(datenum('01-Jan-2008')+iter*deltaT/86400)]);

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
msk = ones(Nsec,Nr);
msk(ZZ<repmat(secI,[1 Nr])) = NaN;
pcolor(LA,-ZZ,msk);
shading interp;
colormap(ax2,icecolor);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
set(ax2,'YLim',[0 1400]);
set(ax2,'YDir','reverse');
set(ax2,'XLim',[startLat (-73)]);
set(ax2,'Color','None')


subplot('Position',axpos(2,:));
pcolor(LA,-ZZ,sqrt(2*secE));
set(gca,'YDir','reverse');
shading interp;
hold on
[C,h] = contour(LA,-ZZ,secS,[34.1:.1:35],'EdgeColor','k');
clabel(C,h,'FontSize',12);
plot(secLats,-secB,'k-','LineWidth',3);
Iidx = find(secI==0,1,'first');
plot(secLats(1:Iidx),-secI(1:Iidx),'k-','LineWidth',3);
hold off;
% caxis([-4 -2]);
caxis([0 0.2]);
set(gca,'YLim',[0 1400]);
set(gca,'XLim',[startLat (-73)]);
colormap(gca,cmocean('amp',15))
cbhandle = colorbar;
set(cbhandle,'Position',cb2_pos);
set(gca,'FontSize',fontsize);
xlabel('Latitude')
ylabel('Depth (m)');
set(gca,'Color',[.8 .8 .8]);
title(cbhandle,'m/s','Fontsize',14,'interpreter','latex');
title('Root-mean-square eddy velocity, 2009-2014 mean');

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
msk = ones(Nsec,Nr);
msk(ZZ<repmat(secI,[1 Nr])) = NaN;
pcolor(LA,-ZZ,msk);
shading interp;
colormap(ax2,icecolor);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
set(ax2,'YLim',[0 1400]);
set(ax2,'YDir','reverse');
set(ax2,'XLim',[startLat (-73)]);
set(ax2,'Color','None')







%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HEAT FLUX PANELS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


kidx = 23;
msk = real((SHELFICEtopo>=0) & (ETA<3.5) & (bathy<0) & (bathy<zz(kidx)));
plotmsk = real(msk);
plotmsk(plotmsk==0) = NaN;
bathycntrs = [0 250 500 1000 2000 3000 4000 5000];
bathy_plot = bathy;
bathy_plot(hFacC(:,:,kidx)==0) = NaN;

iceshelf_msk = real((SHELFICEtopo<0) & (SHELFICEtopo>bathy));
iceshelf_msk(iceshelf_msk==0) = NaN;

wvelth_eddy_plot = wvelth_eddy(:,:,23);
wvelth_eddy_plot(sum(hFacC,3)==0) = NaN;
subplot('Position',axpos(3,:));
pcolor(XC,YC,wvelth_eddy_plot*rho0*Cp);
shading interp;
hold on;
[C,h] = contour(XC,YC,-bathy_plot,bathycntrs,'EdgeColor','k');
clabel(C,h);
hold off;
cbhandle = colorbar;
set(cbhandle,'Position',cb3_pos);
title(cbhandle,'W/m^2')
colormap(gca,cmocean('balance',40));
caxis([-100 100]);
set(gca,'XLim',[-63 -25]);
set(gca,'YLim',[-78.5 -71.3]);
set(gca,'FontSize',fontsize);
xlabel('Longitude');
ylabel('Latitude');
% set(gca,'Color',[.8 .8 .8]);
set(gca,'Color','w');
title('Vertical eddy heat flux, 2009-2014 mean');

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
% hold on;
handle = pcolor(XC,YC,-real(msk));
shading interp;
colormap(ax2,hot);
caxis(ax2,[-1 0]);
set(handle,'FaceAlpha',0.1);
set(ax2,'Color','None');
set(ax2,'XLim',[-63 -25]);
set(ax2,'YLim',[-78.5 -71.3]);
set(ax2,'XTick',[]);
set(ax2,'YTick',[]);
ax3 = axes('Position',get(ax1,'Position'));
handle = pcolor(XC,YC,real(iceshelf_msk));
shading interp;
hold on;
plot(secLons,secLats,'k:','LineWidth',2);
hold off
colormap(ax3,icecolor);
set(handle,'FaceAlpha',1);
set(ax3,'Color','None');
set(ax3,'XLim',[-63 -25]);
set(ax3,'YLim',[-78.5 -71.3]);
set(ax3,'XTick',[]);
set(ax3,'YTick',[]);

sum(sum(rho0.*Cp.*wvelth_eddy(:,:,23).*msk.*RAC)) / sum(sum(msk.*RAC))
nansum(nansum(-tflux_tavg(:,:).*msk.*RAC)) / sum(sum(msk.*RAC))

subplot('Position',axpos(4,:));
plot(squeeze(sum(sum((wvelth_tavg-wvelth_ref).*RAC.*msk)))*rho0*Cp/sum(sum(msk.*RAC)),-squeeze(zz),'LineWidth',1.5);
hold on;
plot(squeeze(sum(sum((wvelth_mean-wvelth_ref).*RAC.*msk)))*rho0*Cp/sum(sum(msk.*RAC)),-squeeze(zz),'LineWidth',1.5);
plot(squeeze(sum(sum(wvelth_eddy.*RAC.*msk)))*rho0*Cp/sum(sum(msk.*RAC)),-squeeze(zz),'LineWidth',1.5);
hold off
set(gca,'XLim',[-2 7]);
xlabel('Vertical heat flux (W/m^2)')
ylabel('Depth (m)');
set(gca,'YDir','reverse');
legend('Total','Mean','Eddy','Location','SouthEast');
set(gca,'YLim',[0 1000]);
% title(['Ice shelf face to \eta=3.5, total surface heat flux = ',num2str(-nansum(nansum(tflux_tavg.*RAC.*msk/1e12))),' TW']);
set(gca,'FontSize',fontsize);
box off;

ax2 = axes('Position',get(gca,'Position'));
plot([-2 7]*sum(sum(msk.*RAC))/1e12,-[zz(kidx) zz(kidx)],'k--','LineWidth',1);
hold on;
plot([0 0],[0 1000],'k:','LineWidth',1);
hold off;
set(ax2,'XLim',[-2 7]*sum(sum(msk.*RAC))/1e12);
set(ax2,'XAxisLocation','top');
set(ax2,'YAxisLocation','right');
set(ax2,'YLim',[0 1000]);
set(ax2,'YTick',[]);
box off;
set(ax2,'Color','None');
set(ax2,'FontSize',fontsize);
set(ax2,'YDir','reverse');
xlabel(ax2,'Vertical heat flux (TW)');

% Create textbox
annotation(gcf,'textbox',...
  [0.0149290060851927 0.709772727272727 0.05 0.05],'String','(a)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.0152332657200811 0.455227272727272 0.05 0.0500000000000002],...
  'String','(b)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.0101622718052738 0 0.05 0.0300000000000001],...
  'String','(c)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
  [0.565943204868154 0 0.0499999999999999 0.0300000000000001],...
  'String','(d)',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',16,...
  'FitBoxToText','off');


% Create textbox
annotation(gcf,'textbox',...
  [0.298174442190669 0.322295454545454 0.0588235294117647 0.0255681818181818],...
  'String',{'Transect'},'EdgeColor','None');

