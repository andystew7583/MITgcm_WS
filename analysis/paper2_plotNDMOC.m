%%%
%%% paper2_plotNDMOC.m
%%%
%%% Plots the overturning circulation in neutral density space.
%%%


%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onesixth_RTOPO2';
loadexp;

%%% Options (see calcOverturning)
calc_psi_eddy = true;
deform_cavity = false;
use_layers = false;
densvar = 'ND1';
psimax = 6;
psistep = 0.25;
% psimax = 2;
% psistep = 0.1;
% ylim = [27.3 28.2];
ylim = [27.3 29];

%%% Construct output file name
outfname = [expname,'_MOC_',densvar];
if (calc_psi_eddy)
  if (use_layers)
    estr = '_layers';
  else
    estr = '_TRM';
  end
else
  estr = '_noeddy';
end
outfname = [outfname,estr];
if (deform_cavity)
  outfname = [outfname,'_deform'];
end
outfname = [outfname,'.mat'];

%%% Load MOC data file
load(fullfile('products',outfname));
psi_plot = psi_mean_stand/1e6;
% psi_plot = (psi_mean_stand+mean(psi_mean_fluc,3))/1e6;
% psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);

%%% For plotting eta grid
bathy_plot = bathy;
bathy_plot(SHELFICEtopo-bathy<=0) = NaN;
ETA_plot = ETA;
ETA_plot(SHELFICEtopo-bathy<=0) = NaN;

%%% For plotting neutral dens
idx = find(XC(:,1)<-38,1,'last');
dens_min = 27.3;
dens_max = 29.0;
dens_plot = dens_tavg;
dens_plot(hFacC==0)= NaN;
dens_plot = squeeze(dens_plot(idx,:,:));
dens_plot(dens_plot<dens_min) = dens_min;
dens_plot(dens_plot>dens_max) = dens_max;

%%% Load pre-computed neutral density variables and compute potential
%%% density
load(fullfile('products',[expname,'_ND1.mat']));
ND1_ref = gg_ref;
PD0 = densjmd95(ss_ref,pt_ref,-gravity*rhoConst*RC(1)/1e4*ones(Nx,Ny,Nr))-1000;

%%% Mask for region in which JM97 neutral density variable is well defined
msk = (YC>-80) & (XC>-64); 
msk = repmat(msk,[1 1 Nr]);
msk = msk*1.0;
msk(msk==0) = NaN;
ND1_ref_msk = ND1_ref.*msk;

%%% Compute stratifications
pp_mid = 0.5*(pp_ref(:,:,1:Nr-1)+pp_ref(:,:,2:Nr));
DRC = repmat(RC(1:Nr-1)-RC(2:Nr),[Nx Ny 1]);
Nsq_true = - gravity/rhoConst * (...
                densjmd95(ss_ref(:,:,1:Nr-1),pt_ref(:,:,1:Nr-1),pp_mid) ...
              - densjmd95(ss_ref(:,:,2:Nr),pt_ref(:,:,2:Nr),pp_mid) ) ...
            ./ DRC;
Nsq_ND1 = -gravity/rhoConst * (ND1_ref(:,:,1:Nr-1)-ND1_ref(:,:,2:Nr)) ./DRC; 

%%% Calculate neutral density-potential density volume bins
NDmin = 26.5;
NDmax = 29;
NDstep = 0.01;
PD0min = 26.9;
PD0max = 28.1;
PD0step = 0.01;
ND1PD0vol = binByVolume(ND1_ref,PD0,ones(size(ND1_ref)), ...
                  NDmin,NDmax,NDstep,PD0min,PD0max,PD0step, ...
                  RAC,DRF,hFacC);
ND1PD0vol(ND1PD0vol==0) = NaN;   
logNsq_min = -11;
logNsq_max = -3;
logNsq_step = 0.05;
Nsq1NsqtrueVol = binByVolume(log10(Nsq_ND1),log10(Nsq_true),ones(size(Nsq_ND1)), ...
                  logNsq_min,logNsq_max,logNsq_step,logNsq_min,logNsq_max,logNsq_step, ...
                  RAC,DRC,ones(Nx,Ny,Nr-1));
Nsq1NsqtrueVol(Nsq1NsqtrueVol==0) = NaN; 
[PPP,NNN] = meshgrid(PD0min:PD0step:PD0max,NDmin:NDstep:NDmax);
[logNsq1,logNsqtrue] = meshgrid(logNsq_min:logNsq_step:logNsq_max,logNsq_min:logNsq_step:logNsq_max);










%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE THE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
fontsize = 14;
framepos = [417    26   791   959];
labelspacing = 200;
axpos = zeros(6,4);
axpos(1,:) = [0.08 0.72 0.84 0.25];
axpos(2,:) = [0.08 0.40 0.37 0.25];
axpos(3,:) = [0.55 0.40 0.37 0.25];
axpos(4,:) = [0.08 0.05 0.84 0.25];
cb1_pos = [0.93 0.72 0.02 0.25];
cb2_pos = [0.93 0.40 0.02 0.25];
cb3_pos = [0.93 0.05 0.02 0.25];
axlabels = {'(a)','(b)','(c)','(d)'};

figure(215);                
clf;
set(gcf,'Position',framepos);                         

%%% Neutral density
subplot('Position',axpos(1,:));
idx = find(XC(:,1)<-38,1,'last');
[ZZ,YY]=meshgrid(RC,YC(1,:));
set(gcf,'Position',framepos);
contourf(YY,-ZZ,dens_plot,[dens_min:.05:dens_max],'EdgeColor','None');
hold on
plot(YY(:,1),-bathy(idx,:),'k-','LineWidth',2);
plot(YY(:,1),-SHELFICEtopo(idx,:),'k-','LineWidth',2);
hold off;
colormap(gca,cmocean('dense',round((dens_max-dens_min)/0.05)-1));
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
set(gca,'YDir','reverse');
set(gca,'Color',[.85 .85 .85]);
axis([-81 YC(1,end-spongethickness) 0 4000]);
cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
title(cbhandle,'kg/m^3');
text(-80.5,3700,'Neutral density, \gamma, at 38W','FontSize',fontsize+2);

%%% ND1 stratification vs actual stratification
subplot('Position',axpos(2,:));
% pcolor(NNN,PPP,log10(ND2PD0vol));
pcolor(logNsq1,logNsqtrue,log10(Nsq1NsqtrueVol));
shading flat;
hold on
plot(logNsq1(1,:),logNsqtrue(:,1)+log10(2),'k--');
colormap(gca,cmocean('amp'));
caxis([10 14.5]);
axis([-8 -3.2 -8 -3.2]);
set(gca,'FontSize',fontsize);
set(gca,'XTick',[-8:1:-3]);
set(gca,'YTick',[-8:1:-3]);
set(gca,'XTickLabel',{'10^-^8','10^-^7','10^-^6','10^-^5','10^-^4','10^-^3'});
set(gca,'YTickLabel',{'10^-^8','10^-^7','10^-^6','10^-^5','10^-^4','10^-^3'});
xlabel('Neutral Density stratification, -(g/\rho_0) \partial\gamma/\partial z');
ylabel('Buoyancy frequency, N^2');
text(-5.8,-3.55,'$N^2= - \frac{1}{2} \frac{g}{\rho_0} \frac{\partial\gamma}{\partial z}$','interpreter','latex','FontSize',fontsize+4);


%%% ND1 vs PD0
subplot('Position',axpos(3,:));
pcolor(NNN,PPP,log10(ND1PD0vol));
shading flat;
colormap(gca,cmocean('amp'));
caxis([10 14.5]);
set(gca,'FontSize',fontsize);
axis([27.2 28.9 27.25 28.1]);
xlabel('Neutral Density, \gamma (kg/m^3)');
ylabel('Potential Density, \sigma_\theta (kg/m^3)');


%%% Add colorbar for scatter plots
cbhandle = colorbar;
set(cbhandle,'Position',cb2_pos);
title(cbhandle,'m^3')
set(cbhandle,'YTick',[10:1:14]);
set(cbhandle,'YTickLabel',{'10^1^0','10^1^1','10^1^2','10^1^3','10^1^4'});


%%% Streamfunction
psi_plot(EE<0 & DD<28.3) = 0; 
subplot('Position',axpos(4,:));
[DD,EE] = meshgrid(dens_levs,eta);
pcolor(EE,DD,psi_plot);
shading interp;
hold on;
[C,h] = contour(EE,DD,psi_plot,[-6 -5 -4 -3 -2 -1 -0.5 -0.25 0.25 0.5 1],'EdgeColor','k');
clabel(C,h);
hold off;
caxis([-psimax psimax]);
set(gca,'YDir','reverse');
set(gca,'YLim',ylim);
set(gca,'XLim',[-7 10.5]);
colormap(gca,cmocean('balance',round(2*psimax/(psistep))));%,'pivot',0));
xlabel('MOC coordinate, \eta');
ylabel('Neutra density, \gamma (kg/m^3)')
set(gca,'FontSize',fontsize);
cbhandle = colorbar;
set(cbhandle,'Position',cb3_pos);
title(cbhandle,'Sv');
text(-6.5,27.4,'Overturning streamfunction (mean component)','FontSize',fontsize+2);

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
% plot(tt_onethird_notides/t1year+2007-18,0*tt_onethird_notides,'k-','LineWidth',0.5);
% hold on;
% plot(tt_onethird_notides/t1year+2007-18,0*tt_onethird_notides,'k--','LineWidth',0.5);
% hold off;
set(ax2,'Color','None');
set(ax2,'XAxisLocation','Top');
set(ax2,'YLim',get(ax1,'YLim'));
set(ax2,'YTick',[]);
box off;
set(ax2,'XLim',get(ax1,'XLim')-78.16);
set(ax2,'FontSize',fontsize);
set(get(ax2,'XLabel'),'String','Pseudo-Latitude');


%%% Panel labels
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.07 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.07 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.07 axpos(4,2)-0.04 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

