%%%
%%% paper2_plotResolutionSensitivity.m
%%%
%%% Plots sensitivity of various overturning metrics to horizontal
%%% resolution.
%%% 

%%% Options
eta_cavity = 0;
eta_shelf = 3.5;
eta_full = 8;
rhoi = 934; %%% Density of ice
FRIS_lonmax = -29.9; %%% Lat/lon range defining the FRIS cavity
FRIS_latmax = -74.5;

%%% List of experiments
expdir = '../experiments';
expname_sens = {'hires_seq_onethird_RTOPO2', ...
            'hires_seq_onethird_notides_RTOPO2', ...
            'hires_seq_onesixth_RTOPO2', ...
            'hires_seq_onesixth_notides_RTOPO2', ...
            'hires_seq_onetwelfth_RTOPO2', ...
            'hires_seq_onetwelfth_notides_RTOPO2'};
Nsens = length(expname_sens);          
          
%%% Construct MOC output file name
calc_psi_eddy = true;
deform_cavity = false;
use_layers = true;
densvar = 'PD0';
outfname_MOC = ['','_MOC_',densvar];
if (calc_psi_eddy)
  if (use_layers)
    estr = '_layers';
  else
    estr = '_TRM';
  end
else
  estr = '_noeddy';
end
outfname_MOC = [outfname_MOC,estr];
if (deform_cavity)
  outfname_MOC = [outfname_MOC,'_deform'];
end
outfname_MOC = [outfname_MOC,'.mat'];

%%% Construct surface flux output file name
outfname_fluxes = ['','_surfFluxes'];
outfname_fluxes = [outfname_fluxes,'.mat'];

%%% Construct T/S output file name
outfname_TS = ['','_TSfluxes'];
outfname_TS = [outfname_TS,'.mat'];

%%% To store metrics
psi_sens = cell(1,Nsens);

psi_cavity_sens = cell(1,Nsens);
psi_shelf_sens = cell(1,Nsens);
psi_full_sens = cell(1,Nsens);
psi_eddy_sens = cell(1,Nsens);
AABWdens_full_sens = cell(1,Nsens);
AABWdens_shelf_sens = cell(1,Nsens);
FRISmelt_sens = cell(1,Nsens);
shelfBuoyLoss_sens = cell(1,Nsens);

psi_full_mean_sens = zeros(1,Nsens);
psi_shelf_mean_sens = zeros(1,Nsens);
psi_cavity_mean_sens = zeros(1,Nsens);
psi_eddy_mean_sens = zeros(1,Nsens);
AABWdens_full_mean_sens = zeros(1,Nsens);
AABWdens_shelf_mean_sens = zeros(1,Nsens);
FRISmelt_mean_sens = zeros(1,Nsens);
shelfBuoyLoss_mean_sens = zeros(1,Nsens);

psi_full_std_sens = zeros(1,Nsens);
psi_shelf_std_sens = zeros(1,Nsens);
psi_cavity_std_sens = zeros(1,Nsens);
psi_eddy_std_sens = zeros(1,Nsens);
AABWdens_full_std_sens = zeros(1,Nsens);
AABWdens_shelf_std_sens = zeros(1,Nsens);
FRISmelt_std_sens = zeros(1,Nsens);
shelfBuoyLoss_std_sens = zeros(1,Nsens);

psimean_full_sens = zeros(1,Nsens);
psimean_shelf_sens = zeros(1,Nsens);
psimean_cavity_sens = zeros(1,Nsens);
psimean_eddy_sens = zeros(1,Nsens);
AABWdensMean_full_sens = zeros(1,Nsens);
AABWdensMean_shelf_sens = zeros(1,Nsens);

%%% Loop through experiments
for m=1:Nsens
  
  %%% Need grids to compute area-integrated fluxes
  expname = expname_sens{m};
  loadexp;
  
  %%% Load pre-computed overturning/flux data
  outfname_tmp = [expname_sens{m},outfname_MOC]
  load(fullfile('products',outfname_tmp));
  Nt = length(times);
  
  %%% Load surface fluxes
  outfname_tmp = [expname_sens{m},outfname_fluxes];
  load(fullfile('products',outfname_tmp));
  
  %%% Load mean T and S
  outfname_tmp = [expname_sens{m},outfname_TS];
  load(fullfile('products',outfname_tmp),'theta_tavg');
  load(fullfile('products',outfname_tmp),'salt_tavg');
  
  %%% Truncate time series to last 8 years
  if (Nt > 96)
    Nt = 96;
    psi_mean = psi_mean(:,:,end-Nt+1:end);
    psi_eddy = psi_eddy(:,:,end-Nt+1:end);
  end
  psi_sens{m} = mean(psi_mean+psi_eddy,3)/1e6;
    
  %%% Mean "sea surface" salinity and temperature for flux calculation
  sst = NaN*ones(Nx,Ny);
  sss = NaN*ones(Nx,Ny);
  ssp = NaN*ones(Nx,Ny);
  for i=1:Nx
    for j=1:Ny
      kmin = find(hFacC(i,j,:)>0,1,'first');
      if (~isempty(kmin))
        sss(i,j) = salt_tavg(i,j,kmin);
        sst(i,j) = theta_tavg(i,j,kmin);
        ssp = -gravity*rhoConst*RC(kmin)/1e4;
      end
    end
  end
  [alpha,beta] = calcAlphaBeta(sss,sst,ssp);
  
  %%% Indices over which to integrate, i.e. defining the FRIS
  xidx = find(XC(:,1)<-29.9);
  yidx = find(YC(1,:)<-74.5);
  [DD,EE] = meshgrid(dens_levs,eta);
  
  %%% Time series of key metrics
  psimax_cavity = zeros(1,Nt);
  psimax_shelf = zeros(1,Nt);
  psimax_full = zeros(1,Nt);
  psimax_eddy = zeros(1,Nt);
  AABWdens_full = zeros(1,Nt);
  AABWdens_shelf = zeros(1,Nt);
  FRISmelt = zeros(1,Nt);
  shelfBuoyLoss = zeros(1,Nt);
  
  %%% Loop through time to compute metrics
  for n=1:Nt
    
    %%% Streamfunction strengths
    tmp = (psi_mean(:,:,n)+psi_eddy(:,:,n))/1e6;
    psimax_full(n) = - min(min(tmp(EE<eta_full)));
    psimax_shelf(n) = -min(min(tmp(EE<eta_shelf)));
    psimax_cavity(n) = max(max(tmp(EE<eta_cavity)));
    
    %%% Transport-weighted density on shelf
    jidx = find(eta>eta_shelf,1,'first');
    kidx = find(tmp(jidx,:)==min(tmp(jidx,:)),1,'last');
    AABWdens_shelf(n) = sum(diff(tmp(jidx,kidx:end),1,2).*(0.5*(dens_levs(kidx:end-1)+dens_levs(kidx+1:end)))) / (tmp(jidx,end)-tmp(jidx,kidx));
    
    %%% Transport-weighted density in open ocean
    jidx = find(eta>eta_full,1,'first');
    kidx = find(tmp(jidx,:)==min(tmp(jidx,:)),1,'last');
    AABWdens_full(n) = sum(diff(tmp(jidx,kidx:end),1,2).*(0.5*(dens_levs(kidx:end-1)+dens_levs(kidx+1:end)))) / (tmp(jidx,end)-tmp(jidx,kidx));
    
    %%% Eddy streamfunction strength
    tmp = psi_eddy(:,:,n)/1e6;
    psimax_eddy(n) = -min(min(tmp((EE<eta_shelf+2) & (EE>eta_shelf-2))));    
    
    %%% Composite surface heat flux
    htFlxDn = tflux(:,:,n);
    tmp = SHIhtFlx(:,:,n);
    htFlxDn(hFacC(:,:,1)==0) = -tmp(hFacC(:,:,1)==0);
    htFlxDn(sum(hFacC,3)==0) = NaN;

    %%% Composite virtual surface salt flux
    sltFlxDn = sflux(:,:,n);
    tmp = SHIfwFlx(:,:,n);
    sltFlxDn(hFacC(:,:,1)==0) = tmp(hFacC(:,:,1)==0).*sss(hFacC(:,:,1)==0);
    sltFlxDn(sum(hFacC,3)==0) = NaN;

    %%% Calculate buoyancy flux
    bfluxT = gravity*(alpha.*htFlxDn/Cp/rhoConst);
    bfluxS = - gravity*(beta.*sltFlxDn/rhoConst);
    bflux = bfluxT + bfluxS;

    %%% Area-integrated buoyancy flux and FRIS melt
    shelfBuoyLoss(n) = -nansum(nansum(bflux.*RAC.*(ETA<eta_shelf)));   
    FRISmelt(n) = -sum(sum(SHIfwFlx(xidx,yidx,n).*RAC(xidx,yidx)));
    
  end
  
  psi_full_sens{m} = psimax_full;
  psi_shelf_sens{m} = psimax_shelf;
  psi_cavity_sens{m} = psimax_cavity;
  psi_eddy_sens{m} = psimax_eddy;
  AABWdens_full_sens{m} = AABWdens_full;
  AABWdens_shelf_sens{m} = AABWdens_shelf;
  shelfBuoyLoss_sens{m} = shelfBuoyLoss;
  FRISmelt_sens{m} = FRISmelt;
  
  psi_full_mean_sens(m) = mean(psimax_full);
  psi_shelf_mean_sens(m) = mean(psimax_shelf);
  psi_cavity_mean_sens(m) = mean(psimax_cavity);
  psi_eddy_mean_sens(m) = mean(psimax_eddy);
  AABWdens_full_mean_sens(m) = mean(AABWdens_full);
  AABWdens_shelf_mean_sens(m) = mean(AABWdens_shelf);
  shelfBuoyLoss_mean_sens(m) = mean(shelfBuoyLoss);
  FRISmelt_mean_sens(m) = mean(FRISmelt);
  
  psi_full_std_sens(m) = std(psimax_full);
  psi_shelf_std_sens(m) = std(psimax_shelf);
  psi_cavity_std_sens(m) = std(psimax_cavity);
  psi_eddy_std_sens(m) = std(psimax_eddy);
  AABWdens_full_std_sens(m) = std(AABWdens_full);
  AABWdens_shelf_std_sens(m) = std(AABWdens_shelf);
  shelfBuoyLoss_std_sens(m) = std(shelfBuoyLoss);
  FRISmelt_std_sens(m) = std(FRISmelt);
  
  %%% Metrics using long-term mean streamfunction
  tmp = mean(psi_mean+psi_eddy,3)/1e6;
  psimean_full_sens(m) = - min(min(tmp(EE<eta_full)));
  psimean_shelf_sens(m) = -min(min(tmp(EE<eta_shelf)));
  psimean_cavity_sens(m) = max(max(tmp(EE<eta_cavity)));

  %%% Transport-weighted density on shelf
  jidx = find(eta>eta_shelf,1,'first');
  kidx = find(tmp(jidx,:)==min(tmp(jidx,:)),1,'last');
  AABWdensMean_shelf_sens(m) = sum(diff(tmp(jidx,kidx:end),1,2).*(0.5*(dens_levs(kidx:end-1)+dens_levs(kidx+1:end)))) / (tmp(jidx,end)-tmp(jidx,kidx));

  %%% Transport-weighted density in open ocean
  jidx = find(eta>eta_full,1,'first');
  kidx = find(tmp(jidx,:)==min(tmp(jidx,:)),1,'last');
  AABWdensMean_full_sens(m) = sum(diff(tmp(jidx,kidx:end),1,2).*(0.5*(dens_levs(kidx:end-1)+dens_levs(kidx+1:end)))) / (tmp(jidx,end)-tmp(jidx,kidx));

  %%% Eddy streamfunction strength
  tmp = mean(psi_eddy,3)/1e6;
  psimean_eddy_sens(m) = -min(min(tmp((EE<eta_shelf+2) & (EE>eta_shelf-2))));    

end





%%%%%%%%%%%%%%%%%%
%%% MAKE PLOTS %%%
%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
bathycntrs = [0 250 500 1000 2000 3000 4000];
axpos = zeros(6,4);
axheight = 0.2;
axwidth = 0.4;
axoff1 = 0.08;
axoff2 = 0.58;
axpos(1,:) = [axoff1 0.79 axwidth axheight];
axpos(2,:) = [axoff2 0.79 axwidth axheight];
axpos(3,:) = [axoff1 0.545 axwidth axheight];
axpos(4,:) = [axoff2 0.545 axwidth axheight];
axpos(5,:) = [axoff1 0.3 axwidth axheight];
axpos(6,:) = [axoff2 0.3 axwidth axheight];
axpos(7,:) = [axoff1 0.055 axwidth axheight];
axpos(8,:) = [axoff2 0.055 axwidth axheight];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
colororder = get(gca,'ColorOrder');
markersize = 10;
linewidth = 1.5;
errbar_width_fac = 30/29;
errbar_linewidth = 0.5;

%%% Grid of resolutions
res_sens = 1./[1/3 1/6 1/12];
res_labels = {'1/3','1/6','1/12'};
res_idx_tides = [1,3,5];
res_idx_notides = [2,4,6];




%%% Set up the figure
figure(213)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[382          55        820         900]);







%%% Overall MOC strength
subplot('Position',axpos(1,:));
semilogx(res_sens,psi_full_mean_sens(res_idx_tides),'o-','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
hold on;
semilogx(res_sens,psi_full_mean_sens(res_idx_notides),'o-','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
semilogx(res_sens,psimean_full_sens(res_idx_tides),'s--','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
semilogx(res_sens,psimean_full_sens(res_idx_notides),'s--','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = psi_full_mean_sens(res_idx_tides(m));
  std_tmp = psi_full_std_sens(res_idx_tides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
end
for m=1:length(res_idx_tides)
  mean_tmp = psi_full_mean_sens(res_idx_notides(m));
  std_tmp = psi_full_std_sens(res_idx_notides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
end
hold off;
set(gca,'YLim',[3 9]);
set(gca,'FontSize',fontsize);
% xlabel('Horizontal grid spacing ($^\circ$)','interpreter','latex');
ylabel('MOC strength (Sv)','interpreter','latex');
set(gca,'XTick',res_sens);
set(gca,'XTickLabel',res_labels);
set(gca,'XLim',[(7/8)*res_sens(1) (8/7)*res_sens(end)]);
legend('Tides','No tides','Location','NorthEast','Orientation','horizontal');

%%% Shelf MOC strength
subplot('Position',axpos(2,:));
semilogx(res_sens,psi_shelf_mean_sens(res_idx_tides),'o-','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
hold on
semilogx(res_sens,psimean_shelf_sens(res_idx_tides),'s--','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
semilogx(res_sens,psimean_shelf_sens(res_idx_notides),'s--','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = psi_shelf_mean_sens(res_idx_tides(m));
  std_tmp = psi_shelf_std_sens(res_idx_tides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
end
semilogx(res_sens,psi_shelf_mean_sens(res_idx_notides),'o-','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = psi_shelf_mean_sens(res_idx_notides(m));
  std_tmp = psi_shelf_std_sens(res_idx_notides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
end
hold off;
set(gca,'FontSize',fontsize);
% xlabel('Horizontal grid spacing ($^\circ$)','interpreter','latex');
ylabel('Shelf MOC strength (Sv)','interpreter','latex');
set(gca,'XTick',res_sens);
set(gca,'XTickLabel',res_labels);
set(gca,'XLim',[(7/8)*res_sens(1) (8/7)*res_sens(end)]);

%%% Cavity MOC strength
subplot('Position',axpos(3,:));
semilogx(res_sens,psi_cavity_mean_sens(res_idx_tides),'o-','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
hold on;
semilogx(res_sens,psimean_cavity_sens(res_idx_tides),'s--','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
semilogx(res_sens,psimean_cavity_sens(res_idx_notides),'s--','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = psi_cavity_mean_sens(res_idx_tides(m));
  std_tmp = psi_cavity_std_sens(res_idx_tides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
end
semilogx(res_sens,psi_cavity_mean_sens(res_idx_notides),'o-','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = psi_cavity_mean_sens(res_idx_notides(m));
  std_tmp = psi_cavity_std_sens(res_idx_notides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
end
hold off;
set(gca,'FontSize',fontsize);
% xlabel('Horizontal grid spacing ($^\circ$)','interpreter','latex');
ylabel('Cavity MOC strength (Sv)','interpreter','latex');
set(gca,'XTick',res_sens);
set(gca,'XTickLabel',res_labels);
set(gca,'XLim',[(7/8)*res_sens(1) (8/7)*res_sens(end)]);

%%% Eddy MOC strength
subplot('Position',axpos(4,:));
semilogx(res_sens,psi_eddy_mean_sens(res_idx_tides),'o-','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
hold on;
semilogx(res_sens,psimean_eddy_sens(res_idx_tides),'s--','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
semilogx(res_sens,psimean_eddy_sens(res_idx_notides),'s--','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = psi_eddy_mean_sens(res_idx_tides(m));
  std_tmp = psi_eddy_std_sens(res_idx_tides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
end
semilogx(res_sens,psi_eddy_mean_sens(res_idx_notides),'o-','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = psi_eddy_mean_sens(res_idx_notides(m));
  std_tmp = psi_eddy_std_sens(res_idx_notides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
end
hold off;
% set(gca,'YLim',[0 1.5]);
set(gca,'FontSize',fontsize);
% xlabel('Horizontal grid spacing ($^\circ$)','interpreter','latex');
ylabel('Eddy MOC strength (Sv)','interpreter','latex');
set(gca,'XTick',res_sens);
set(gca,'XTickLabel',res_labels);
set(gca,'XLim',[(7/8)*res_sens(1) (8/7)*res_sens(end)]);

%%% AABW density
subplot('Position',axpos(5,:));
semilogx(res_sens,AABWdens_full_mean_sens(res_idx_tides),'o-','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
hold on;
semilogx(res_sens,AABWdensMean_full_sens(res_idx_tides),'s--','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
semilogx(res_sens,AABWdensMean_full_sens(res_idx_notides),'s--','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = AABWdens_full_mean_sens(res_idx_tides(m));
  std_tmp = AABWdens_full_std_sens(res_idx_tides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
end
semilogx(res_sens,AABWdens_full_mean_sens(res_idx_notides),'o-','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = AABWdens_full_mean_sens(res_idx_notides(m));
  std_tmp = AABWdens_full_std_sens(res_idx_notides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
end
hold off;
set(gca,'FontSize',fontsize);
% xlabel('Horizontal grid spacing ($^\circ$)','interpreter','latex');
ylabel('WSBW density (kg/m$^3$)','interpreter','latex');
set(gca,'XTick',res_sens);
set(gca,'XTickLabel',res_labels);
set(gca,'XLim',[(7/8)*res_sens(1) (8/7)*res_sens(end)]);


%%% Shelf AABW density
subplot('Position',axpos(6,:));
semilogx(res_sens,AABWdens_shelf_mean_sens(res_idx_tides),'o-','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
hold on;
semilogx(res_sens,AABWdensMean_shelf_sens(res_idx_tides),'s--','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
semilogx(res_sens,AABWdensMean_shelf_sens(res_idx_notides),'s--','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = AABWdens_shelf_mean_sens(res_idx_tides(m));
  std_tmp = AABWdens_shelf_std_sens(res_idx_tides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
end
semilogx(res_sens,AABWdens_shelf_mean_sens(res_idx_notides),'o-','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = AABWdens_shelf_mean_sens(res_idx_notides(m));
  std_tmp = AABWdens_shelf_std_sens(res_idx_notides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
end
hold off;
set(gca,'YLim',[27.87 28.01]);
set(gca,'FontSize',fontsize);
% xlabel('Horizontal grid spacing ($^\circ$)','interpreter','latex');
ylabel('DSW density (kg/m$^3$)','interpreter','latex');
set(gca,'XTick',res_sens);
set(gca,'XTickLabel',res_labels);
set(gca,'XLim',[(7/8)*res_sens(1) (8/7)*res_sens(end)]);

%%% FRIS melt
subplot('Position',axpos(7,:));
semilogx(res_sens,FRISmelt_mean_sens(res_idx_tides)/1e12*t1year,'o-','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
hold on;
for m=1:length(res_idx_tides)
  mean_tmp = FRISmelt_mean_sens(res_idx_tides(m))/1e12*t1year;
  std_tmp = FRISmelt_std_sens(res_idx_tides(m))/1e12*t1year;
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
end
semilogx(res_sens,FRISmelt_mean_sens(res_idx_notides)/1e12*t1year,'o-','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = FRISmelt_mean_sens(res_idx_notides(m))/1e12*t1year;
  std_tmp = FRISmelt_std_sens(res_idx_notides(m))/1e12*t1year;
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
end
hold off;
set(gca,'YLim',[100 220]);
set(gca,'FontSize',fontsize);
xlabel('Horizontal grid spacing ($^\circ$)','interpreter','latex');
ylabel('FRIS melt (Gt/yr)','interpreter','latex');
set(gca,'XTick',res_sens);
set(gca,'XTickLabel',res_labels);
set(gca,'XLim',[(7/8)*res_sens(1) (8/7)*res_sens(end)]);

%%% Shelf buoyancy loss
subplot('Position',axpos(8,:));
semilogx(res_sens,shelfBuoyLoss_mean_sens(res_idx_tides),'o-','Color',colororder(1,:),'MarkerSize',markersize,'LineWidth',linewidth);
hold on;
for m=1:length(res_idx_tides)
  mean_tmp = shelfBuoyLoss_mean_sens(res_idx_tides(m));
  std_tmp = shelfBuoyLoss_std_sens(res_idx_tides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(1,:),'LineWidth',errbar_linewidth);
end
semilogx(res_sens,shelfBuoyLoss_mean_sens(res_idx_notides),'o-','Color',colororder(2,:),'MarkerSize',markersize,'LineWidth',linewidth);
for m=1:length(res_idx_tides)
  mean_tmp = shelfBuoyLoss_mean_sens(res_idx_notides(m));
  std_tmp = shelfBuoyLoss_std_sens(res_idx_notides(m));
  semilogx([res_sens(m) res_sens(m)],[mean_tmp+std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp+std_tmp mean_tmp+std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
  semilogx([res_sens(m)/errbar_width_fac res_sens(m)*errbar_width_fac],[mean_tmp-std_tmp mean_tmp-std_tmp],'Color',colororder(2,:),'LineWidth',errbar_linewidth);
end
hold off;
set(gca,'YLim',[-1000 16000]);
set(gca,'FontSize',fontsize);
xlabel('Horizontal grid spacing ($^\circ$)','interpreter','latex');
ylabel('Shelf buoyancy loss (m$^4$/s$^3$)','interpreter','latex');
set(gca,'XTick',res_sens);
set(gca,'XTickLabel',res_labels);
set(gca,'XLim',[(7/8)*res_sens(1) (8/7)*res_sens(end)]);

%%% Add panel labels
annotation('textbox',[axpos(1,1)-0.04 axpos(1,2)-0.04 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.04 axpos(2,2)-0.04 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.04 axpos(3,2)-0.04 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(4,1)-0.04 axpos(4,2)-0.04 0.03 0.03],'String',axlabels{4},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(5,1)-0.04 axpos(5,2)-0.04 0.03 0.03],'String',axlabels{5},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(6,1)-0.04 axpos(6,2)-0.04 0.03 0.03],'String',axlabels{6},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(7,1)-0.04 axpos(7,2)-0.04 0.03 0.03],'String',axlabels{7},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(8,1)-0.04 axpos(8,2)-0.04 0.03 0.03],'String',axlabels{8},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

  
  
  
  
  
  
  
  

