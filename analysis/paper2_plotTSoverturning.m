%%%
%%% paper2_plotTSoverturning.m
%%%
%%% Plots overturning circulation in T/S space.
%%%

%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%

%%% Load experiment
expdir = '../experiments';
expname = 'hires_seq_onetwelfth_notides_RTOPO2';
loadexp;

%%% Load streamfunction
calc_psi_eddy = true; %%% Set true to compute eddy-induced transports
outfname = [expname,'_MOC_TS'];
if (calc_psi_eddy)  
  estr = '_TRM';  
else
  estr = '_noeddy';
end
outfname = [outfname,estr];
outfname = [outfname,'.mat'];
load(fullfile('products',outfname));
NS = length(SS);
NT = length(TT);

%%% Load mean T/S and potential density
outfname = [expname,'_TSfluxes'];
outfname = [outfname,'.mat'];
load(fullfile('products',outfname),'theta_tavg');
load(fullfile('products',outfname),'salt_tavg');
pp_ref = -gravity*rhoConst*repmat(RC,[Nx Ny 1])/1e4; %%% Pressure is just Boussinesq hydrostatic reference pressure
theta_tavg(hFacC==0) = NaN; %%% Remove topography
salt_tavg(hFacC==0) = NaN;
pd_tavg = densjmd95(salt_tavg,theta_tavg,-gravity*rhoConst*repmat(RC(1),[Nx Ny Nr])/1e4); %%% Potential density

%%% Grids
[TTT,SSS]=meshgrid(TT,SS);
DDD = densjmd95(SSS,TTT,-RC(1)*gravity*rhoConst/1e4*ones(size(SSS))) - 1000;









%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE THE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
fontsize = 14;
framepos = [417    26   691   959];
labelspacing = 200;
axpos = zeros(6,4);
axpos(1,:) = [0.09 0.78 0.83 0.2];
axpos(2,:) = [0.09 0.53 0.83 0.2];
axpos(3,:) = [0.09 0.05 0.83 0.42];
cb1_pos = [0.93 0.78 0.02 0.2];
cb2_pos = [0.93 0.53 0.02 0.2];
cb3_pos = [0.93 0.05 0.02 0.42];
axlabels = {'(a)','(b)','(c)'};
% psimax = 2;
% psistep = 0.1;
psimax = 7;
psistep = 0.25;
colorcntrs = [-psimax:psistep:psimax];
linecntrs = [-8 -6 -4 -3 -2 -1 -0.25 0.25 0.5 1];
% linecntrs = [];

figure(205);                
clf;
set(gcf,'Position',framepos);            
                

%%% Temperature
subplot('Position',axpos(1,:));
idx = find(XC(:,1)<-38,1,'last');
[ZZ,YY]=meshgrid(RC,YC(1,:));
set(gcf,'Position',framepos);
contourf(YY,-ZZ,squeeze(theta_tavg(idx,:,:)),[-3:0.2:1],'EdgeColor','None');
hold on
plot(YY(:,1),-bathy(idx,:),'k-','LineWidth',2);
plot(YY(:,1),-SHELFICEtopo(idx,:),'k-','LineWidth',2);
hold off;
colormap(gca,cmocean('thermal',20));
caxis([-3 1]);
set(gca,'FontSize',fontsize);
%xlabel('Latitude');
ylabel('Depth (m)');
set(gca,'YDir','reverse');
set(gca,'Color',[.85 .85 .85]);
axis([-81 YC(1,end-spongethickness) 0 4000]);
cbhandle = colorbar;
set(cbhandle,'Position',cb1_pos);
title(cbhandle,'^oC');
text(-80.5,3700,'Potential temperature \theta at 38W','FontSize',fontsize+2,'fontweight','bold');

%%% Salinity
subplot('Position',axpos(2,:));
idx = find(XC(:,1)<-38,1,'last');
salt_plot = squeeze(salt_tavg(idx,:,:));
salt_plot(salt_plot<34.2) = 34.2;
[ZZ,YY]=meshgrid(RC,YC(1,:));
set(gcf,'Position',framepos);
contourf(YY,-ZZ,salt_plot,34.2:0.02:34.8,'EdgeColor','None');
hold on
plot(YY(:,1),-bathy(idx,:),'k-','LineWidth',2);
plot(YY(:,1),-SHELFICEtopo(idx,:),'k-','LineWidth',2);
hold off;
colormap(gca,cmocean('haline',20));
caxis([34.2 34.8]);
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
set(gca,'YDir','reverse');
set(gca,'Color',[.85 .85 .85]);
axis([-81 YC(1,end-spongethickness) 0 4000]);
cbhandle = colorbar;
set(cbhandle,'Position',cb2_pos);
title(cbhandle,'g/kg');
text(-80.5,3700,'Salinity S at 38W','FontSize',fontsize+2,'fontweight','bold');

%%% Streamfunction
subplot('Position',axpos(3,:));
psi_plot = mean(-psi_TS_mean_intT-psi_TS_eddy_intT,3)/1e6;
psi_plot = sign(psi_plot).*min(abs(psi_plot),psimax);
ss_bdry = [squeeze(salt_tavg(:,Ny,:)) ; squeeze(salt_tavg(Nx,:,:))];
pt_bdry = [squeeze(theta_tavg(:,Ny,:)) ; squeeze(theta_tavg(Nx,:,:))];
contourf(SSS,TTT,psi_plot,colorcntrs,'EdgeColor','None');
hold on;
[C,h] = contour(SSS,TTT,psi_plot,linecntrs,'EdgeColor','k');
clabel(C,h,'Color','k');
[C,h] = contour(SSS,TTT,DDD,[27:.1:29],'EdgeColor',[0.5 0.5 0.5]);
clabel(C,h,'Color',[0.5 0.5 0.5]);
plot(ss_bdry(:),pt_bdry(:),'o','MarkerSize',1,'MarkerFaceColor',[.7 .7 .7],'Color',[.7 .7 .7]);
hold off;
colormap(redblue(round(2*psimax/psistep)));
caxis([-psimax psimax]);
colorbar;
xlabel('S (g/kg)');
ylabel('\theta (^oC)');
axis([34 34.9 -2.5 1.2]);
set(gca,'FontSize',fontsize);
cbhandle = colorbar;
set(cbhandle,'Position',cb3_pos);
title(cbhandle,'Sv');
text(34.03,1,'T/S streamfunction','FontSize',fontsize+2,'fontweight','bold');
text(34.8,-1.8,'HSSW','FontSize',fontsize);
text(34.71,0.6,'WDW','FontSize',fontsize);
text(34.75,-2.3,'ISW','FontSize',fontsize);
text(34.4,-2,'WW','FontSize',fontsize);
text(34.1,-0.5,'AASW','FontSize',fontsize);

%%% Panel labels
annotation('textbox',[axpos(1,1)-0.07 axpos(1,2)-0.05 0.03 0.03],'String',axlabels{1},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(2,1)-0.07 axpos(2,2)-0.05 0.03 0.03],'String',axlabels{2},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[axpos(3,1)-0.07 axpos(3,2)-0.05 0.03 0.03],'String',axlabels{3},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


