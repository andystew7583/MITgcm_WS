%%%%plotting BSF, TS, and SHImelt of reference run
setExpname
loadexp
load SavedFiles/bsf_controldata.mat
load SavedFiles/TSanalysis4panel.mat
load shimelt_control.mat
figure(1)
clf
s1=subplot(2,2,[1:2]);
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


contourfm(YG,XG,Psi_plot,40,'EdgeColor','none');
colormap (flipud(pmkmp(40)));  

h = colorbar;
set(gca,'FontSize',10);
set(h,'Position',[0.84 0.55 0.01 .4])
% title(h,'Sv','Fontsize',20);
% cbfreeze(h);
% 
hold on
contourm(YG,XG,Psi_plot,0:1:30,'EdgeColor',[.5,.5,.5],'LineWidth',.6)
hold on


[cs,C] = contourm(YC,XC,bathy,[-5000:1000:-1000 -500 -200 -100],'EdgeColor','black'); 
hh = clabelm(cs,C);
set (hh,'fontsize',10,'BackgroundColor','none','Edgecolor','none')
         
hold on        
contourm(YG,XG,topog_msk,1,'EdgeColor','black','ShowText','off','LineWidth',1)
        
hh = clabelm(cs,C);
set (hh,'fontsize',10,'BackgroundColor','none')

%     hh = clabelm(cs,C);
%     set (hh,'fontsize',7,'BackgroundColor','none')

caxis([0 30])
ax2 = gca;
outerpos = ax2.OuterPosition;
ti = ax2.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax2_width = outerpos(3) - ti(1) - ti(3);
ax2_height = outerpos(4) - ti(2) - ti(4);
ax2.Position = [left bottom-.01 ax2_width ax2_height];

xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
h = title('Barotropic Stream Function (Sv)','interpreter','latex','Fontsize',18);
set(h,'Position',[0 -.94 0])
text(.28,-1.3,'\textbf{a}','fontsize',18,'interpreter','latex');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom-.03 ax_width ax_height];



%%%%%%%%%%%
%%%%%%%%%%% next plot...
s2=subplot(2,2,3);
cmap=jet;
m=pcolor(s,t,new_v_3445c');
shading interp
m = colorbar;
colormap(s2,jet);
l = (min(min(new_v_3445c)));
caxis([22 30]);
ylabel(m,'Log (Volume) m$^3$','interpreter','latex','fontsize',18)
xlim([33.5 35]);
xticks(33.5:.2:36);
% caxis([18 22]);
set(gca,'OuterPosition',[.01,.01,.45,.45]);

title('Temperature/Salinity', 'interpreter','latex','FontSize',18);
xlabel('Salinity (psu)', 'interpreter','latex','FontSize',18);
ylabel('Temperature ($^{o}$C)', 'interpreter','latex','FontSize',18);
b = gca; legend(b,'off');



pt_max =1+.6;
pt_min = -3;
ss_max = 36;
ss_min = 33-1;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/100;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/100;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% Calculate potential density
pd = densmdjwf(SS_grid,PT_grid,0)-1000;
hold on
[C,f] = contour(SS_grid,PT_grid,pd,27.0:.1:33.0,'EdgeColor','k');
clabel(C,f)
% 
% 
ylim([-3 1.5]);
xlim([33.5 35]);
text(34.8,-3.6,'\textbf{b}','fontsize',18,'interpreter','latex');
set(gca,'fontsize',12)


hold on

subplot(2,2,4);
c3=plot(tt(1:300),-melttot_exp1(1:300),'b','linewidth',1);
hold on 
% n=plot(tt_3(1:meltlen_2),-melttot_exp2(1:meltlen_2),'g');
% w=plot(tt_2(1:meltlen_3),-melttot_exp3(1:meltlen_3),'r');

% 
% 
% % 
xlabel('t (years)','interpreter','latex');
ylabel('Integrated Shelf Ice Melt (Gt/yr)','interpreter','latex');
title('Integrated FRIS Melt', 'interpreter','latex','fontsize',16);
% % 
hold on
% % 

%%%shade up to 1 years
h1 = line([0 0],[0 2000]);
h2 = line([2 2],[0 2000]);
c = patch([0 2 2 0],[0 0 2000 2000],[.5 .5 .5]);

c.FaceAlpha=.1;
ylim([0 2000]);
xlim([0 25]);
hold on
z = 124*ones(40,1);
k3 = plot(0:39,z,'k:','linewidth',2);
hold on
% % 
legend([k3],{'Observed Melt'},'location','northeast','interpreter','latex','fontsize',13)
% 

text(22,-125,'\textbf{c}','fontsize',18,'interpreter','latex');
set(gca,'OuterPosition',[.54,.01,.45,.45]);

set(gca,'fontsize',12);
