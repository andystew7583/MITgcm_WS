%%%% plotting tracer time series
load SavedFiles/tracertimeseries.mat
load SavedFiles/melt343445data.mat

figure(1);
clf;
set(gcf,'Position',[406   347   862   635]);
% axes('FontSize',16);
subplot(3,1,1)

cn=plot(tt(1:300),thetatot_exp1(1:300),'b','linewidth',1);
hold on
w=plot(tt(1:300),thetatot_exp2(1:300),'r','linewidth',1);
hold on
% % n= plot(tt(1:228),thetatot_exp3(1:228),'g','linewidth',1);
xlim([0 25])
% xlabel('t (years)','interpreter','latex');
ylabel('$\theta_{\mathrm{cavity}}$ ($^\circ$C)','interpreter','latex','fontsize',15);
title('Melt-Weighted Bottom Cavity Temperature ','interpreter','latex','fontsize',16);
h1 = line([0 0],[-3 0]);
h2 = line([2 2],[-3 0]);
c = patch([0 2 2 0],[-3 -3 0 0],[.5 .5 .5]);
c.FaceAlpha=.3;
axis([0 25 -2.5 -.5]);

z = 155*ones(40,1);
k3 = plot(0:39,z,'k:','linewidth',2);
hold on
handle1 = legend([cn,w],{'REF','FRESH'},'location','east','interpreter','latex','fontsize',14);
set(handle1,'Position',[0.7607    0.7604    0.1323    0.0776]);
text(22,-2.7,'\textbf{a}','fontsize',14,'interpreter','latex','fontweight','bold');
set(gca,'FontSize',14);


subplot(3,1,2)

cn=plot(tt(1:300),salttot_exp1(1:300),'b','linewidth',1);
hold on
w=plot(tt(1:300),salttot_exp2(1:300),'r','linewidth',1);
hold on
% % n= plot(tt(1:228),thetatot_exp3(1:228),'g','linewidth',1);
xlim([0 25])
% xlabel('t (years)','interpreter','latex');
ylabel('S$_{\mathrm{cavity}}$ (psu)','interpreter','latex','fontsize',15);
title('Melt-Weighted Bottom Cavity Salinity','interpreter','latex','fontsize',16);
hold on
h1 = line([0 0],[33 35]);
h2 = line([2 2],[33 35]);
c = patch([0 2 2 0],[33 33 35 35],[.5 .5 .5]);
c.FaceAlpha=.3;
axis([0 25 34.3 34.55]);
z = 155*ones(40,1);
k3 = plot(0:39,z,'k:','linewidth',2);
hold on
% legend([cn,w],{'34.45 psu','34 psu'},'location','southeast','interpreter','latex','fontsize',16);
text(22,34.25,'\textbf{b}','fontsize',14,'interpreter','latex','fontweight','bold');
set(gca,'FontSize',14);

subplot(3,1,3)
c1=plot(tt(1:300),-melttot_exp1(1:300),'b','linewidth',1);
hold on
w=plot(tt(1:300),-melttot_exp3(1:300),'r');

% 
% 
% % 
xlabel('t (years)','interpreter','latex','fontsize',17);
ylabel('Gt/yr','interpreter','latex','fontsize',15);
title('Integrated FRIS Melt', 'interpreter','latex','fontsize',16);


%%%shade up to 1 years
h1 = line([0 0],[1 2000]);
h2 = line([2 2],[1 2000]);
c = patch([0 2 2 0],[1 1 2000 2000],[.5 .5 .5]);


c.FaceAlpha=.3;
axis([0 25 0 1200]);

z = 124*ones(40,1);
k3 = plot(0:39,z,'k:','linewidth',2);
hold on
% % 
handle = legend([k3],{'Observed Melt'},'location','south','interpreter','latex','fontsize',13);
set(handle,'Position',[0.7242    0.1646    0.1698    0.0304]);
% 

text(22,-150,'\textbf{c}','fontsize',14,'interpreter','latex','fontweight','bold');
set(gca,'FontSize',14);



% %%%%%find R2 values
% 
% melt_warm = melttot_exp3(2:301);
% melt_cold = melttot_exp1(1:300);
% thetatot_exp1=thetatot_exp1(1:300);
% thetatot_exp2 = thetatot_exp2(2:301);
% w = fitlm(melt_warm,thetatot_exp2);
% c = fitlm(thetatot_exp1,melt_cold);
% 
% meanwarmtemp = mean(thetatot_exp2);
% meancoldtemp = mean(thetatot_exp1);
% 
% meancoldsalt = nanmean(salttot_exp1);
% meanwarmsalt = nanmean(salttot_exp2);
