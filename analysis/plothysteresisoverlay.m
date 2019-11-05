%%%% plot hyst overlay

load SavedFiles/salthyst_forplotting.mat
load SavedFiles/temp_hyst_forplotting.mat
load SavedFiles/shimelt_hyst_forplotting.mat

c_s = [.5,.6,.7,.8,.9,.95,1,1.05, 1.1,1.15, 1.2];
w_s = [.5,.6,.7,.8,.9,.95,1,1.05,1.1,1.15,1.2];


dataw = [Cavity_salt_wminus50n,Cavity_salt_wminus40n,Cavity_salt_wminus30n,Cavity_salt_wminus20n,Cavity_salt_wminus10n,Cavity_salt_wminus5n,Cavity_salt_wcontroln,Cavity_salt_wplus5n,Cavity_salt_wplus10n,Cavity_salt_wplus15n,Cavity_salt_wplus20n];
datac = [Cavity_salt_cminus50n,Cavity_salt_cminus40n,Cavity_salt_cminus30n,Cavity_salt_cminus20n,Cavity_salt_cminus10n,Cavity_salt_cminus5n,Cavity_salt_controln,Cavity_salt_cplus5n,Cavity_salt_cplus10n,Cavity_salt_cplus15n,Cavity_salt_cplus20n];

datawtemp = [Cavity_temp_wminus50n,Cavity_temp_wminus40n,Cavity_temp_wminus30n,Cavity_temp_wminus20n,Cavity_temp_wminus10n,Cavity_temp_wminus5n,Cavity_temp_wcontroln,Cavity_temp_wplus5n,Cavity_temp_wplus10n,Cavity_temp_wplus15n,Cavity_temp_wplus20n];
datactemp = [Cavity_temp_cminus50n,Cavity_temp_cminus40n,Cavity_temp_cminus30n,Cavity_temp_cminus20n,Cavity_temp_cminus10n,Cavity_temp_cminus5n,Cavity_temp_controln,Cavity_temp_cplus5n,Cavity_temp_cplus10n,Cavity_temp_cplus15n,Cavity_temp_cplus20n];


datawshimelt = [SHImelt_control_wminus50,SHImelt_control_wminus40,SHImelt_control_wminus30,SHImelt_control_wminus20,SHImelt_control_wminus10,SHImelt_control_wminus5, SHImelt_control_warm,SHImelt_control_wplus5,SHImelt_control_wplus10,SHImelt_control_wplus15,SHImelt_control_wplus20];
datacshimelt = [SHImelt_control_minus50,SHImelt_control_minus40, SHImelt_control_minus30,SHImelt_control_minus20,SHImelt_control_minus10, SHImelt_control_minus5,SHImelt_control,SHImelt_control_plus5,SHImelt_control_plus10,SHImelt_control_plus15,SHImelt_control_plus20];



cmap=jet(15);
w1 = cmap(10,:);
w2 = cmap(5,:);
% 
% 
% p2 = plot(w_s,datawshimelt,'-ok','LineWidth',2,'MarkerSize',10,'Markerfacecolor','r');
% hold on
% p4 = plot(w_s,dataw,'-o','Color',w1,'linewidth',2,'MarkerSize',10,'Markerfacecolor','r');
% hold on
% p6 = plot(w_s,datawtemp,'-o','Color',w2,'linewidth',2,'MarkerSize',10,'Markerfacecolor','r');
% hold on


purp = [0.4940, 0.1840, 0.5560]	;
green = [73,156,40]/256;
%%%% [plotting with 3 yaxes]
figure(1)
clf
one = subplot(2,1,1)
%%%%changing line properties (colors)
b=[0, 0.4470, 0.7410];
b2=[0.3010, 0.7450, 0.9330]	;
r1 = [0.8500, 0.3250, 0.0980];
r2=[0.6350, 0.0780, 0.1840]	;

%%%%plot warm 
p5 = plot(w_s,datawtemp,'-o','Color',purp,'LineWidth',2,'MarkerSize',10,'Markerfacecolor','r');
hold on
plot(w_s(7),datawtemp(7),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');

%%%% plot cold
p6 = plot(c_s,datactemp,'-o','Color',purp,'LineWidth',2,'MarkerSize',10,'Markerfacecolor','b');
hold on
plot(c_s(7),datactemp(7),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');
ylabel('\begin{tabular}{c} Melt-Weighted \\ T$_{Cavity}$\end{tabular}','interpreter','latex','fontsize',14,'color',purp);
set(gca,'Ycolor',purp);
hold on

yyaxis right

%%%%plot warm 
p3 = plot(w_s,datawshimelt,'-o','Color',green,'LineWidth',2,'MarkerSize',10,'Markerfacecolor','r');
hold on
plot(w_s(7),datawshimelt(7),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');

%%%% plot cold
p4 = plot(c_s,datacshimelt,'-o','Color',green,'LineWidth',2,'MarkerSize',10,'Markerfacecolor','b');
hold on
plot(c_s(7),datacshimelt(7),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');
ylabel('Integrated FRIS Melt','interpreter','latex','fontsize',14,'Color',green);
set(gca,'Ycolor',green);

b=title('\begin{tabular}{c} Wind Perturbation ($\chi$) vs Melt-Weighted FRIS Bottom Cavity Temperature ($^\circ$C) \\ and Integrated Shelf Ice Melt (Gt) \end{tabular}','FontSize',15,'interpreter','latex');

% xlabel('Wind Perturbation ($\chi$)','interpreter','latex','fontsize',14);
hold on;
set(gca,'fontsize',12)
text(1.15,1100,'\textbf{a}','fontsize',18,'interpreter','latex');
% h = zeros(2, 1);
% h(1) = plot(NaN,NaN,'or','MarkerFaceColor','r');
% h(2) = plot(NaN,NaN,'ob','MarkerFaceColor','b');
% legend(h, {'Initial Warm FRIS Cavity','Initial Cold FRIS Cavity'},'location','north','interpreter','latex','fontsize',13);

subplot(2,1,2)
% hlines{3}(1).LineStyle = '--';
% hlines{3}(1).Marker = 'o';
% hlines{3}(1).MarkerFaceColor = 'b';
% hlines{3}(1).Color = green;
% hlines{3}(1).LineWidth = 2.0;
% hlines{3}(1).MarkerSize = 10;
% hlines{3}(2).LineStyle = '--';
% hlines{3}(2).Marker = 'o';
% hlines{3}(2).MarkerFaceColor = 'r';
% hlines{3}(2).Color = green;
% hlines{3}(2).LineWidth = 2.0;
% hlines{3}(2).MarkerSize = 10;
  %%%%plot warm 
p1 = plot(w_s,dataw,'-or','LineWidth',2,'MarkerSize',10,'Markerfacecolor','r');
hold on
plot(w_s(7),dataw(7),'-sk','MarkerSize',15,'Markerfacecolor','k');
hold on
%%%% plot cold
p2 = plot(c_s,datac,'-ob','LineWidth',2,'MarkerSize',10,'Markerfacecolor','b');
hold on
plot(c_s(7),datac(7),'-sk','MarkerSize',15,'Markerfacecolor','k');

%%%changing y-axis properties
% ax(1).YColor='k';
% ax(2).YColor = purp;
% % ax(3).YColor = green;
% % ax(3).FontSize=12;
% ax(2).FontSize=12;
% ax(1).FontSize=12;

%%%legend
hold on;

title('Wind Perturbation ($\chi$) vs Melt-Weighted FRIS Bottom Cavity Salinity (psu)','interpreter','latex','FontSize',20);
ylabel('\begin{tabular}{c} Melt-Weighted \\ S$_{Cavity}$ \end{tabular}','interpreter','latex','fontsize',14,'color','k');
xlabel('Wind Perturbation ($\chi$)','interpreter','latex','fontsize',14);
text(1.15,34.3,'\textbf{b}','fontsize',18,'interpreter','latex');

legend([p1,p2],{'Initial FRESH FRIS Cavity','Initial REF FRIS Cavity'},'location','north','interpreter','latex','fontsize',13);
set(gca,'fontsize',13)
