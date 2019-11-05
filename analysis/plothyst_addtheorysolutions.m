load SavedFiles/salthyst_forplotting.mat
load SavedFiles/temp_hyst_forplotting.mat
load SavedFiles/shimelt_hyst_forplotting.mat

setExpname
loadexp 



c_s = [.5,.6,.7,.8,.9,.95,1,1.05, 1.1,1.15, 1.2];
w_s = [.5,.6,.7,.8,.9,.95,1,1.05,1.1,1.15,1.2];
c_s_1 = [.7,.8,.9,.95,1,1.05, 1.1];
w_s_1 = [.7,.8,.9,.95,1,1.05,1.1];

hssw_vals = [34.3554   34.3959   34.4364   34.4566   34.4769   34.4971   34.5174];
CDW_temp = -1*ones(1,7);
CDW_salt= 34.35*ones(1,7);
HSSW_temp = -2.1*ones(1,7);
shimelt_hi = 846*ones(1,11);
shimelt_lo = 172*ones(1,11);
 
Vwind_lo = .43*ones(1200,1);
Vwind_hi = 1.16*ones(1200,1);
dataw = [Cavity_salt_wminus30n,Cavity_salt_wminus20n,Cavity_salt_wminus10n,Cavity_salt_wminus5n,Cavity_salt_wcontroln,Cavity_salt_wplus5n,Cavity_salt_wplus10n];
datac = [Cavity_salt_cminus30n,Cavity_salt_cminus20n,Cavity_salt_cminus10n,Cavity_salt_cminus5n,Cavity_salt_controln,Cavity_salt_cplus5n,Cavity_salt_cplus10n];

datawtemp = [Cavity_temp_wminus30n,Cavity_temp_wminus20n,Cavity_temp_wminus10n,Cavity_temp_wminus5n,Cavity_temp_wcontroln,Cavity_temp_wplus5n,Cavity_temp_wplus10n];
datactemp = [Cavity_temp_cminus30n,Cavity_temp_cminus20n,Cavity_temp_cminus10n,Cavity_temp_cminus5n,Cavity_temp_controln,Cavity_temp_cplus5n,Cavity_temp_cplus10n];


datawshimelt = [SHImelt_control_wminus50,SHImelt_control_wminus40,SHImelt_control_wminus30,SHImelt_control_wminus20,SHImelt_control_wminus10,SHImelt_control_wminus5, SHImelt_control_warm,SHImelt_control_wplus5,SHImelt_control_wplus10,SHImelt_control_wplus15,SHImelt_control_wplus20];
datacshimelt = [SHImelt_control_minus50,SHImelt_control_minus40, SHImelt_control_minus30,SHImelt_control_minus20,SHImelt_control_minus10, SHImelt_control_minus5,SHImelt_control,SHImelt_control_plus5,SHImelt_control_plus10,SHImelt_control_plus15,SHImelt_control_plus20];



figure(1)
subplot(3,1,1)

%%%%plot warm 
p5 = plot(w_s_1,datawtemp,'-o','Color','r','LineWidth',2,'MarkerSize',10,'Markerfacecolor','r');
hold on
plot(w_s_1(5),datawtemp(5),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');

%%%% plot cold
p6 = plot(c_s_1,datactemp,'-o','Color','b','LineWidth',2,'MarkerSize',10,'Markerfacecolor','b');
hold on
plot(c_s_1(5),datactemp(5),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');
ylabel('\begin{tabular}{c} Melt-Weighted \\ T$_{Cavity}$ ($^\circ$C)\end{tabular}','interpreter','latex','fontsize',14);
hold on
a=plot(w_s_1,CDW_temp,'r--');
hold on
b=plot(c_s_1,HSSW_temp,'b--');

t=title('Wind Perturbation ($\chi$) vs Melt-Weighted FRIS Bottom Cavity Temperature','FontSize',15,'interpreter','latex');
legend([p5,p6,a,b],{'Initial FRESH','Initial REF','$\theta_{CDW}$ (Assumed)','$\theta_{HSSW}$ (Assumed)'},'location','northwest','interpreter','latex','fontsize',13);
xlim([.4 1.3]);


hold on
subplot(3,1,2)
%%%%plot warm 
p7 = plot(w_s_1,dataw,'-o','Color','r','LineWidth',2,'MarkerSize',10,'Markerfacecolor','r');
hold on
plot(w_s_1(5),dataw(5),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');

%%%% plot cold
p8 = plot(c_s_1,datac,'-o','Color','b','LineWidth',2,'MarkerSize',10,'Markerfacecolor','b');
hold on
plot(c_s_1(5),datac(5),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');
ylabel('\begin{tabular}{c} Melt-Weighted \\ S$_{Cavity}$(psu)\end{tabular}','interpreter','latex','fontsize',14);
hold on
ere1=plot(w_s_1,CDW_salt,'r--');
hold on
ere2=plot(c_s_1,hssw_vals,'b--');

t=title('Wind Perturbation ($\chi$) vs Melt-Weighted FRIS Bottom Cavity Salinity (psu)','FontSize',15,'interpreter','latex');
xlim([.4 1.3]);
legend([ere1,ere2],{'$S_{CDW}$ (Assumed)','$S_{HSSW}$ (Predicted)'},'location','northwest','interpreter','latex','fontsize',13);


subplot(3,1,3)

p3 = plot(w_s,datawshimelt,'-o','Color','r','LineWidth',2,'MarkerSize',10,'Markerfacecolor','r');
hold on
plot(w_s(7),datawshimelt(7),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');

%%%% plot cold
p4 = plot(c_s,datacshimelt,'-o','Color','b','LineWidth',2,'MarkerSize',10,'Markerfacecolor','b');
hold on
plot(c_s(7),datacshimelt(7),'-sk','LineWidth',3,'MarkerSize',15,'Markerfacecolor','k');
ylabel('Integrated Melt (Gt)','interpreter','latex','fontsize',14,'Color','k');

t=title('Wind Perturbation ($\chi$) vs Integrated Shelf Ice Melt (Gt)','FontSize',15,'interpreter','latex');


hold on
text(.8,1020,'$\Sigma_{melt}=\Sigma_{melt_{0}}+C_{m}(\theta_{CDW}-\theta_{f})$','interpreter','latex','fontsize',12);
hold on
text(.8,100,'$\Sigma_{melt}=\Sigma_{melt_{0}}+C_{m}(\theta_{HSSW}-\theta_{f})$','interpreter','latex','fontsize',12);
hold on
vlow = text(.45,50,'V$_{wind}$ = V$_{lower}$','interpreter','latex','fontsize',13);
vhi = text(1.22,50,'V$_{wind}$ = V$_{upper}$','interpreter','latex','fontsize',13);
set(vhi,'Rotation',-90);

set(vlow,'Rotation',90);
plot(w_s,shimelt_hi,'r--');
hold on
plot(c_s,shimelt_lo,'b--');
hold on

plot(Vwind_hi,1:1200,'k--');
hold on
plot(Vwind_lo,1:1200,'k--');
hold on

axis([.4 1.3 0 1200]);


xlabel('Wind Perturbation ($\chi$)','interpreter','latex','fontsize',18);
% legend([p3,p4],{'Initial FRESH','Initial REF'},'location','north','interpreter','latex','fontsize',13);





