%%%% plotting HSSW and melt/temp analysis
hysteresisHSSW
load cavtempanal.mat


figure(1)
clf
p1=subplot(1,2,1)
datanewwarm = [SALT_winter_tot_wminus50,SALT_winter_tot_cminus50,SALT_winter_tot_wminus40,SALT_winter_tot_cminus40, SALT_winter_tot_wminus30,SALT_winter_tot_wminus20,SALT_winter_tot_wminus10,SALT_winter_tot_wminus5,SALT_winter_tot_1,SALT_winter_tot_wplus5,SALT_winter_tot_wplus10];
w_s=[.5,.5,.6,.6,.7,.8,.9,.95,1,1.05,1.1];



datanewcold= [ SALT_winter_tot_cminus30,SALT_winter_tot_cminus20,SALT_winter_tot_cminus10,SALT_winter_tot_cminus5,SALT_winter_tot_control,SALT_winter_tot_cplus5,SALT_winter_tot_cplus10,SALT_winter_tot_cplus15,SALT_winter_tot_wplus15,SALT_winter_tot_cplus20,SALT_winter_tot_wplus20];
w_c=[.7,.8,.9,.95,1,1.05,1.1,1.15,1.15,1.2,1.2];


dataw = [SALT_winter_tot_wminus50,SALT_winter_tot_wminus40, SALT_winter_tot_wminus30,SALT_winter_tot_wminus20,SALT_winter_tot_wminus10, SALT_winter_tot_wminus5,SALT_winter_tot_1,SALT_winter_tot_wplus5,SALT_winter_tot_wplus10,SALT_winter_tot_wplus15,SALT_winter_tot_wplus20];
datac = [SALT_winter_tot_cminus50,SALT_winter_tot_cminus40, SALT_winter_tot_cminus30,SALT_winter_tot_cminus20,SALT_winter_tot_cminus10, SALT_winter_tot_cminus5,SALT_winter_tot_control,SALT_winter_tot_cplus5,SALT_winter_tot_cplus10,SALT_winter_tot_cplus15,SALT_winter_tot_cplus20];

hold on
errorbar(w_s(1),datanewwarm(1),err_exp4w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(2),datanewwarm(2),err_exp13w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(3),datanewwarm(3),err_exp12w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(4),datanewwarm(4),err_exp8w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(5),datanewwarm(5),err_exp7w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(6),datanewwarm(6),err_exp9w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(7),datanewwarm(7),err_exp1w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on
errorbar(w_s(8),datanewwarm(8),err_exp11w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on

errorbar(w_s(9),datanewwarm(9),err_exp5w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on

errorbar(w_s(10),datanewwarm(10),err_exp6w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on

errorbar(w_s(11),datanewwarm(11),err_exp2w,'-ro','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
hold on




%%%% plot cold
% c = plot(w_c,datac,'-ob','LineWidth',2);
hold on
errorbar(w_c(1),datanewcold(1),err_exp1,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(2),datanewcold(2),err_exp12,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');

hold on
errorbar(w_c(3),datanewcold(3),err_exp11,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(4),datanewcold(4),err_exp9,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(5),datanewcold(5),err_exp10,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(6),datanewcold(6),err_exp7,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(7),datanewcold(7),err_exp5,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(8),datanewcold(8),err_exp4,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(9),datanewcold(9),err_exp2,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(10),datanewcold(10),err_exp13,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
errorbar(w_c(11),datanewcold(11),err_exp3,'-bo','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on



b=title('\begin{tabular}{c} JJA Salt Flux $(g/s)$, \\ Ronne Polynya \end{tabular}','FontSize',15,'interpreter','latex');
xlabel('Wind Perturbation ($\chi$)','Fontsize',16,'interpreter','latex')
ylabel('Integrated Salt Flux $(g/s)$','Fontsize',16,'interpreter','latex')


plot(w_s,f,'--r')


hold on
plot(w_c,fw,'--b')


text(1.08,1.1e8,'R$^2$=.97','fontsize',18,'color','b','interpreter','latex');
text(1.08,1.5e8,'R$^2$=.96','fontsize',18,'color','r','interpreter','latex');
text(.45,.25e8,'a','fontsize',18,'color','k','interpreter','latex');
xlim([.4 1.3]);
ylim([0 5e8]);
hold on

set(gca,'Position',[.05 .1 .4 .8]);
set(p1,'Box','on');
set(gca,'fontsize',14);
pbaspect([1 1 1])



subplot(1,2,2)


plot(Temp_diff_control_n,melttot_control_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_control_n+.001,melttot_control_n-.005,'C')

hold on
plot(Temp_diff_cplus5_n,melttot_c_meridplus5_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cplus5_n+.001,melttot_c_meridplus5_n-.005,'+5%')

hold on
plot(Temp_diff_cplus10_n,melttot_c_meridplus10_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cplus5_n+.001,melttot_c_meridplus5_n-.005,'+10%')

hold on
plot(Temp_diff_cplus15_n,melttot_c_meridplus15_n,'-bo','markerfacecolor','blue','markersize',10);

hold on
plot(Temp_diff_cplus10_n,melttot_c_meridplus10_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cplus10_n+.001,melttot_c_meridplus5_n-.005,'+10%')

hold on
plot(Temp_diff_cplus20_n,melttot_c_meridplus20_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cplus20_n+.001,melttot_c_meridplus20_n-.005,'+20%')
hold on

plot(Temp_diff_cminus50_n,melttot_c_meridminus50_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_cplus20_n+.001,melttot_c_meridplus20_n-.005,'+20%')

hold on
% plot(Temp_diff_cminus2_n,melttot_c_meridminus2_n,'-bo','markerfacecolor','blue');
% hold on
plot(Temp_diff_cminus5_n,melttot_c_meridminus5_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cminus5_n+.001,melttot_c_meridminus5_n-.005,'-5%')

hold on
plot(Temp_diff_cminus10_n,melttot_c_meridminus10_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cminus10_n+.001,melttot_c_meridminus10_n-.005,'-10%')

hold on
plot(Temp_diff_cminus20_n,melttot_c_meridminus20_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_cminus20_n+.001,melttot_c_meridminus20_n-.005,'-5%')'markersize',10

hold on
plot(Temp_diff_cminus40_n ,melttot_c_meridminus40_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_cminus15_n+.001,melttot_c_meridminus15_n-.005,'-15%')

hold on




plot(Temp_diff_w_control_n,melttot_w_control_n,'-ro','markerfacecolor','red','markersize',10);
hold on
% text(Temp_diff_w_control_n+.001,melttot_w_control_n-.005,'Control')
% 
plot(Temp_diff_w_plus5_n,melttot_w_meridplus5_n,'-ro','markerfacecolor','red','markersize',10);
% % text(Temp_diff_w_plus5_n,melttot_w_meridplus5_n,'+5%')

hold on
plot(Temp_diff_w_plus10_n,melttot_w_meridplus10_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_w_plus10_n+.001,melttot_w_meridplus10_n+.001,'+10%')

hold on
% plot(Temp_diff_w_plus2_n,melttot_w_meridplus2_n,'-ro','markerfacecolor','red');
% text(Temp_diff_w_plus2_n+.001,melttot_w_meridplus2_n+.001,'+2%')
% hold on

plot(Temp_diff_w_plus20_n,melttot_w_meridplus20_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_w_plus20_n+.001,melttot_w_meridplus20_n+.001,'+20%')

hold on
plot(Temp_diff_w_minus5_n,melttot_w_meridminus5_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_w_minus5_n+.001,melttot_w_meridminus5_n+.001,'-5%')

hold on
plot(Temp_diff_w_minus10_n,melttot_w_meridminus10_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_w_minus10_n+.0001,melttot_w_meridminus10_n+.0001,'-10%')

hold on
plot(Temp_diff_w_minus20_n,melttot_w_meridminus20_n,'-ro','markerfacecolor','red','markersize',10);
% text(Temp_diff_w_minus20_n+.0001,melttot_w_meridminus20_n+.0001,'-20%')


hold on
plot(Temp_diff_w_plus15_n,melttot_w_meridplus15_n,'-bo','markerfacecolor','blue','markersize',10);
% text(Temp_diff_w_plus15_n+.001,melttot_w_meridplus15_n+.001,'+15%')

hold on
plot(Temp_diff_w_minus40_n,melttot_w_meridminus40_n,'-ro','markerfacecolor','red','markersize',10);
% xlim([.15 .25]);

hold on
plot(Temp_diff_w_minus50_n,melttot_w_meridminus50_n,'-ro','markerfacecolor','red','markersize',10);

xlim([0 1.5])



hold on
plot(totc,z,'--b','linewidth',3);


hold on
plot(totw,h,'--r','linewidth',3);



text(.2,-12e8,'R$^2$=.92','fontsize',19,'color','b','interpreter','latex')
text(.2,-11e8,'R$^2$=.98','fontsize',19,'color','r','interpreter','latex')
text(.15,-14.25e8,'b','fontsize',18,'color','k','interpreter','latex');

%%%line of best fit

b=title('\begin{tabular}{c} FRIS Freezing Bottom Temperature Anomaly($^\circ$C) vs. \\ Salt Flux from Freshwater $(g/s)$ \end{tabular}','FontSize',15,'interpreter','latex');
ylabel('Integrated Salt Flux from Freshwater $(g/s)$','FontSize',16,'interpreter','latex')
xlabel('Bottom Temperature Anomaly from Freezing ($^o$C)','FontSize',16,'interpreter','latex');
ylim([-1.5e9 -1e8]);
set(gca,'fontsize',14);

set(gca,'Position',[.55 .1 .4 .8]);

pbaspect([1 1 1])