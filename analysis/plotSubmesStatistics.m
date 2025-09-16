%%%
%%% plotSubmesStatistics.m
%%%
%%% Plots statistics relevant to submesoscale motions.
%%%

%%% Need scripts from analysis base directory
addpath ..;

%%% Read experiment data
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
expdir = '../experiments/';
% loadexp;

%%% Load pre-computed statistics
load(fullfile('./products',[expname,'_SubmesStats.mat']));
binsize = binedges(2) - binedges(1); %%% Assumes regular bin width

%%% Monthly statistics
days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
vorthist_monthly = zeros(length(days_per_month),size(vorthist,2));
divhist_monthly = 0*vorthist_monthly;
cntr = 0;
for n=1:length(days_per_month)
  len = 2*days_per_month(n);
  idx = cntr+(1:len);
  vorthist_monthly(n,:) = mean(vorthist(idx,:),1);
  divhist_monthly(n,:) = mean(divhist(idx,:),1);  
  cntr = cntr + len;
end


%%% Quarterly statistics
vorthist_quarterly = zeros(4,size(vorthist,2));
divhist_quarterly = 0*vorthist_quarterly;
for n=1:4
  vorthist_quarterly(n,:) = mean(vorthist_monthly((3*(n-1)+1):3*n,:),1);
  divhist_quarterly(n,:) = mean(divhist_monthly((3*(n-1)+1):3*n,:),1);
end


figure(71);
semilogy(binmid,mean(vorthist,1),'k-','LineWidth',2);
hold on;
semilogy(binmid,vorthist_monthly,'LineWidth',1);
hold off;
legend('Annual','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');

figure(72);
semilogy(binmid,mean(vorthist,1),'k-','LineWidth',2);
hold on;
semilogy(binmid,vorthist_monthly([11 12 1 2],:),'LineWidth',1);
hold off;
legend('Annual','Nov','Dec','Jan','Feb');

figure(73);
semilogy(binmid,mean(vorthist,1),'k-','LineWidth',2);
hold on;
semilogy(binmid,vorthist_monthly([3 4 5 6],:),'LineWidth',1);
hold off;
legend('Annual','Mar','Apr','May','Jun');

figure(74);
semilogy(binmid,mean(vorthist,1),'k-','LineWidth',2);
hold on;
semilogy(binmid,vorthist_monthly([7 8 9 10],:),'LineWidth',1);
hold off;
legend('Annual','Jul','Aug','Sep','Oct');

figure(75);
set(gcf,'Position',[288         299        1537         988]);
semilogy(binmid,mean(vorthist,1)/sum(mean(vorthist,1),2)/binsize,'k-','LineWidth',2);
hold on;
semilogy(binmid,vorthist_quarterly./sum(vorthist_quarterly,2)/binsize,'LineWidth',1);
hold off;
legend('Annual','JFM','AMJ','JAS','OND');
grid on;
xlabel('$\zeta/|f|$','interpreter','latex');
ylabel('$P(\zeta/|f|)$','interpreter','latex');
set(gca,'YLim',[3e-4 2e1]);
set(gca,'XLim',[-2 2]);
set(gca,'Color','w');
set(gca,'Position',[.07 .07 .9 .9]);
set(gca,'FontSize',16);
print('-dpng','-r300','Figures/submes/VorticityPDF.png');

figure(76);
set(gcf,'Position',[288         299        1537         988]);
semilogy(binmid,mean(divhist,1)/sum(mean(divhist,1),2)/binsize,'k-','LineWidth',2);
hold on;
semilogy(binmid,divhist_quarterly./sum(divhist_quarterly,2)/binsize,'LineWidth',1);
hold off;
legend('Annual','JFM','AMJ','JAS','OND');
grid on;
xlabel('$\delta/|f|$','interpreter','latex');
ylabel('$P(\delta/|f|)$','interpreter','latex');
set(gca,'YLim',[3e-4 2e1]);
set(gca,'XLim',[-2 2]);
set(gca,'Color','w');
set(gca,'Position',[.07 .07 .9 .9]);
set(gca,'FontSize',16);
print('-dpng','-r300','Figures/submes/DivergencePDF.png');