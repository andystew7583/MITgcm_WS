%%%
%%% plotDepthIntegralTimeSeries
%%%
%%%


  
%%% Options
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
loadexp;



%%% Load pre-computed output
fname = [expname,'_DepthIntegrals.mat'];
load(fullfile('products',fname));

%%% Compute volume averages
EKE_avg = zeros(1,length(tt));
S_avg = zeros(1,length(tt));
msk = (SHELFICEtopo>=0) & (bathy >= -1000) & (YC<YC(1,end-spongethickness+1)) & (XC<XC(end-spongethickness+1,1));
volW = sum(sum(sum(RAW.*msk.*DRF.*hFacW)));
volS = sum(sum(sum(RAS.*msk.*DRF.*hFacS)));
volC = sum(sum(sum(RAC.*msk.*DRF.*hFacC)));
for n=1:length(tt)
  EKE_avg(n) = sum(sum(0.5.*usq_eddy_int(:,:,n).*RAW.*msk)) / volW  ...
             + sum(sum(0.5.*vsq_eddy_int(:,:,n).*RAS.*msk)) / volS;
  S_avg(n) = 2*sum(sum(salt_int(:,:,n).*RAC.*msk)) / volC;        %%% TODO factor of 2 accounts for spurious factor of 0.5 in calculation script
end



%%% Set up the figure
figure(202)
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[417   591   829   369]);
% fontsize = 18;
fontsize = 12;
defaultcolororder = get(gca,'ColorOrder');

%%% Make the plot
tyears = datenum('2008-01-01')+tt-t1month/86400/2;
h1 = plot(tyears,EKE_avg);
ax(1) = gca;
ax(2) = axes('Position',get(ax(1),'Position'));
h2 = plot(ax(2),tyears,S_avg,'Color',defaultcolororder(2,:));
set(ax(1),'YAxisLocation','Left');
set(ax(2),'YAxisLocation','Right');
set(ax(2),'XAxisLocation','Top');
set(ax(2),'Color','None');
set(h1,'LineWidth',1.5);
set(h2,'LineWidth',1.5);
hold on;
area(ax(2),[datenum('2008-01-01'),datenum('2009-01-01')],[34.45 34.62;34.45,34.62],'FaceColor',[.5 .5 .5],'FaceAlpha',0.5);
area(ax(2),[datenum('2011-01-01'),datenum('2012-01-01')],[34.45,34.62;34.45,34.62],'FaceColor','y','FaceAlpha',0.25);
% area(ax(1),[datenum('2008-01-01'),datenum('2009-01-01'),datenum('2009-01-01'),datenum('2008-01-01'),datenum('2008-01-01')],[0 0 3e-3 3e-3 0],'FaceColor',[.5 .5 .5],'FaceAlpha',0.5);
% area(ax(1),[datenum('2011-01-01'),datenum('2012-01-01'),datenum('2012-01-01'),datenum('2011-01-01'),datenum('2011-01-01')],[0 0 3e-3 3e-3 0],'FaceColor','y','FaceAlpha',0.25);
hold off;
set(ax(1),'XLim',[datenum('2008-01-01') datenum('2015-01-01')])
set(ax(2),'XLim',[datenum('2008-01-01') datenum('2015-01-01')])
ylim(ax(2),[34.45 34.62]);
set(ax(2),'YTick',[34.46:0.02:34.62]);
xlabel(ax(1),'Year');
datetick(ax(1),'x');
set(ax(2),'XTick',[]);
set(ax,'Box','off');
set(ax,'FontSize',fontsize);
set(ax,'Position',[0.06 0.12 0.85 0.8]);
ylabel(ax(1),'Shelf-averaged EKE (m^2/s^2)');
ylabel(ax(2),'Shelf-averaged salinity (g/kg)');
text(ax(2),datenum('2008-02-01'),34.46,'Spin-up period','FontSize',fontsize);
text(ax(2),datenum('2011-02-01'),34.46,'12-hourly output','FontSize',fontsize);
ax(1).YColor = defaultcolororder(1,:);
ax(2).YColor = defaultcolororder(2,:);