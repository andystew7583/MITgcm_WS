%%%
%%% plotFRISMeltRate_WCbatch.m
%%%
%%% Plots FRIS melt rate and continental shelf buoyancy loss for all experiments in our Weddell
%%% Catastrophe batch.
%%%

%%% Pointer to experiment directory
expdir = '../experiments';

%%% Pointer to storage directory for output .mat files
proddir = './products_WCbatch';

%%% Set true to overwrite existing files
overwrite = false;


%%% Read Google Drive spreadsheet with experiment data
url = "https://docs.google.com/spreadsheets/d/19mLumSIpCtc6AsPfuzLBkuPxVXoyXuaD3yWbCmxe8FY/export?format=csv&gid=0";
T = readtable(url);

%%% Narrow down to production experiments
T = T(T.production_=="Y",:);

%%% Narrow down to experiments that are complete and downloaded
T = T(T.completed_=="Y" & T.downloaded_=="Y",:);

%%% Restrict to a specific batch of experiments
T = T(T.batch >= 1 & T.batch<=10,:);

%%% Pre-allocate storage
strat = zeros(1,size(T,1));
batchnum = zeros(1,size(T,1));
dpyc = zeros(1,size(T,1));
FRISmelt_batch = zeros(1,size(T,1));
initFRISmelt_batch = zeros(1,size(T,1));
SIprod_batch = zeros(1,size(T,1));
has_WC = zeros(1,size(T,1));

%%% Loop through experiments and compute MOC
Nexps = size(T,1);
for n = 1:Nexps
  expname = T.Name(n);
  expname = expname{1};
  tmax = T.EndTime_yr_(n) + 0.05;
  tmin = tmax - 9; %%% Last 9 years

  %%% Record upstream stratification and pycnocline depth
  strat(n) = T.strat(n);
  dpyc(n) = 400+T.Dpyc(n);
  batchnum(n) = T.batch(n);    
  if (T.WeddellCatastrophe_(n)=="Y")
    has_WC(n) = 1;
  elseif (T.WeddellCatastrophe_(n)=="N")
    has_WC(n) = 0;
  else
    has_WC(n) = -1;
  end

  %%% Load pre-computed melt rates
  datafname = [expname,'_FRISMeltRate.mat'];
  load(fullfile(proddir,datafname));
  tt = tt / t1year;

  %%% Compute time-mean FRIS melt rate in Gt/yr
  FRISmelt_batch(n) = -mean(SHImelt((tt>tmin) & (tt<=tmax)))*t1year/1e12;
  if (has_WC(n) == 0)
    initFRISmelt_batch(n) = FRISmelt_batch(n);
  else
    initFRISmelt_batch(n) = -mean(SHImelt((tt>0.05) & (tt<=3.05)))*t1year/1e12;
  end    
  SIprod_batch(n) = mean(SIprod((tt>tmin) & (tt<=tmax)))*t1year/1e12;

end

beta = 8e-4;
g = 9.81;
f0 = 2*(2*pi*366/365/86400)*sind(65);
Sref = 34.67;
rhofresh = 1000;
Sstrat = strat./(g*beta);
FWupstream = (1/26)*(g.*beta.*dpyc.^4.*Sstrat.^2)/(12*f0*Sref) * rhofresh * t1year / 1e12;
% FWnet = FWupstream + FRISmelt_batch(1) - SIprod_batch;
FWnet = FWupstream + initFRISmelt_batch(1) - SIprod_batch;

figure(1);
clf;
plot(FWnet(has_WC==1),FRISmelt_batch(has_WC==1),'ro');
hold on;
plot(FWnet(has_WC==0),FRISmelt_batch(has_WC==0),'bo');
plot(FWnet(has_WC==-1),FRISmelt_batch(has_WC==-1),'go');
hold off;
xlabel('$F_{\mathrm{upstream}} + F_{\mathrm{melt}}^{\mathrm{ref}} - F_{\mathrm{polynya}}$ (Gt/yr)','interpreter','latex','FontSize',16);
ylabel('$F_{\mathrm{melt}}$ (Gt/yr)','interpreter','latex','FontSize',16);
set(gca,'FontSize',16);
legend('Weddell Catastrophe','No Weddell Catastrophe','Intermediate State','Location','NorthWest');

legstr = {'250m','300m','350m','400m','450m','500m','550m','600m','650m','700m'};
figure(2);
clf;
% legstr = {};
for n=[1:10]
  % legstr = {legstr{:},['Batch ',num2str(n)]};
  plot(FWnet(batchnum==n),FRISmelt_batch(batchnum==n),'o-');
  if (n == 1)
    hold on;
  end
end
plot([0 0],[0 800],'k--');
set(gca,'YLim',[0 800]);
hold off;
xlabel('$F_{\mathrm{upstream}} + F_{\mathrm{melt}}^{\mathrm{ref}} - F_{\mathrm{polynya}}$ (Gt/yr)','interpreter','latex','FontSize',16);
ylabel('$F_{\mathrm{melt}}$ (Gt/yr)','interpreter','latex','FontSize',16);
set(gca,'FontSize',16);
leghandle = legend(legstr,'Location','NorthWest');
set(leghandle,'Position',[ 0.1482    0.3643    0.1384    0.4369]);


batchnum_plot = 0;

figure(3);
clf;
plot(strat(batchnum==batchnum_plot),FRISmelt_batch(batchnum==batchnum_plot),'o-');
ylabel('$F_{\mathrm{melt}}$ (Gt/yr)','interpreter','latex','FontSize',16);
xlabel('Upstream stratification (s$^{-2}$)','interpreter','latex','FontSize',16);
set(gca,'FontSize',16);

figure(4);
clf;
plot(SIprod_batch(batchnum==batchnum_plot),FRISmelt_batch(batchnum==batchnum_plot),'o-');
ylabel('$F_{\mathrm{melt}}$ (Gt/yr)','interpreter','latex','FontSize',16);
xlabel('$F_{\mathrm{polynya}}$ (Gt/yr)','interpreter','latex','FontSize',16);
set(gca,'FontSize',16);

figure(5);
clf;
plot(FWnet(batchnum==batchnum_plot),FRISmelt_batch(batchnum==batchnum_plot),'o-');
xlabel('$F_{\mathrm{upstream}} + F_{\mathrm{melt}}^{\mathrm{ref}} - F_{\mathrm{polynya}}$ (Gt/yr)','interpreter','latex','FontSize',16);
ylabel('$F_{\mathrm{melt}}$ (Gt/yr)','interpreter','latex','FontSize',16);
set(gca,'FontSize',16);



figure(6)
for n = 1:Nexps
  expname = T.Name(n);
  expname = expname{1};
  tmax = T.EndTime_yr_(n) + 0.05;
  tmin = tmax - 3; %%% Last 3 years

  %%% Load pre-computed melt rates
  datafname = [expname,'_FRISMeltRate.mat'];
  load(fullfile(proddir,datafname));
  tt = tt / t1year;

  plot (tt,smooth(-SHImelt*t1year/1e12,24));
  if (n == 1)
    hold on;
  end

end
hold off;
xlabel('Time (years)')
ylabel('$F_{\mathrm{melt}}$ (Gt/yr)','interpreter','latex','FontSize',16);
set(gca,'FontSize',16);
