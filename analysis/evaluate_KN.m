%%%
%%% evaluate_KN.m
%%%
%%% Compares model output against Kapp Norvegia data.
%%%

%%% Load Kapp Norvegia data
load ../data/KN_section/KN_section.mat;

%%% Load experiment configuration
loadexp;

%%% Load experiment data
expdir = '../experiments';
% expname = 'WC_seq_onethird_notides_RTOPO2';
% tmin = 4.05;
% tmax = 8.05;
% expname = 'WC_seq_onethird_RTOPO2_dpyc200_alphaV0.5';
% tmin = 4.05;
% tmax = 6.05;
expname = 'WC_seq_onethird_RTOPO2_dpyc200';
tmin = 1.05;
tmax = 2.05;

%%% Simulation output iteration numbers
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Index of 17W
xidx = find(XC(:,1)>-17,1);

%%% To store monthly climatologies of T and S
theta = zeros(Ny,Nr,12);
salt = zeros(Ny,Nr,12);
navg = zeros(1,12);

%%% Time-averaging loop
for n=1:length(dumpIters)

  t = dumpIters(n)*deltaT/t1year;
  
  if ((t > tmin) && (t < tmax))

    %%% Read monthly-mean data
    T = rdmdsWrapper(fullfile(exppath,'results','THETA'),dumpIters(n));
    S = rdmdsWrapper(fullfile(exppath,'results','SALT'),dumpIters(n));

    %%% Ignore missing data
    if (isempty(T) || isempty (S))
      continue;
    end

    %%% Remove dry points
    T(hFacC==0) = NaN;
    S(hFacC==0) = NaN;

    %%% Restrict to KN section
    T = squeeze(T(xidx,:,:));
    S = squeeze(S(xidx,:,:));
   
    %%% Add to average for corresponding month (N.B. this assumes that all
    %%% simulations start in January)
    monidx = mod((n - 1),12) + 1;
    theta(:,:,monidx) = theta(:,:,monidx) + T;
    salt(:,:,monidx) = salt(:,:,monidx) + S;
    navg(monidx) = navg(monidx) + 1;

  end

end

%%% Divide by number of iterations to compute average
for m=1:length(navg)
  theta(:,:,m) = theta(:,:,m) / navg(m);
  salt(:,:,m) = salt(:,:,m) / navg(m);
end

%%% Plotting options
framepos = [1000         142        1018        1072];
fontsize = 13;
[ZZ,YY] = meshgrid(RC,YC(1,:));

%%% Generate model plots
for m = 1:length(navg)

  jrange = 1:Ny;
  
  figure(1);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  pcolor(YY(jrange,:),-ZZ(jrange,:)/1000,salt(jrange,:,m));
  shading interp;

  set(gca,'YDir','reverse');
  hold on;
    h = plot(yy,SHELFICEtopo(xidx,:)/1000,'k','LineWidth',3);  
    h = plot(yy,bathy(xidx,:)/1000,'k','LineWidth',3);  
  hold off;
  caxis([33.6 34.7]);
  if (m >= 10)
    xlabel('Latitude','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  set(gca,'XLim',[-74 -70]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('haline',20));
    set(handle,'Position',[.93 .05 .01 .9]);
  end

  figure(2);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  pcolor(YY(jrange,:),-ZZ(jrange,:)/1000,theta(jrange,:,m));
  shading interp;
  set(gca,'YDir','reverse');
  hold on;
    h = plot(yy,SHELFICEtopo(xidx,:)/1000,'k','LineWidth',3);  
    h = plot(yy,bathy(xidx,:)/1000,'k','LineWidth',3);  
  hold off;
  caxis([-2.2 1.2])
  if (m >= 10)
    xlabel('Latitude','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  set(gca,'XLim',[-74 -70]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('thermal',20));
    set(handle,'Position',[.93 .05 .01 .9]);
  end

  jrange = 1:find(YC(1,:)>-71,1);
  krange = 1:find(squeeze(RC)<-1000,1);

  figure(3);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  scatter(reshape(salt(jrange,krange,m),1,[]),reshape(theta(jrange,krange,m),1,[])); 
  if (m >= 10)
    xlabel('Salinity (g/kg)','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Temperature ($^\circ$C)','interpreter','latex','fontsize',fontsize);
  end
  axis([33 35 -2.2 1.2])

end







[ZZ,YY] = meshgrid(KN_section.pressure,KN_section.distance);

%%% Generate KN plots
for m = 1:length(navg)
  
  m_KN = mod(m + 5,12)+1;

  
  figure(4);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  pcolor(YY(:,:)/1000,ZZ(:,:)/1000,KN_section.salt(:,:,m_KN)');
  shading interp;
  set(gca,'YDir','reverse');
  caxis([33.6 34.7]);
  if (m >= 10)
    xlabel('Latitude','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  set(gca,'XLim',[0 100]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('haline',20));
    set(handle,'Position',[.93 .05 .01 .9]);
  end

  figure(5);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  pcolor(YY(:,:)/1000,ZZ(:,:)/1000,KN_section.temperature(:,:,m_KN)');
  shading interp;
  set(gca,'YDir','reverse');
  caxis([-2.2 1.2])
  if (m >= 10)
    xlabel('Latitude','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  set(gca,'XLim',[0 100]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('thermal',20));
    set(handle,'Position',[.93 .05 .01 .9]);
  end

  jrange = 1:26; %%% Within 100km of coast
  krange = 1:51; %%% Above 1000m depth

  figure(6);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  scatter(reshape(KN_section.salt(krange,jrange,m_KN),1,[]),reshape(KN_section.temperature(krange,jrange,m_KN),1,[])); 
  if (m >= 10)
    xlabel('Salinity (g/kg)','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Temperature ($^\circ$C)','interpreter','latex','fontsize',fontsize);
  end
  axis([33 35 -2.2 1.2])

end

