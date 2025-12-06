%%%
%%% evaluate_KN.m
%%%
%%% Compares model output against Kapp Norvegia data.
%%%

%%% Load Kapp Norvegia data
load ../data/KN_section/KN_section.mat;

%%% Set true to plot Kapp Norvegia data
plot_KN_data = true;
plot_KN_hyd = false;


%%% Load experiment data
expdir = '../experiments';
% expname = 'WC_seq_onethird_notides_RTOPO2';
% tmin = 4.05;
% tmax = 8.05;
% expname = 'WC_seq_onethird_RTOPO2_unmodEB_alphaV0.5';
% tmin = 6.05;
% tmax = 9.05;
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_strat4e-5';
% tmin = 6.05;
% tmax = 9.05;
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_strat2e-5';
% tmin = 6.05;
% tmax = 9.05;
% expname = 'hires_seq_onethird_RTOPO2';
% tmin = 24.05;
% tmax = 27.05;
% expname = 'hires_seq_onetwelfth_RTOPO2';
% tmin = 6.05;
% tmax = 9.05;
% expname = 'WC_seq_onethird_RTOPO2_dpyc100';
% tmin = 15.05;
% tmax = 18.05;
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_Smin34';
% tmin = 6.05;
% tmax = 9.05;
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_gAlphaU0.5_gAlphaV0.5';
% tmin = 6.05;
% tmax = 9.05;
% expname = 'WC_seq_onethird_RTOPO2_dpyc200_SvaryBCs';
% tmin = 4.05;
% tmax = 7.05;
% expname = 'WC_seq_onethird_RTOPO2_dpyc100_Smin34_Sf34.3';
% tmin = 6.05;
% tmax = 9.05;
expname = 'WC_seq_onethird_RTOPO2_dpyc100_Sf34.3_offshorestrat7e-6';
tmin = 6.05;
tmax = 9.05;

%%% To store diagnostic figures
figdir = fullfile('KN_evaluation',expname);
mkdir(figdir);


%%% Load experiment configuration
loadexp;

%%% Simulation output iteration numbers
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Index of 17W
xidx = find(XC(:,1)>-10,1);

%%% Subsample north of ice shelf cavity, within 1.5 degrees of ice shelf front for TS plots
jmin_TS = find(hFacC(xidx,:,1)>0,1);
lat_icefront = YC(1,jmin_TS);
lat_min_plot = lat_icefront - 1.5;
lat_max_plot = lat_icefront + 2.5;
jmax_TS = find(YC(1,:)>lat_icefront+1.5,1);
jrange_TS = jmin_TS:jmax_TS;
krange_TS = 1:find(squeeze(RC)<-1000,1);

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




%%% Compute thermocline depth
pyc_depth = zeros(Ny,12);
for m = 1:12
  dz_pyc = squeeze(DRC(2:Nr))';
  zz_pyc = squeeze(RF(2:Nr))';
  dTdz = - (theta(:,1:Nr-1,m) - theta(:,2:Nr,m)) ./ dz_pyc;
  dTdz(dTdz<0) = 0;
  for j=1:Ny
    pyc_depth_idx = find((zz_pyc>-1000) & (zz_pyc<-50) & ~isnan(dTdz(j,:)));
    pyc_depth(j,m) = sum(dTdz(j,pyc_depth_idx).*zz_pyc(pyc_depth_idx).*dz_pyc(pyc_depth_idx),2) ./ sum(dTdz(j,pyc_depth_idx).*dz_pyc(pyc_depth_idx),2);
  end
end






%%% Plotting options
framepos = [1000         142        1018        1072];
fontsize = 13;

%%% Generate model plots
for m = 1:length(navg)

  jrange = 1:Ny;
  [ZZ,YY] = meshgrid(RC,YC(1,:));
  
  %%% Straightforward latitude/depth plot of salinity
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
    h = plot(yy,-SHELFICEtopo(xidx,:)/1000,'k','LineWidth',3);  
    h = plot(yy,-bathy(xidx,:)/1000,'k','LineWidth',3);  
  hold off;
  caxis([33.6 34.7]);
  if (m >= 10)
    xlabel('Latitude','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  set(gca,'XLim',[lat_min_plot lat_max_plot]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('haline',20));
    set(handle,'Position',[.93 .05 .01 .9]);
    title(['Years ',num2str(ceil(tmin)),'-',num2str(floor(tmax))],'FontSize',fontsize);
    print('-dpng','-r150',fullfile(figdir,'Salt_lat.png'));
  end



  %%% Straightforward latitude/depth plot of temperature
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
    h = plot(yy,-SHELFICEtopo(xidx,:)/1000,'k','LineWidth',3);  
    h = plot(yy,-bathy(xidx,:)/1000,'k','LineWidth',3);  
    plot(yy,-pyc_depth(:,m)/1000,'-','Color',[.3 .3 .3]);
  hold off;
  caxis([-2.2 1.2])
  if (m >= 10)
    xlabel('Latitude','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  set(gca,'XLim',[lat_min_plot lat_max_plot]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('thermal',20));
    set(handle,'Position',[.93 .05 .01 .9]);
    title(['Years ',num2str(ceil(tmin)),'-',num2str(floor(tmax))],'FontSize',fontsize);
    print('-dpng','-r150',fullfile(figdir,'Temp_lat.png'));
  end


  jidx = find((YC(xidx,:)>lat_icefront) & (bathy(xidx,:)>-4750)); %%% To coincide with KN section
  [ZB,YB] = meshgrid(RC,-bathy(xidx,jidx));

  %%% Bathymetric depth/depth plot of salinity
  figure(6);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  pcolor(YB,-ZB/1000,salt(jidx,:,m));
  shading interp;
  set(gca,'YDir','reverse');  
  caxis([33.6 34.7]);
  if (m >= 10)
    xlabel('Bathymetric depth (m)','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('haline',20));
    set(handle,'Position',[.93 .05 .01 .9]);
    title(['Years ',num2str(ceil(tmin)),'-',num2str(floor(tmax))],'FontSize',fontsize);
    print('-dpng','-r150',fullfile(figdir,'Salt_bath.png'));
  end

  %%% Bathymetric depth/depth plot of temperature
  figure(8);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  pcolor(YB,-ZB/1000,theta(jidx,:,m));
  shading interp;
  set(gca,'YDir','reverse');  
  hold on;
  plot(YB(:,1),-pyc_depth(jidx,m)/1000,'-','Color',[.3 .3 .3])
  hold off;
  caxis([-2.2 1.2])
  if (m >= 10)
    xlabel('Bathymetric depth (m)','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
  end
  set(gca,'YLim',[0 1]);
  if (m == 12)
    handle=colorbar;
    set(handle,'FontSize',fontsize);
    colormap(cmocean('thermal',20));
    set(handle,'Position',[.93 .05 .01 .9]);
    title(['Years ',num2str(ceil(tmin)),'-',num2str(floor(tmax))],'FontSize',fontsize);
    print('-dpng','-r150',fullfile(figdir,'Temp_bath.png'));
  end



  figure(3);
  if (m == 1)
    clf;
    set(gcf,'Position',framepos);
    drawnow;
  end
  subplot(4,3,m);
  scatter(reshape(salt(jrange_TS,krange_TS,m),1,[]),reshape(theta(jrange_TS,krange_TS,m),1,[])); 
  if (m >= 10)
    xlabel('Salinity (g/kg)','interpreter','latex','fontsize',fontsize);
  end
  if (mod(m-1,3)==0)
    ylabel('Temperature ($^\circ$C)','interpreter','latex','fontsize',fontsize);
  end
  axis([33 35 -2.2 1.2])
  if ((m == 12) && ~plot_KN_data)
    title(['Years ',num2str(ceil(tmin)),'-',num2str(floor(tmax))],'FontSize',fontsize);
    print('-dpng','-r150',fullfile(figdir,'TS.png'));
  end

end




if (plot_KN_data)
  
  
  %%% Generate KN plots
  for m = 1:length(navg)
    
    m_KN = mod(m + 5,12)+1;

    [ZZ,YY] = meshgrid(KN_section.pressure,KN_section.distance);
  
    if (plot_KN_hyd)
      
      %%% Straightforward cross-slope distance/depth plot of salinity
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
      clim([33.6 34.7]);
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
        print('-dpng','-r150',fullfile('KN_evaluation','KN_Salt_lat.png'));
      end
    

      %%% Straightforward cross-slope distance/depth plot of temperature
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
        print('-dpng','-r150',fullfile('KN_evaluation','KN_Temp_lat.png'));
      end


      [ZB,YB] = meshgrid(KN_section.pressure,KN_section.bottom_depth);

      %%% Plot salinity in bottom depth/depth space
      figure(7);
      if (m == 1)
        clf;
        set(gcf,'Position',framepos);
        drawnow;
      end
      subplot(4,3,m);
      pcolor(YB(:,:)/1000,ZB(:,:)/1000,KN_section.salt(:,:,m_KN)');
      shading interp;
      set(gca,'YDir','reverse');
      clim([33.6 34.7]);
      if (m >= 10)
        xlabel('Bathymetric depth (m)','interpreter','latex','fontsize',fontsize);
      end
      if (mod(m-1,3)==0)
        ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
      end
      set(gca,'YLim',[0 1]);
      % set(gca,'XLim',[0 100]);
      if (m == 12)
        handle=colorbar;
        set(handle,'FontSize',fontsize);
        colormap(cmocean('haline',20));
        set(handle,'Position',[.93 .05 .01 .9]);
        print('-dpng','-r150',fullfile('KN_evaluation','KN_Salt_bath.png'));
      end

      %%% Plot temperature in bottom depth/depth space
      figure(9);
      if (m == 1)
        clf;
        set(gcf,'Position',framepos);
        drawnow;
      end
      subplot(4,3,m);
      pcolor(YB(:,:)/1000,ZB(:,:)/1000,KN_section.temperature(:,:,m_KN)');
      shading interp;
      set(gca,'YDir','reverse');
      caxis([-2.2 1.2]);
      if (m >= 10)
        xlabel('Bathymetric depth (m)','interpreter','latex','fontsize',fontsize);
      end
      if (mod(m-1,3)==0)
        ylabel('Depth (km)','interpreter','latex','fontsize',fontsize);
      end
      set(gca,'YLim',[0 1]);
      % set(gca,'XLim',[0 100]);
      if (m == 12)
        handle=colorbar;
        set(handle,'FontSize',fontsize);
        colormap(cmocean('thermal',20));
        set(handle,'Position',[.93 .05 .01 .9]);
        print('-dpng','-r150',fullfile('KN_evaluation','KN_Temp_bath.png'));
      end

    end
  
    jrange = 1:26; %%% Within 100km of coast
    krange = 1:51; %%% Above 1000m depth
  
    figure(3);
    % if (m == 1)
    %   clf;
    %   set(gcf,'Position',framepos);
    %   drawnow;
    % end
    subplot(4,3,m);
    hold on;
    scatter(reshape(KN_section.salt(krange,jrange,m_KN),1,[]),reshape(KN_section.temperature(krange,jrange,m_KN),1,[])); 
    hold off;
    % if (m >= 10)
    %   xlabel('Salinity (g/kg)','interpreter','latex','fontsize',fontsize);
    % end
    % if (mod(m-1,3)==0)
    %   ylabel('Temperature ($^\circ$C)','interpreter','latex','fontsize',fontsize);
    % end
    axis([33 35 -2.2 1.2])
    if (m == 12)
      legend('Model','Kapp Norvegia');
      title(['Years ',num2str(ceil(tmin)),'-',num2str(floor(tmax))],'FontSize',fontsize);
      print('-dpng','-r150',fullfile(figdir,'TS.png'));
    end

  end

end