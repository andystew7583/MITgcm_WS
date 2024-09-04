%%%
%%% calcInstShelfHeatBudgetTimeSeries.m
%%%
%%% Computes time series of quantities related to continental shelf heat/salt budgets.
%%%

%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% loadexp;

%%% Index of the upper grid cell face dividing the upper and lower portions
%%% of the water column
zidx_icefront = 25;

%%% Reference surface freezing temperature
theta0 = -1.9;

%%% Reference salinity
salt0 = 34.6;

%%% Set true to deform coordinates in the cavity
deform_cavity = false;

%%% Set true to use barotropic streamfunction as the coordinate system
use_PsiBT = false;

%%% Set true to use depth-averaged temperature as the coordinate system
use_meanT = false;

%%% Define coordinate system for integrating to compute heatfunction
if (use_PsiBT)

  infname = [expname,'_TSfluxes'];
  load(fullfile('products',infname),'uvel_tavg');

  %%% Calculate depth-averaged zonal velocity
  UU = sum(uvel_tavg.*repmat(DRF,[Nx Ny 1]).*hFacW,3);
  clear('uvel_tavg');
  
  %%% Calculate barotropic streamfunction
  Psi = zeros(Nx+1,Ny+1);
  Psi(2:Nx+1,2:Ny+1) = -cumsum(UU.*DYG,2);
  Psi = Psi(1:Nx,1:Ny);
  
  %%% Interpolate to cell centers
  ETA = 0.25*(Psi(1:Nx,1:Ny)+Psi([2:Nx 1],1:Ny)+Psi(1:Nx,[2:Ny 1])+Psi([2:Nx 1],[2:Ny 1]))/1e6;
  
  %%% Streamunction grid for flux calculation
  eta = -2:.1:10;
  Neta = length(eta);

else

  if (use_meanT)

    %%% Load time-mean temperature
    outfname = [expname,'_TSfluxes.mat'];
    load(fullfile('./products',outfname),'theta_tavg');
    
    %%% Coordinate
    ETA = sum(theta_tavg.*DRF.*hFacC,3)./sum(DRF.*hFacC,3);
    eta = -2.6:0.025:0;
    Neta = length(eta);

  else

    ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
    eta = -9:.1:11;
    Neta = length(eta);

    %%% Bounds and horizontal mask for heat budget volume
    eta_icefront = -1.1;
    eta_shelfbreak = 3.5;             

  end

end

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
%%% For daily/12-hourly outputs
dumpStart = 1578240;
dumpStep = 86400/2/60;
nDumps = 731;
dumpIters = dumpStart:dumpStep:dumpStart+(nDumps-1)*dumpStep;
itersToRead = dumpIters;
times = dumpIters*deltaT;
Ntime = length(itersToRead);


%%% Calendar information
days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays = [0 cumsum(days_per_month)];
Nmonths = length(days_per_month);





%%% Write to output file
outfname = [expname,'_InstShelfHeatBudget.mat'];
load(fullfile('products',outfname));

startyear = 2011;
rho0 = 1000;
Cp = 4000;
times = (((1:12) - 0.5) / 12)*t1year;
figure(106);
plot(startyear+times/t1year,(hflux_icefront-mflux_icefront*theta0)*rho0*Cp/1e12,'LineWidth',1.5);
hold on;
plot(startyear+times/t1year,-(hflux_shelfbreak-mflux_shelfbreak*theta0)*rho0*Cp/1e12,'LineWidth',1.5);
plot(startyear+times/t1year,-(wt-wm*theta0)*rho0*Cp/1e12,'LineWidth',1.5);
plot(startyear+times/t1year,ttend*rho0*Cp/1e12,'LineWidth',1.5);
plot(startyear+times/t1year,(hflux_icefront-hflux_shelfbreak-wt - ttend)*rho0*Cp/1e12,'k--','LineWidth',1.5);
plot(startyear+times/t1year,0*times,'k:','LineWidth',1.5);
hold off
legend('Q_c_a_v_i_t_y','Q_s_h_e_l_f','Q_v_e_r_t','Tendency','Sum');
xlabel('Year');
ylabel('Heat flux (TW)');
set(gca,'Position',[0.05 0.1 0.9 0.85]);
set(gca,'FontSize',14);


startyear = 2011;
rho0 = 1000;
Cp = 4000;
times = (((1:12) - 0.5) / 12)*t1year;
figure(107);
plot(startyear+times/t1year,(sflux_icefront-mflux_icefront*salt0)*rho0/1e9,'LineWidth',1.5);
hold on;
plot(startyear+times/t1year,-(sflux_shelfbreak-mflux_shelfbreak*salt0)*rho0/1e9,'LineWidth',1.5);
plot(startyear+times/t1year,-(ws-wm*salt0)*rho0/1e9,'LineWidth',1.5);
plot(startyear+times/t1year,stend*rho0/1e9,'LineWidth',1.5);
plot(startyear+times/t1year,(sflux_icefront-sflux_shelfbreak-ws - stend)*rho0/1e9,'k--','LineWidth',1.5);
plot(startyear+times/t1year,0*times,'k:','LineWidth',1.5);
hold off
legend('Q_c_a_v_i_t_y','Q_s_h_e_l_f','Q_v_e_r_t','Tendency','Sum');
xlabel('Year');
ylabel('Salt flux (Gg/s)');
set(gca,'Position',[0.05 0.1 0.9 0.85]);
set(gca,'FontSize',14);