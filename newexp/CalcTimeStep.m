%%%
%%% CalcTimeStep.m
%%%
%%% Estimates the time step required for MITgcm runs.
%%%


%%% Load model grid
defineGrid
inputpath = './DEFAULTS/input'; 
% inputpath = './Expanded'; 




%%% GSW scripts
addpath gsw_matlab_v3_05_4/ 
addpath gsw_matlab_v3_05_4/library/
addpath gsw_matlab_v3_05_4/thermodynamics_from_t/



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DEFORMATION RADIUS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Load initial theta and salt
ModelTwistedTheta = zeros(Nx,Ny,Nr);
fid = fopen(fullfile(inputpath,hydrogThetaFile),'r','b');
for k=1:Nr
  ModelTwistedTheta(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);  

TwistedSaltPretzal = zeros(Nx,Ny,Nr);
fid = fopen(fullfile(inputpath,hydrogSaltFile),'r','b');
for k=1:Nr
  TwistedSaltPretzal(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);

%%% Read bathymetry and ice shelf depth
fid = fopen(fullfile(inputpath,bathyFile),'r','b');
h = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);

fid = fopen(fullfile(inputpath,SHELFICEtopoFile),'r','b');
hice = fread(fid,[Nx Ny],'real*8'); 
fclose(fid);


%%% Check Brunt-Vaisala frequency using full EOS  
Cig = zeros(Nx,Ny);
Rd = zeros(Nx,Ny);

for i=1:Nx
  for j=1:Ny
    
    %%% Ignore dry columns
    if (hice(i,j)-h(i,j) <= 0)
      continue;
    end

    %%% Calculate N^2
    pp = - zz;
    ss = squeeze(TwistedSaltPretzal(i,j,:))';
    pt = squeeze(ModelTwistedTheta(i,j,:))';        
    SA = gsw_SA_from_SP(ss,pp,xmc(i),ymc(j));  
    CT = gsw_CT_from_pt(SA,pt);
    [N2 pp_mid] = gsw_Nsquared(SA,CT,pp);
    dz_mid = zz(1:end-1)-zz(2:end);

    %%% Calculate internal wave speed and first Rossby radius of deformation
    N = sqrt(N2);
    N(N2<0) = 0;
    Cig(i,j) = 0;
    for k=1:length(dz_mid);
%       Cig = 0   
      if (zz(k) > h(i,j) && zz(k) < hice(i,j))         
         Cig(i,j) = Cig(i,j) + N(k)*dz_mid(k);
      end
    end
    f0 = 2*Omega*sind(ymc(j));
    Rd(i,j) = Cig(i,j)./(pi*abs(f0));

  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATE TIME STEP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% These estimates are in no way complete, but they give at least some
%%% idea of the time step needed to keep things stable. In complicated 
%%% simulations, preliminary tests may be required to estimate the
%%% parameters used to calculate these time steps.        

%%% Upper bound for absolute horizontal fluid velocity (m/s)
%%% At the moment this is just an estimate
Umax = 2.0;  
deltaT = NaN*ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny

    %%% Not clear that this is a good idea - grid spacing is smallest under the
   %%% ice
%     if ((hice(i,j)-h(i,j) <= 0) || (hice(i,j)<0)) 

    %%% Ignore dry columns
    if ((hice(i,j)-h(i,j) <= 0))
      continue;
    end
    
    %%% Grid spacings
    dx = dmxg(i)*2*pi/360*Rp*cosd(ymc(j));
    dy = dmyg(j)*2*pi/360*Rp;

    %%% Gravity wave CFL

    %%% Max gravity wave speed 
    cmax = Cig(i,j);
    %%% Max gravity wave speed using total ocean depth
    cgmax = Umax + cmax;
    %%% Advective CFL
    deltaT_adv = min([0.5*dx/cmax,0.5*dy/cmax]);
    %%% Gravity wave CFL
    deltaT_gw = min([0.5*dx/Umax,0.5*dy/Umax]);
    %%% CFL time step based on full gravity wave speed
    deltaT_fgw = min([0.5*dx/cgmax,0.5*dy/cgmax]);

    %%% Time step size  
    deltaT(i,j) = min([deltaT_fgw deltaT_gw deltaT_adv]);

  end
end

%%% Theoretical minimum time step
deltaT = round(nanmin(nanmin(deltaT)));

%%% Add a safety factor to the actual time step
deltaT = .8*deltaT;

%%% Write to file
fid = fopen('TIME_STEP','w');
fprintf(fid,'%d',deltaT);
fclose(fid);