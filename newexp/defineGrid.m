%%%
%%% defineGrid.m
%%%
%%% Defines a Weddell Sea model grid for MITgcm.
%%%

%%% Choose resolution
res_fac = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MPI parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% N.B. nPx_max and nPy_max are used to define the grid size, ensuring
%%% that the grid can be split into nPx_max x nPy_max processors. This
%%% allows smaller nPx and nPy (divisors of nPx_max and nPy_max) to be
%%% selected without modifying the grid size.
switch (res_fac)
  case 3
    nPx_max = 16; %%% max. no. of processors in x-direction
    nPy_max = 8; %%% max. no. of processors in y-direction
    nPx = 16; %%% no. of processors in x-direction
    nPy = 8; %%% no. of processors in y-direction
  case 6
    nPx_max = 36; %%% max. no. of processors in x-direction
    nPy_max = 48; %%% max. no. of processors in y-direction
%     nPx = 18; %%% no. of processors in x-direction
%     nPy = 12; %%% no. of processors in y-direction
    nPx = 36; %%% no. of processors in x-direction
    nPy = 24; %%% no. of processors in y-direction
%     nPx = 9; %%% no. of processors in x-direction
%     nPy = 6; %%% no. of processors in y-direction
  case 12
    nPx_max = 36; %%% max. no. of processors in x-direction
    nPy_max = 48; %%% max. no. of processors in y-direction
    nPx = 18; %%% no. of processors in x-direction
    nPy = 12; %%% no. of processors in y-direction
%     nPx = 36; %%% no. of processors in x-direction
%     nPy = 48; %%% no. of processors in y-direction
%     nPx = 9; %%% no. of processors in x-direction
%     nPy = 6; %%% no. of processors in y-direction
  case 24
    nPx_max = 72; %%% max. no. of processors in x-direction
    nPy_max = 48; %%% max. no. of processors in y-direction
    nPx = 18; %%% no. of processors in x-direction
    nPy = 12; %%% no. of processors in y-direction
  case 30
    nPx_max = 36; %%% max. no. of processors in x-direction
    nPy_max = 48; %%% max. no. of processors in y-direction
    nPx = 12; %%% no. of processors in x-direction
    nPy = 6; %%% no. of processors in y-direction
  case 32
    nPx_max = 36; %%% max. no. of processors in x-direction
    nPy_max = 48; %%% max. no. of processors in y-direction
    nPx = 12; %%% no. of processors in x-direction
    nPy = 6; %%% no. of processors in y-direction
  case 48
    nPx_max = 144; %%% max. no. of processors in x-direction
    nPy_max = 96; %%% max. no. of processors in y-direction
    nPx = 18; %%% no. of processors in x-direction
    nPy = 12; %%% no. of processors in y-direction
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Data format parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ieee='b';
prec='real*8';
realdigits = 8;
realfmt=['%.',num2str(realdigits),'e'];
months = 12; 

%%%% Option for finer grid %%
FinerVerticalGrid = 0;

%%%% which directory to use
gendir = '/data3';

%%% Set true for open southern/western boundaries
use_OB_SW = false;

%%% Set true to set SSH on boundaries (ONLY IF USING NONLINEAR FREE
%%% SURFACE)
use_OB_eta = false;

%%% By default we have active ice shelves
use_shelfice = true;

%%% Set true if we are nested within a lower-resolution simulation
obcs_nest = false;

%%% Name under which to store this grid
switch (res_fac)
  case 3
    grid_name = 'one_third';
  case 6
    grid_name = 'one_sixth';
  case 12
    grid_name = 'one_twelfth';
  case 24
    grid_name = 'one_twentyfourth';
    obcs_nest = true;
  case 30
    grid_name = 'one_thirtieth';
    use_OB_SW = true;
%     use_shelfice = false;
    obcs_nest = true;
  case 32
    grid_name = 'one_thirtysecond';
    use_OB_SW = true;
%     use_shelfice = false;
    obcs_nest = true;
  case 48
    grid_name = 'one_fortyeighth';
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DOMAIN LIMITS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Maximum depth
% H = 5400;
H = 5000;
% H = 4500;


xmin = -83;
ymin = -83.5;

switch (res_fac)
  case 24    
    %%% Southern WS subdomain
    xmax = -20;
    ymax = -70;
  case 30    
    %%% Filchner overflow domain
    xmax = -26;
    ymax = -71;
    xmin = -48;
    ymin = -77;
  case 32    
    %%% Riiser-Larsen domain
    xmax = -12;
    ymax = -69;
    xmin = -32;
    ymin = -75;
  case 48
    %%% Southern WS subdomain
    xmax = -20;
    ymax = -70;
  otherwise
    %%% Full MITgcm_WS domain
    xmax = 21;
    ymax = -64;
end

EXF_xmin = -83;
EXF_ymin = -83.5;
EXF_xmax = 21;
EXF_ymax = -64;

zmin = 0;
zmax = -H;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ATMOSPHERIC FORCING CYCLE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
%%% Seconds in one hour
t1min = 60;
%%% Seconds in one hour
t1hour = 60*t1min;
%%% hours in one day
t1day = 24*t1hour;
%%% Metres in one kilometre
m1km = 1000;


%%% High-res subdomain needs to start from 2008 because we don't have all
%%% lower-res data required for 2007
switch (res_fac)
  case 24
    start_year = 2008;
    endyr = 2014;
  case 30
    start_year = 2008;
    endyr = 2014;
  case 32
    start_year = 2008;
    endyr = 2009;
  case 48
    start_year = 2011;
    endyr = 2012;
  otherwise
    start_year = 2007;
    endyr = 2015;
end
base_year = 2006;

start_month = 1;
end_month = 12;

%%%%  total years, months
Nmon = 12;
Nyears = endyr-start_year+1;
Nmonths = Nmon*Nyears;

%%%% Finding if current year is a leap year
%%% N.B. only works because we happen to be using the period 2007-2015
Ndays = 0;
is_leap_year = zeros(Nyears,1);
for i=1:Nyears
  currentyear = start_year + i - 1;
  if (mod(currentyear,4)==0) %%% Only works because of the particular period we're looking at
    is_leap_year(i) = 1;
    Ndays = Ndays + 366;
  else
    Ndays = Ndays + 365;
  end
end
simTime = Ndays*86400;
t1year = simTime/Nyears;
t1month = t1year/12;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Grid size of the desired run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nx = (xmax-xmin)*res_fac;
Nx = floor(Nx/nPx_max)*nPx_max;

%%% Physical parameters
g = 9.81; %%% Gravity
Omega = 2*pi*366/365/86400; %%% Earth rotation rate
Rp = 6400000; %%% Planetary radius (approx.)
rho0 = 1027.5; %%% Reference density
Pa1dbar = 1e4; %%% Pascals in 1 decibar

%%% Grid cell thickness constraints
hFacMin = 0.1;
hFacMinDr = 10;

%%% Domain lengths
Lx = xmax-xmin;
Ly = ymax-ymin;
Lz = zmax-zmin;

EXF_Lx = EXF_xmax-EXF_xmin;
EXF_Ly = EXF_ymax-EXF_ymin;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Zonal grid %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dmxg = Lx/Nx*ones(Nx,1);
xmc = xmin+dmxg/2:dmxg:xmax-dmxg/2; 
xmc = xmc';
xvmg = xmc;
xumg = xmin:dmxg:xmax-dmxg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SOSE GRID FOR TESTING PURPOSES%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% xmin = -83;
% xmax = -3;
% Nx = 480;
% dx = (Lx/Nx*(1/6))*ones(Nx,1);
% Ny = 48;
% dy = Ly/Ny*(1/6)*ones(Ny,1);

% 
% ymin = -77.9583;
% ymax = -24.71;
% Ny = 320;
% dy = (Ly/Ny*(1/6))*ones(Ny,1);

% dz = [10.0 11.0 12.0 13.0 14.0 16.0 18.0 20.0 23.0 26.0 ...
%      29.0 33.0 37.0 42.0 48.0 55.0 63.0 72.0 82.0 ...
%      94.0 108.0 124.0 142.0 163.0 187.0 215.0 247.0 ...
%      284.0 262.0 250.0 250.0 250.0 250.0 250.0 250.0 ...
%      250.0 250.0 250.0 250.0 250.0 250.0 250.0];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Latitudinal grid %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% First pass - create approximate dmyg
yval = ymin;
dmyg = [];
while (yval < ymax)
  dyval = cosd(yval)*dmxg(1);
  dmyg = [dmyg dyval];
  yval = yval + dyval;
end

%%% Second pass - adjust number of grid points to be a multiple of nPy
Ny = floor((length(dmyg)/nPy_max))*nPy_max;
dmyg = dmyg(1:Ny);

%%% Third pass - resize dmyg so that its total length is equal to the
%%% domain size
dmyg = dmyg*(Ly/sum(dmyg));

%%% Other y-grids
yvmg = ymin + cumsum([0 dmyg(1:Ny-1)]);
ymc = yvmg + dmyg/2;
dmyc = ymc(1:Ny) - ymc([Ny 1:Ny-1]);
yumg = ymc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TILE GRID SIZES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


sNx = (Nx)/nPx;  %%% no. of x-gridpoints per tile
sNy = (Ny)/nPy;  %%% no. of y-gridpoints per tile



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Grid for EXF Interpolation (Default will be from 1/3rd degree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Grid sizes
EXF_nPx = 16; %%% Need to be defined to make sure EXF grid matches the grid 
EXF_nPy = 8; %%% used to generate the surface forcing dataset
EXF_Nx = (EXF_xmax-EXF_xmin)*3;
EXF_Nx = floor(EXF_Nx/EXF_nPx)*EXF_nPx;

%%% Time frame
EXF_start_year = 2007;
EXF_endyr = 2015;
EXF_Nyears = EXF_endyr-EXF_start_year+1;
EXF_Nmonths = Nmon*Nyears;

%%%% Finding if current year is a leap year
%%% N.B. only works because we happen to be using the period 2007-2015
EXF_Ndays = 0;
is_leap_year = zeros(EXF_Nyears,1);
for i=1:EXF_Nyears
  currentyear = EXF_start_year + i - 1;
  if (mod(currentyear,4)==0) %%% Only works because of the particular period we're looking at
    is_leap_year(i) = 1;
    EXF_Ndays = EXF_Ndays + 366;
  else
    EXF_Ndays = EXF_Ndays + 365;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Zonal grid %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EXF_dmxg = EXF_Lx/EXF_Nx*ones(EXF_Nx,1);
EXF_xmc = EXF_xmin+EXF_dmxg/2:EXF_dmxg:EXF_xmax-EXF_dmxg/2; 
EXF_xmc = EXF_xmc';
EXF_xvmg = EXF_xmc;
EXF_xumg = EXF_xmin:EXF_dmxg:EXF_xmax-EXF_dmxg;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Latitudinal grid %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% First pass - create approximate dmyg
EXF_yval = EXF_ymin;
EXF_dmyg = [];
while (EXF_yval < EXF_ymax)
  EXF_dyval = cosd(EXF_yval)*EXF_dmxg(1);
  EXF_dmyg = [EXF_dmyg EXF_dyval];
  EXF_yval = EXF_yval + EXF_dyval;
end

%%% Second pass - adjust number of grid points to be a multiple of nPy
EXF_Ny = floor((length(EXF_dmyg)/EXF_nPy))*EXF_nPy;
EXF_dmyg = EXF_dmyg(1:EXF_Ny);

%%% Third pass - resize dmyg so that its total length is equal to the
%%% domain size
EXF_dmyg = EXF_dmyg*(EXF_Ly/sum(EXF_dmyg));

%%% Other y-grids
EXF_yvmg = ymin + cumsum([0 EXF_dmyg(1:EXF_Ny-1)]);
EXF_ymc = EXF_yvmg + EXF_dmyg/2;
EXF_dmyc = EXF_ymc(1:EXF_Ny) - EXF_ymc([EXF_Ny 1:EXF_Ny-1]);
EXF_yumg = EXF_ymc;




%%% Meshgrid for interpolating surface fields

[EXF_XMC,EXF_YMC] = meshgrid(EXF_xmc,EXF_ymc);






   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Vertical grid%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FinerVerticalGrid == 1

    z0 = 0;
    z1 = 250;
    z2 = 1800;
    z3 = 4000;
    z4 = 5000;
    dz0 = 1;
    dz1 = 5; 
    dz2 = 30;
    dz3 = 80;
    dz4 = 150;
    dz5 = 220;
    N0 = 1;
    N1 = 10; 
    N2 = 90;
    N3 = 35;
    N4 = 10;
    N5 = 4;
    nn_c = cumsum([N0 N1 N2 N3 N4 N5]);
end
if FinerVerticalGrid == 0
    z0 = 0;
    z1 = 500;
%     z1 = 250;
    z2 = 1800;
    z3 = 4000;
    z4 = 5000;
    dz0 = 3;
%     dz0 = 1;
    dz1 = 15; 
    dz2 = 30;
    dz3 = 80;
    dz4 = 150;
    dz5 = 220;
%     N1 = 10;
    N1 = 20;
    N2 = 59;
%     N2 = 69;
    N3 = 35;
    N4 = 10;
    N5 = 4;
    nn_c = cumsum([1 N1 N2 N3 N4 N5]);
end
%%% 1/3 degree configuration

% z0 = 0;
% z1 = 250;
% z2 = 1800;
% z3 = 4000;
% z4 = 5000;
% dz0 = 10;
% dz1 = 30;
% dz2 = 60;
% dz3 = 160;
% dz4 = 220;
% dz5 = 320;
% N1 = 10;
% N2 = 30;
% N3 = 20;
% N4 = 5;
% N5 = 2;




% nn_c = cumsum([N0 N1 N2 N3 N4 N5]);
% nn_c = cumsum([1 N1 N2 N3 N4 N5]);

dz_c = [dz0 dz1 dz2 dz3 dz4 dz5];
nn = 1:(N1+N2+N3+N4+N5+1);
dz = interp1(nn_c,dz_c,nn,'pchip');

zz = -cumsum((dz+[0 dz(1:end-1)])/2);
Nz = length(zz);

%%% Heights of cell vertical faces
zzf = -cumsum([0 dz]);



%%%% vertical coordinate of center of cell
mrc = zz;
mrc = mrc';

%%% Vertical coordinate of cell faces
mrf = zzf';

%%% For convenience
Nr = Nz;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% TWIST FUNCTION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creating meshgrid of tracer lat/lons for our defined grid
    [XMC,YMC] = meshgrid(xmc,ymc);

% %Meshgrid for U velocity
    [XUMG,YUMG] = meshgrid(xumg,yumg);
%  
% %Meshgrid for V Velocity
     [XVMG,YVMG] = meshgrid(xvmg,yvmg);
%     

% %Creating meshgrid for interpolation of new grid vertical level (tracers)
    [Xm3C,Ym3C,RM3C] = meshgrid(xmc,ymc,mrc);
%  
% %Creating meshgrid for interp of new vertical grid for U vel
%     [XUmG,YUmG,RM3U] = meshgrid(xumg,yumg,mrc);
% 
% %Creating meshgrid for interp of new vertical grid for V vel
%     [XVmG,YVmG,RM3V] = meshgrid(xvmg,yvmg,mrc);
% 


%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT FILES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

bathyFile = 'bathyFile.bin';
hydrogThetaFile = 'hydrogThetaFile.bin';
hydrogSaltFile = 'hydrogSaltFile.bin';
pSurfInitFile = 'pSurfInitFile.bin';
uVelInitFile = 'uVelInitFile.bin';
vVelInitFile = 'vVelInitFile.bin';
SHELFICEtopoFile = 'SHELFICEtopoFile.bin';
SHELFICEloadAnomalyFile = 'SHELFICEloadAnomalyFile.bin';



%%%%%%%% SEA ICE

HeffFile = 'HeffFile.bin';
AreaFile = 'AreaFile.bin';
HsaltFile = 'HsaltFile.bin';
HsnowFile = 'HsnowFile.bin';
uIceFile = 'uIceFile.bin';
vIceFile = 'vIceFile.bin';

%%%%%%%% SFC FORCINGS

zwind = 'uwindfile.bin';
mwind = 'vwindfile.bin';
aTemp = 'atempfile.bin';
aPrecip = 'precipfile.bin';
aLW = 'lwdownfile.bin';
aSW = 'swdownfile.bin';
anewAQ = 'aqhfile.bin';
pressure = 'apressurefile.bin';

%%%%%%%% OBCS 

OBNtFile = 'OBNtFile.bin';
OBEtFile = 'OBEtFile.bin';
OBWtFile = 'OBWtFile.bin';
OBStFile = 'OBStFile.bin';

OBNsFile = 'OBNsFile.bin';
OBEsFile = 'OBEsFile.bin';
OBWsFile = 'OBWsFile.bin';
OBSsFile = 'OBSsFile.bin';

OBNuFile = 'OBNuFile.bin';
OBEuFile = 'OBEuFile.bin';
OBWuFile = 'OBWuFile.bin';
OBSuFile = 'OBSuFile.bin';

OBNvFile = 'OBNvFile.bin';
OBEvFile = 'OBEvFile.bin';
OBWvFile = 'OBWvFile.bin';
OBSvFile = 'OBSvFile.bin';

OBNetaFile = 'OBNetaFile.bin';
OBEetaFile = 'OBEetaFile.bin';
OBWetaFile = 'OBWetaFile.bin';
OBSetaFile = 'OBSetaFile.bin';

OBNaFile = 'OBNaFile.bin';
OBEaFile = 'OBEaFile.bin';
OBWaFile = 'OBWaFile.bin';
OBSaFile = 'OBSaFile.bin';

OBNhFile = 'OBNhFile.bin';
OBEhFile = 'OBEhFile.bin';
OBWhFile = 'OBWhFile.bin';
OBShFile = 'OBShFile.bin';

OBNsnFile = 'OBNsnFile.bin';
OBEsnFile = 'OBEsnFile.bin';
OBWsnFile = 'OBWsnFile.bin';
OBSsnFile = 'OBSsnFile.bin';

OBNuiceFile = 'OBNuiceFile.bin';
OBEuiceFile = 'OBEuiceFile.bin';
OBWuiceFile = 'OBWuiceFile.bin';
OBSuiceFile = 'OBSuiceFile.bin';

OBNviceFile = 'OBNviceFile.bin';
OBEviceFile = 'OBEviceFile.bin';
OBWviceFile = 'OBWviceFile.bin';
OBSviceFile = 'OBSviceFile.bin';

OBEslFile = 'OBEslFile.bin';
OBWslFile = 'OBWslFile.bin';
OBNslFile = 'OBNslFile.bin';
OBSslFile = 'OBSslFile.bin';


 
AngleSN = zeros(size(XMC));
AngleCS = ones(size(XMC));



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% WRITE DATA %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% General Input Directory %%%%%%%

%Input file path

inputfolder = './DEFAULTS/input';
inputconfigdir = ['./InputConfigs/' grid_name];
if (~exist(inputconfigdir))
  mkdir(inputconfigdir);
end








