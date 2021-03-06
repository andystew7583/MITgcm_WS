%%%
%%% loadexp.m
%%%
%%% Loads an experiment's parameters into memory. The variables 'expname'
%%% and 'expdir' must be set prior to running this script. This is a bit
%%% clumsy, but with so many variables to load it's the only practical way
%%% of doing this.
%%%

%%%
%%% Example names/paths of experiments to read from:
%%%
%%% expname = 'OT_w0.5_c7';
%%% expdir = './';
%%%
%%% expname = 'OT_w0.5_c7';
%%% expdir = './OT_batch1';
%%%


%%% Set up experiment paths
exppath = fullfile(expdir,expname);
inputpath = fullfile(exppath,'input');
resultspath = fullfile(exppath,'results');

%%% Path to MITgcm matlab scripts
addpath ../utils/matlab/

%%% Load parameters used for this experiment
run(fullfile(inputpath,'params.m'));



%%%
%%% GRID
%%%

%%% Grid dimensions (not specified explicitly in params.m)
Nx = length(delX);
Ny = length(delY);
Nr = length(delR);

%%% Domain dimensions
Lx = sum(delX);
Ly = sum(delY);
H = sum(delR);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ATMOSPHERIC FORCING CYCLE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Beginning Data from 2007  ( first full year we have)
base_year = 2006;
start_year = 2007;
endyr = 2015;

start_month = 1;
end_month = 12;

%%%%  total years, months
Nmon = 12;
Nyears = 9;
Nmonths = 108;

%%%% Finding if current year is a leap year
%%% N.B. only works because we happen to be using the period 2007-2015
Ndays = 0;
is_leap_year = zeros(Nyears,1);
for i=1:Nyears
  if (mod(i+2,4)==0) %%% At i=2, the year is 2008, e.g.
    is_leap_year(i) = 1;
    Ndays = Ndays + 366;
  else
    Ndays = Ndays + 365;
  end
end

%%% Useful measures of time
t1day = 86400;
t1year = Ndays*t1day/Nyears;
t1month = t1year/12;

days_start = 1;
days_end = Ndays;





%%% Modify to ensure we catch all the time steps in the event that the
%%% simulation has been restarted
startTime = nIter0*deltaT;
nTimeSteps = ceil((endTime-startTime)/deltaT);






%%% Other physical parameters
rho0 = 1000;
Cp = 4e3;


atmos_forcing = 0;
open_boundaries = 0;
initialconditions = 0;
tides = 0;


%%% Load data files

fid = fopen(fullfile(inputpath,'bathyFile.bin'),'r','b');
    bathy = fread(fid,[Nx Ny],'real*8');
    fclose(fid);

    fid = fopen(fullfile(inputpath,SHELFICEtopoFile),'r','b');
    SHELFICEtopo = fread(fid,[Nx Ny],'real*8');
    fclose(fid);
    
    
    
    

if open_boundaries ==1

    OBEt = zeros(Ny,Nr,12);
    fid = fopen(fullfile(inputpath,OBEtFile),'r','b');
    for k=1:12
        OBEt(:,:,k) = fread(fid,[Ny Nr],'real*8');
    end
    fclose(fid);
% 
    OBEv = zeros(Ny,Nr,12);
    fid = fopen(fullfile(inputpath,OBEvFile),'r','b');
    for k=1:12
     OBEv(:,:,k) = fread(fid,[Ny Nr],'real*8');
    end
    fclose(fid);
% 

    OBEs = zeros(Ny,Nr,12);
    fid = fopen(fullfile(inputpath,OBEsFile),'r','b');
    for k=1:12
     OBEs(:,:,k) = fread(fid,[Ny Nr],'real*8');
    end
    fclose(fid);


    OBNu = zeros(Nx,Nr,12);
    fid = fopen(fullfile(inputpath,OBNuFile),'r','b');
    for k=1:12
     OBNu(:,:,k) = fread(fid,[Nx Nr],'real*8');
    end
    fclose(fid);
    
    
    OBNv = zeros(Nx,Nr,12);
    fid = fopen(fullfile(inputpath,OBNvFile),'r','b');
    for k=1:12
        OBNv(:,:,k) = fread(fid,[Nx Nr],'real*8');
    end
    fclose(fid);

    OBEu = zeros(Ny,Nr,12);
    fid = fopen(fullfile(inputpath,OBEuFile),'r','b');
    for k=1:12
        OBEu(:,:,k) = fread(fid,[Ny Nr],'real*8');
    end
    fclose(fid);

    OBNt = zeros(Nx,Nr,24);
    fid = fopen(fullfile(inputpath,OBNtFile),'r','b');
    for k=1:12
        OBNt(:,:,k) = fread(fid,[Nx Nr],'real*8');
    end
    fclose(fid);
    
       
    OBNs = zeros(Nx,Nr,12);
    fid = fopen(fullfile(inputpath,OBNsFile),'r','b');
    for k=1:12
        OBNs(:,:,k) = fread(fid,[Nx Nr],'real*8');
    end
    fclose(fid);
    
    
    OBNa = zeros(Nx,12);
    fid = fopen(fullfile(inputpath,OBNaFile),'r','b');
    for k=1:12
        OBNa(:,k) = fread(fid,[Nx],'real*8');
    end
    fclose(fid);
    
    OBEa = zeros(Ny,12);
    fid = fopen(fullfile(inputpath,OBEaFile),'r','b');
    for k=1:12
        OBEa(:,k) = fread(fid,Ny,'real*8');
    end
    fclose(fid);
    
    OBNh = zeros(Nx,12);
    fid = fopen(fullfile(inputpath,OBNhFile),'r','b');
    for k=1:12
        OBNh(:,k) = fread(fid,[Nx],'real*8');
    end
    fclose(fid);
    
    OBEh = zeros(Ny,12);
    fid = fopen(fullfile(inputpath,OBEhFile),'r','b');
    for k=1:12
        OBEh(:,k) = fread(fid,Ny,'real*8');
    end
    fclose(fid);
    
    OBNuice = zeros(Nx,12);
    fid = fopen(fullfile(inputpath,OBNuiceFile),'r','b');
    for k=1:12
        OBNuice(:,k) = fread(fid,[Nx],'real*8');
    end
    fclose(fid);
    
    OBEuice = zeros(Ny,12);
    fid = fopen(fullfile(inputpath,OBEuiceFile),'r','b');
    for k=1:12
        OBEuice(:,k) = fread(fid,Ny,'real*8');
    end
    fclose(fid);
    
    OBNvice = zeros(Nx,12);
    fid = fopen(fullfile(inputpath,OBNviceFile),'r','b');
    for k=1:12
        OBNvice(:,k) = fread(fid,[Nx],'real*8');
    end
    fclose(fid);
    
    OBEvice = zeros(Ny,12);
    fid = fopen(fullfile(inputpath,OBEviceFile),'r','b');
    for k=1:12
        OBEvice(:,k) = fread(fid,Ny,'real*8');
    end
    fclose(fid);
    
    OBEam = zeros(Ny,10);
    fid = fopen(fullfile(inputpath,OBEamFile),'r','b');
    for k=1:10
        OBEam(:,k) = fread(fid,Ny,'real*8');
    end
    fclose(fid);
%     
%     OBNeta = zeros(Nx,12);
%     fid = fopen(fullfile(inputpath,OBNetaFile),'r','b');
%     for k=1:12
%         OBNeta(:,k) = fread(fid,[Nx],'real*8');
%     end
%     fclose(fid);
%     
%     OBEeta = zeros(Ny,12);
%     fid = fopen(fullfile(inputpath,OBEetaFile),'r','b');
%     for k=1:12
%         OBEeta(:,k) = fread(fid,Ny,'real*8');
%     end
%     fclose(fid);
end


if initialconditions ==1
    

    fid = fopen(fullfile(inputpath,HeffFile),'r','b');
    heff = fread(fid,[Nx Ny],'real*8');
    fclose(fid);
    
    fid = fopen(fullfile(inputpath,uIceFile),'r','b');
    uIce = fread(fid,[Nx Ny],'real*8');
    fclose(fid);   


    hydrogTheta = zeros(Nx,Ny,Nr);
    fid = fopen(fullfile(inputpath,hydrogThetaFile),'r','b');
    for k=1:Nr
     hydrogTheta(:,:,k) = fread(fid,[Nx Ny],'real*8');
    end
    fclose(fid);
% 
    hydrogSalt = zeros(Nx,Ny,Nr);
    fid = fopen(fullfile(inputpath,hydrogSaltFile),'r','b');
    for k=1:Nr
        hydrogSalt(:,:,k) = fread(fid,[Nx Ny],'real*8');
    end
    fclose(fid);

    fid = fopen(fullfile(inputpath,'uIceFile.bin'),'r','b');
    SIuice = fread(fid,[Nx Ny],'real*8');
    fclose(fid);


    uVelInit = zeros(Nx,Ny,Nr);
    fid = fopen(fullfile(inputpath,uVelInitFile),'r','b');
    for k=1:Nr
        uVelInit(:,:,k) = fread(fid,[Nx Ny],'real*8');
    end
    fclose(fid);

    vVelInit = zeros(Nx,Ny,Nr);
    fid = fopen(fullfile(inputpath,vVelInitFile),'r','b');
    for k=1:Nr
         vVelInit(:,:,k) = fread(fid,[Nx Ny],'real*8');
    end
    fclose(fid);
    
    fid = fopen(fullfile(inputpath,pSurfInitFile),'r','b');
    pSurfInit = fread(fid,[Nx Ny],'real*8');
    fclose(fid);

    fid = fopen(fullfile(inputpath,SHELFICEloadAnomalyFile),'r','b');
    SHELFICEloadAnomaly = fread(fid,[Nx Ny],'real*8');
    fclose(fid);
    
    fid = fopen(fullfile(inputpath,HeffFile),'r','b');
    SIHeffInit = fread(fid,[Nx Ny],'real*8');
    fclose(fid);
    
    fid = fopen(fullfile(inputpath,AreaFile),'r','b');
    SIAreaInit = fread(fid,[Nx Ny],'real*8');
    fclose(fid);
    
    fid = fopen(fullfile(inputpath,uIceFile),'r','b');
    SIUvelInit = fread(fid,[Nx Ny],'real*8');
    fclose(fid);
    
    fid = fopen(fullfile(inputpath,vIceFile),'r','b');
    SIVvelInit = fread(fid,[Nx Ny],'real*8');
    fclose(fid);
end


if atmos_forcing ==1
  
    EXF_Nx = atemp_nlon;
    EXF_Ny = atemp_nlat;
    
    uwind = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
    fid = fopen(fullfile(inputpath,uwindfile),'r','b');
    for k=1:days_start:days_end
        uwind(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    end
    fclose(fid);

    vwind = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
    fid = fopen(fullfile(inputpath,vwindfile),'r','b');
    for k=days_start:days_end
        vwind(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    end
    fclose(fid);

    atemp = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
    fid = fopen(fullfile(inputpath,atempfile),'r','b');
    for k=days_start:days_end
        atemp(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    end
    fclose(fid);

    lw = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
    fid = fopen(fullfile(inputpath,lwdownfile),'r','b');
    for k=days_start:days_end      
        lw(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    end
    fclose(fid);

    sw = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
    fid = fopen(fullfile(inputpath,swdownfile),'r','b');
    for k=days_start:days_end
        sw(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    end
    fclose(fid);

    precip = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
    fid = fopen(fullfile(inputpath,precipfile),'r','b');
    for k=days_start:days_end
        precip(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    end
    fclose(fid);

    apressure = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
    fid = fopen(fullfile(inputpath,apressurefile),'r','b');
    for k=days_start:days_end
        apressure(:,:,k-days_start+1) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    end
    fclose(fid);

    aq = zeros(EXF_Nx,EXF_Ny,length(days_start:days_end));
    fid = fopen(fullfile(inputpath,aqhfile),'r','b');
    for k=1:days_start:days_end
        aq(:,:,k) = fread(fid,[EXF_Nx EXF_Ny],'real*8');
    end
    fclose(fid); 
    
end


    
if tides ==1    

    fid = fopen(fullfile(inputpath,'OBEamFile'),'r','b');
    OBE = fread(fid,[Ny 10],'real*8');
    fclose(fid);
end


rhoShelf = 917;




XC = rdmds(fullfile(resultspath,'XC'));
XG = rdmds(fullfile(resultspath,'XG'));
YC = rdmds(fullfile(resultspath,'YC'));
YG = rdmds(fullfile(resultspath,'YG'));
DXC = rdmds(fullfile(resultspath,'DXC'));
DXG = rdmds(fullfile(resultspath,'DXG'));
DYC = rdmds(fullfile(resultspath,'DYC'));
DYG = rdmds(fullfile(resultspath,'DYG'));
RC = rdmds(fullfile(resultspath,'RC'));
RF = rdmds(fullfile(resultspath,'RF'));
RAC = rdmds(fullfile(resultspath,'RAC'));
RAS = rdmds(fullfile(resultspath,'RAS'));
RAW = rdmds(fullfile(resultspath,'RAW'));
RAZ = rdmds(fullfile(resultspath,'RAZ'));
DRF = rdmds(fullfile(resultspath,'DRF'));
DRC = rdmds(fullfile(resultspath,'DRC'));
hFacS = rdmds(fullfile(resultspath,'hFacS'));
hFacW = rdmds(fullfile(resultspath,'hFacW'));
hFacC = rdmds(fullfile(resultspath,'hFacC'));



zz = RC;
xx = XC(:,1);
yy = YC(1,:);

% run ../newexp/defineGrid;
% xx = XMC(:,1);
% yy = YMC(1,:);

% 

% 
% % 
% % 
% % 
