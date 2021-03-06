%%%
%%% setParams.m
%%%
%%% Sets basic MITgcm parameters plus parameters for included packages, and
%%% writes out the appropriate input files.,
%%%


function nTimeSteps = setParams (inputpath,codepath,listterm)  
  

  %%%%%%%%%%%%%%%%%%
  %%%%% SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%
  
  %%% Some commonly varied parameters
  
  %   nIter0 = 0; %%% Initial iteration for brand new runs
  nIter0 = 1; %%% Initial iteration for pickup runs
%   nIter0 = 1183320;
%   nIter0 = 1774980;
  use_extended_diags = true;  
  use_layers = true;
  use_tides = true;
  
  
  
  %%% GSW scripts
  addpath gsw_matlab_v3_05_4/ 
  addpath gsw_matlab_v3_05_4/library/
  addpath gsw_matlab_v3_05_4/thermodynamics_from_t/
  
  %%% If set true, plots of prescribed variables will be shown
  showplots = true;      
  fignum = 1;
  
  %%% Get parameter type definitions
  paramTypes;
  
  %%% To store parameter names and values
  parm01 = parmlist;
  parm02 = parmlist;
  parm03 = parmlist;
  parm04 = parmlist;
  parm05 = parmlist;
  PARM={parm01,parm02,parm03,parm04,parm05}; 

  %%% Model grid is defined externally
  run defineGrid.m   
  
 
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% FIXED PARAMETER VALUES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 


  %%% Diffusion parameters
  viscAh = 0; %%% Horizontal viscosity    
  viscA4 = 0; %%% Biharmonic viscosity
  viscAhGrid = 0; %%% Grid-dependent viscosity
  viscA4Grid = 0; %%% Grid-dependent biharmonic viscosity
  viscC4smag = 4; %%% Smagorninsky hyperviscosity parameter
  %%%%viscC4Leith=2.15;
  %%%%viscC4Leithd=2.15; %%%%% Modified Leith non-dimensional viscosity factor 
  viscAr = 3e-4; %%% Vertical viscosity
  diffKhT = 0; %%% Horizontal temp diffusion
  diffKrT = 0; %%% Vertical temp diffusion   
  
  %%% PARM01
  %%% momentum scheme
  %%% viscosity  
  
  parm01.addParm('implicSurfPress',0.6,PARM_REAL);
  parm01.addParm('implicDiv2DFlow',0.6,PARM_REAL);
  parm01.addParm('viscAr',viscAr,PARM_REAL);
  parm01.addParm('viscA4',viscA4,PARM_REAL);
  parm01.addParm('viscAh',viscAh,PARM_REAL);
  parm01.addParm('viscA4Grid',viscA4Grid,PARM_REAL);
  parm01.addParm('viscAhGrid',viscAhGrid,PARM_REAL);
  parm01.addParm('viscA4GridMax',0.5,PARM_REAL);
  parm01.addParm('viscAhGridMax',1,PARM_REAL);
  parm01.addParm('useAreaViscLength',false,PARM_BOOL);
  parm01.addParm('useFullLeith',true,PARM_BOOL);
  parm01.addParm('viscC4leith',0,PARM_REAL);
  parm01.addParm('viscC4leithD',0,PARM_REAL);  
  parm01.addParm('viscC4smag',viscC4smag,PARM_REAL);  
%   parm01.addParm('viscC4leith',2.15,PARM_REAL);
%   parm01.addParm('viscC4leithD',2.15,PARM_REAL);
  parm01.addParm('viscC2leith',0,PARM_REAL);
  parm01.addParm('viscC2leithD',0,PARM_REAL);  
  %%% diffusivity
  parm01.addParm('diffKrT',diffKrT,PARM_REAL);
  parm01.addParm('diffKhT',diffKhT,PARM_REAL);
  parm01.addParm('diffK4T',0,PARM_REAL);  
  %%% advection schemes
  parm01.addParm('tempAdvScheme',7,PARM_INT);
  parm01.addParm('saltAdvScheme',7,PARM_INT);
  parm01.addParm('multiDimAdvection',true,PARM_BOOL);
  parm01.addParm('tempStepping',true,PARM_BOOL);
  parm01.addParm('saltStepping',true,PARM_BOOL);
  parm01.addParm('staggerTimeStep',true,PARM_BOOL);
%   parm01.addParm('convertFW2Salt',-1,PARM_REAL);
  %%% equation of state
  parm01.addParm('eosType','JMD95Z',PARM_STR);   
  %%% boundary conditions
  parm01.addParm('no_slip_sides',false,PARM_BOOL);
%   parm01.addParm('no_slip_sides',true,PARM_BOOL);
  parm01.addParm('no_slip_bottom',false,PARM_BOOL);
%   parm01.addParm('no_slip_bottom',true,PARM_BOOL);
  parm01.addParm('bottomDragLinear',0,PARM_REAL);
  parm01.addParm('bottomDragQuadratic',2.1e-3,PARM_REAL);
  %%% physical parameters
  parm01.addParm('gravity',g,PARM_REAL);
  parm01.addParm('rhonil',rho0,PARM_REAL);
  parm01.addParm('rhoConst',rho0,PARM_REAL);
 
  %%% implicit diffusion and convective adjustment  
  parm01.addParm('ivdc_kappa',0,PARM_REAL);
  parm01.addParm('implicitDiffusion',true,PARM_BOOL);
  parm01.addParm('implicitViscosity',true,PARM_BOOL);
  %%% exact volume cxsonservation
  parm01.addParm('exactConserv',true,PARM_BOOL);
  %%% C-V scheme for Coriolis term
  parm01.addParm('useCDscheme',false,PARM_BOOL);
  %%% partial cells for smooth topography
  parm01.addParm('hFacMin',hFacMin,PARM_REAL); 
  parm01.addParm('hFacMinDr',hFacMinDr,PARM_REAL);

  %%% file IO stuff
  parm01.addParm('readBinaryPrec',64,PARM_INT);
  parm01.addParm('useSingleCpuIO',true,PARM_BOOL);
%   parm01.addParm('debugLevel',-1,PARM_INT);
  parm01.addParm('debugLevel',2,PARM_INT);

%%%%% Vertical advection
  parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL);
%   parm01.addParm('momImplVertAdv',true,PARM_BOOL);
%   parm01.addParm('tempImplVertAdv',true,PARM_BOOL);
%   parm01.addParm('saltImplVertAdv',true,PARM_BOOL);


  %%% Wet-point method at boundaries - may improve boundary stability
  parm01.addParm('useJamartWetPoints',true,PARM_BOOL);
  parm01.addParm('useJamartMomAdv',true,PARM_BOOL);

  %%% PARM02
  parm02.addParm('useSRCGSolver',true,PARM_BOOL);  
  parm02.addParm('cg2dMaxIters',1000,PARM_INT);  
  parm02.addParm('cg2dTargetResidual',1e-12,PARM_REAL);
 
  %%% PARM03
  %parm03.addParm('alph_AB',1/2,PARM_REAL);
  %parm03.addParm('beta_AB',5/12,PARM_REAL);
  parm03.addParm('momDissip_In_AB',false,PARM_BOOL);
  
  %%%%%%% Added to make SEAICE work
  %parm03.addParm('momForcingOutAB',1,PARM_INT);
  parm03.addParm('tracForcingOutAB',1,PARM_INT);


  %%%%%%%
  parm03.addParm('nIter0',nIter0,PARM_INT);
  parm03.addParm('abEps',0.1,PARM_REAL);
%   parm03.addParm('chkptFreq',0.01*t1year,PARM_REAL);
%   parm03.addParm('chkptFreq',0.1*t1year,PARM_REAL);
%   parm03.addParm('pChkptFreq',1*t1year,PARM_REAL);
  parm03.addParm('chkptFreq',t1month,PARM_REAL);
  parm03.addParm('pChkptFreq',t1year,PARM_REAL);
%   parm03.addParm('pChkptFreq',t1month,PARM_REAL);
  parm03.addParm('taveFreq',0,PARM_REAL);
  parm03.addParm('dumpFreq',0,PARM_REAL);
  parm03.addParm('monitorFreq',t1year,PARM_REAL);
  parm03.addParm('cAdjFreq',0,PARM_REAL);
  parm03.addParm('dumpInitAndLast',true,PARM_BOOL);
  parm03.addParm('pickupStrictlyMatch',false,PARM_BOOL);

  %%% PARM04
  parm04.addParm('usingCartesianGrid',false,PARM_BOOL);
  parm04.addParm('usingSphericalPolarGrid',true,PARM_BOOL);    
  
  %%% Store grid spacings
  parm04.addParm('delX',dmxg,PARM_REALS);
  parm04.addParm('delY',dmyg,PARM_REALS);

  
  parm04.addParm('delR',dz,PARM_REALS);      
  parm04.addParm('ygOrigin',ymin,PARM_REAL);
  parm04.addParm('xgOrigin',xmin,PARM_REAL);
  
 
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% BATHYMETRY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
   
    
  parm05.addParm('bathyFile',bathyFile,PARM_STR); 
  
  fid = fopen(fullfile(inputpath,bathyFile),'r','b');
  h = fread(fid,[Nx Ny],'real*8'); 
  fclose(fid);
  

   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Make RIGID WALL WITH TOPOGRAPHY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  fid = fopen(fullfile(inputpath,SHELFICEtopoFile),'r','b');
  SHELFICEtopo = fread(fid,[Nx Ny],'real*8'); 
  fclose(fid);
                
  %%% Ensure ocean depth is zero at western and southern boundaries
  h(:,1) = 0;
  h(1,:) = 0;
  SHELFICEtopo(:,1) = 0;
  SHELFICEtopo(1,:) = 0;
  
  %%% Eliminate any spurious openings at the northern boundary
  idx_obcs_n = find(h(:,end)>=0,1,'last');
  h(1:idx_obcs_n,end) = 0;
  idx_obcs_e = find(h(end,:)>=0,1,'last');
 
  %%% Overwrite bathymetry data file
  writeDataset(h,fullfile(inputpath,bathyFile),ieee,prec);
  clear h

  
  %%% Overwrite ice draft data file  
  writeDataset(SHELFICEtopo,fullfile(inputpath,SHELFICEtopoFile),ieee,prec);
  clear SHELFICEtopo
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL DATA %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%

     
  %%% Initial conditions set externally - just tell MITgcm where the
  %%% files are
  parm05.addParm('hydrogThetaFile',hydrogThetaFile,PARM_STR);
  parm05.addParm('hydrogSaltFile',hydrogSaltFile,PARM_STR);
  parm05.addParm('pSurfInitFile',pSurfInitFile,PARM_STR);
  parm05.addParm('uVelInitFile',uVelInitFile,PARM_STR);
  parm05.addParm('vVelInitFile',vVelInitFile,PARM_STR);
      
  
  
  %%%%%%%%%%%%%%%%%%%%%
  %%%%% TIME STEP %%%%%
  %%%%%%%%%%%%%%%%%%%%%
  
  %%% Overall minimum time step
  fid = fopen('TIME_STEP','r');
  deltaT = round(fscanf(fid,'%f'));  
  fclose(fid);
  nTimeSteps = ceil(simTime/deltaT);
  simTimeAct = nTimeSteps*deltaT;
  
  %%% Write end time time step size  
  parm03.addParm('endTime',nIter0*deltaT+simTimeAct,PARM_INT);
  parm03.addParm('deltaT',deltaT,PARM_REAL); 
    
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SURFACE WIND FORCING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  %%% Save as a parameter  
%   writeDataset(tau_mat,fullfile(inputpath,'zonalWindFile.bin'),ieee,prec); 
%   parm05.addParm('zonalWindFile','zonalWindFile.bin',PARM_STR);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SURFACE HEAT/SALT FLUXES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% Save as parameters
%   writeDataset(heat_flux,fullfile(inputpath,'surfQfile.bin'),ieee,prec);
%   parm05.addParm('surfQfile','surfQfile.bin',PARM_STR);  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data' file
  write_data(inputpath,PARM,listterm,realfmt);
  
 
  
  
  
  
  
  
  
  
  
  
  
  
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%%%%%%% SEA ICE $ %%%%%%%%%%%
%   
%   
%   % to store parameter names and values
  seaice_parm01 = parmlist;
  SEAICE_PARM = {seaice_parm01};
% %   
% %   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% SEA ICE  %%%%%%%%%%%%
    %%%%%%%% PARAMETERS %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oringinal albedos from llc_1/48th  other values from llc_2160 or 
%%% Revised albedos from SOSE
% 
  SEAICEwriteState   = true;
  SEAICEuseDYNAMICS  = true;
  SEAICE_multDim     = 7;
%   SEAICE_dryIceAlb   = 0.8783;
%   SEAICE_dryIceAlb   = 0.8509;
  SEAICE_dryIceAlb   = 0.92;
%   SEAICE_wetIceAlb   = 0.7869;
%   SEAICE_wetIceAlb   = 0.7284;
  SEAICE_wetIceAlb   = 0.80;
%   SEAICE_drySnowAlb  = 0.9482;
%   SEAICE_drySnowAlb  = 0.7754;
  SEAICE_drySnowAlb  = 0.96;
%   SEAICE_wetSnowAlb  = 0.8216;
%   SEAICE_wetSnowAlb  = 0.7753;
  SEAICE_wetSnowAlb  = 0.83;
  SEAICE_waterDrag   = 5.5399;
  SEAICE_drag        = 0.002;
%   HO                 = 0.1;
  HO                 = .5;
%%%% test .1,.5

  SEAICE_no_slip          = false;
%   SEAICE_no_slip          = true;

  SEAICEadvScheme         = 7;
%   SEAICEadvScheme         = 33;


  %%%SOSEdoesn't have a seaice dataset for salinity, they used this value
  %%%in their estimate
  
  LSR_ERROR               = 1.0e-4;
%   LSR_ERROR               = 2.0e-4; 
  MIN_ATEMP               = -40;
  MIN_TICE                = -40;
  SEAICE_area_reg         = 0.15;
  SEAICE_hice_reg         = 0.1;
  IMAX_TICE               = 10;
  SEAICE_EPS		      = 1.0e-8;
%   SEAICE_EPS              = 2.0e-9;
  SEAICE_doOpenWaterMelt  = true;
  SEAICE_areaLossFormula  = 1;
  SEAICE_wetAlbTemp       = 0.0;
  SEAICE_saltFrac         = 0.0;
  SEAICE_frazilFrac       = 0.003;
%  SEAICE_frazilFrac       = 0.01;
%   SEAICE_frazilFrac       = 1.0;
  SEAICEscaleSurfStress   = true;
  
  seaice_parm01.addParm('LSR_ERROR',LSR_ERROR,PARM_REAL);
  seaice_parm01.addParm('SEAICEwriteState',SEAICEwriteState,PARM_BOOL);
  seaice_parm01.addParm('SEAICEuseDYNAMICS',SEAICEuseDYNAMICS,PARM_BOOL);
  seaice_parm01.addParm('SEAICE_multDim',SEAICE_multDim,PARM_INT);
  seaice_parm01.addParm('SEAICE_dryIceAlb',SEAICE_dryIceAlb,PARM_REAL);
  seaice_parm01.addParm('SEAICE_wetIceAlb',SEAICE_wetIceAlb,PARM_REAL);
  seaice_parm01.addParm('SEAICE_drySnowAlb',SEAICE_drySnowAlb,PARM_REAL);
  seaice_parm01.addParm('SEAICE_wetSnowAlb',SEAICE_wetSnowAlb,PARM_REAL);
  seaice_parm01.addParm('SEAICE_waterDrag',SEAICE_waterDrag,PARM_REAL);
  seaice_parm01.addParm('SEAICE_drag',SEAICE_drag,PARM_REAL);
  seaice_parm01.addParm('HO',HO,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_dryIceAlb_south',SEAICE_dryIceAlb_south,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_wetIceAlb_south',SEAICE_wetIceAlb_south,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_drySnowAlb_south',SEAICE_drySnowAlb_south,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_wetSnowAlb_south',SEAICE_wetSnowAlb_south,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_waterDrag_south',SEAICE_waterDrag_south,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_drag_south',SEAICE_drag_south,PARM_REAL);
  seaice_parm01.addParm('SEAICE_no_slip',SEAICE_no_slip,PARM_BOOL);
%   seaice_parm01.addParm('SEAICE_salinity',SEAICE_salinity,PARM_REAL);
  seaice_parm01.addParm('SEAICEadvScheme',SEAICEadvScheme,PARM_INT);
  seaice_parm01.addParm('MIN_ATEMP',MIN_ATEMP,PARM_REAL);
  seaice_parm01.addParm('MIN_TICE',MIN_TICE,PARM_REAL);
  seaice_parm01.addParm('SEAICE_area_reg',SEAICE_area_reg,PARM_REAL);
  seaice_parm01.addParm('SEAICE_hice_reg',SEAICE_hice_reg,PARM_REAL);
  seaice_parm01.addParm('IMAX_TICE',IMAX_TICE,PARM_INT);
  seaice_parm01.addParm('SEAICE_EPS',SEAICE_EPS,PARM_REAL);
  seaice_parm01.addParm('SEAICE_doOpenWaterMelt',SEAICE_doOpenWaterMelt,PARM_BOOL);
  seaice_parm01.addParm('SEAICE_areaLossFormula',SEAICE_areaLossFormula,PARM_INT);
  seaice_parm01.addParm('SEAICE_wetAlbTemp',SEAICE_wetAlbTemp,PARM_REAL);
  seaice_parm01.addParm('SEAICE_saltFrac',SEAICE_saltFrac,PARM_REAL);
  seaice_parm01.addParm('SEAICE_frazilFrac',SEAICE_frazilFrac,PARM_REAL);
  seaice_parm01.addParm('SEAICEscaleSurfStress',SEAICEscaleSurfStress,PARM_BOOL);


  
  seaice_parm01.addParm('HeffFile',HeffFile,PARM_STR);
  seaice_parm01.addParm('AreaFile',AreaFile,PARM_STR);
  seaice_parm01.addParm('HsnowFile',HsnowFile,PARM_STR);
  seaice_parm01.addParm('HsaltFile',HsaltFile,PARM_STR);
  seaice_parm01.addParm('uIceFile',uIceFile,PARM_STR);
  seaice_parm01.addParm('vIceFile',vIceFile,PARM_STR);
  
  
     
%   
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%%% WRITE THE 'data.seaice' FILE %%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   
  write_data_seaice(inputpath,SEAICE_PARM,listterm,realfmt);  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%% SHELF ICE  %%%%%%%%%%
  %%%%     PARAMETERS       %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  % to store parameter names and values
  shelfice_parm01 = parmlist;
  SHELFICE_PARM = {shelfice_parm01};
 
 
  SHELFICEloadAnomalyFile = 'SHELFICEloadAnomalyFile.bin';
  SHELFICEtopoFile = 'SHELFICEtopoFile.bin';
  SHELFICEuseGammaFrict = true;
  SHELFICEboundaryLayer = false;
  SHELFICEconserve = false;
%   SHELFICEheatTransCoeff = .0005;
%   SHELFICEheatTransCoeff = .0001;
  SHELFICEheatTransCoeff = 0;

  
  
  shelfice_parm01.addParm('SHELFICEloadAnomalyFile',SHELFICEloadAnomalyFile,PARM_STR);
  shelfice_parm01.addParm('SHELFICEtopoFile',SHELFICEtopoFile,PARM_STR);
  shelfice_parm01.addParm('SHELFICEuseGammaFrict',SHELFICEuseGammaFrict,PARM_BOOL);
  shelfice_parm01.addParm('SHELFICEboundaryLayer',SHELFICEboundaryLayer,PARM_BOOL);
  shelfice_parm01.addParm('SHELFICEconserve',SHELFICEconserve,PARM_BOOL);
  shelfice_parm01.addParm('SHELFICEheatTransCoeff',SHELFICEheatTransCoeff,PARM_REAL);


  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.shelfice' FILE %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  write_data_shelfice(inputpath,SHELFICE_PARM,listterm,realfmt);  

  
  
  
  

  
  
    
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% OBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  obcs_parm01 = parmlist;
  obcs_parm02 = parmlist;
  obcs_parm03 = parmlist;
  obcs_parm04 = parmlist;
  obcs_parm05 = parmlist;
  OBCS_PARM = {obcs_parm01,obcs_parm02,obcs_parm03,obcs_parm04,obcs_parm05};  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFINE OPEN BOUNDARY TYPES (OBCS_PARM01) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%%%%% tides ^%%%%%%%% 
  useOBCStides = use_tides;

  tidalPeriod=[44714, 43200, 45570, 43082, 86164, 92950, 86637, 96726,1180300,2380706];
%   tidalPeriod=[44714, 43200, 45570, 43082, 86164, 92950, 86637];

  OBNamFile= 'OBNamFile.bin';
  OBNphFile= 'OBNphFile.bin';
  OBEamFile= 'OBEamFile.bin';
  OBEphFile= 'OBEphFile.bin';
  
  
  
  obcs_parm01.addParm('useOBCStides',useOBCStides,PARM_BOOL);
  
  obcs_parm01.addParm('tidalPeriod',tidalPeriod,PARM_INTS);    

  obcs_parm01.addParm('OBNamFile',OBNamFile,PARM_STR);  
  obcs_parm01.addParm('OBNphFile',OBNphFile,PARM_STR); 
  obcs_parm01.addParm('OBEamFile',OBEamFile,PARM_STR); 
  obcs_parm01.addParm('OBEphFile',OBEphFile,PARM_STR); 

  
  
  
  %%% Enables an Orlanski radiation condition at the northern boundary
  
  
    useOrlanskiNorth = false;

    OB_Ieast = -1*ones(1,Ny);
%     OB_Iwest = 1*ones(1,Ny);
    OB_Jnorth = -1*ones(1,Nx);
  OB_Jnorth(1:idx_obcs_n) = 0;

     
  obcs_parm01.addParm('useOrlanskiNorth',useOrlanskiNorth,PARM_BOOL);
  obcs_parm01.addParm('OB_Jnorth',OB_Jnorth,PARM_INTS);    
  obcs_parm01.addParm('OB_Ieast',OB_Ieast,PARM_INTS);    
%   obcs_parm01.addParm('OB_Iwest',OB_Iwest,PARM_INTS);    






    useOBCSsponge = true;
    useSeaiceSponge = true;
    
%     useOBCSsponge = false;

    useOBCSprescribe = true;
    
    
    OBNtFile = 'OBNtFile.bin';
    OBEtFile = 'OBEtFile.bin';
%     OBWtFile = 'OBWtFile.bin';

    OBNsFile = 'OBNsFile.bin';
    OBEsFile = 'OBEsFile.bin';
%     OBWsFile = 'OBWsFile.bin';
    
    OBNuFile = 'OBNuFile.bin';
    OBEuFile = 'OBEuFile.bin';
%     OBWuFile = 'OBWuFile.bin';

    OBNvFile = 'OBNvFile.bin';
    OBEvFile = 'OBEvFile.bin';
%     OBWvFile = 'OBWvFile.bin';

    OBNetaFile = 'OBNetaFile.bin';
    OBEetaFile = 'OBEetaFile.bin';
%     OBWetaFile = 'OBWetaFile.bin';

    OBNaFile = 'OBNaFile.bin';
    OBEaFile = 'OBEaFile.bin';
%     OBWaFile = 'OBWaFile.bin';

  fid = fopen(fullfile(inputpath,OBNaFile),'r','b');
  OBNa = fread(fid,[Nx 12],'real*8'); 
  fclose(fid);


    OBNhFile = 'OBNhFile.bin';
    OBEhFile = 'OBEhFile.bin';
%     OBWhFile = 'OBWhFile.bin';

    OBNsnFile = 'OBNsnFile.bin';
    OBEsnFile = 'OBEsnFile.bin';
%     OBWsnFile = 'OBWsnFile.bin';

    OBNuiceFile = 'OBNuiceFile.bin';
    OBEuiceFile = 'OBEuiceFile.bin';
% %     OBWuiceFile = 'OBWuiceFile.bin';
% 
    OBNviceFile = 'OBNviceFile.bin';
    OBEviceFile = 'OBEviceFile.bin';
%     OBWviceFile = 'OBWviceFile.bin';


  obcs_parm01.addParm('useOBCSsponge',useOBCSsponge,PARM_BOOL);
  obcs_parm01.addParm('useSeaiceSponge',useSeaiceSponge,PARM_BOOL);
  obcs_parm01.addParm('useOBCSprescribe',useOBCSprescribe,PARM_BOOL);  

  obcs_parm01.addParm('OBNtFile',OBNtFile,PARM_STR);  
  obcs_parm01.addParm('OBEtFile',OBEtFile,PARM_STR);  
%   obcs_parm01.addParm('OBWtFile',OBWtFile,PARM_STR);  

  obcs_parm01.addParm('OBNsFile',OBNsFile,PARM_STR);  
  obcs_parm01.addParm('OBEsFile',OBEsFile,PARM_STR);  
%   obcs_parm01.addParm('OBWsFile',OBWsFile,PARM_STR);

  obcs_parm01.addParm('OBNuFile',OBNuFile,PARM_STR);  
  obcs_parm01.addParm('OBEuFile',OBEuFile,PARM_STR);  
%   obcs_parm01.addParm('OBWuFile',OBWuFile,PARM_STR);

  obcs_parm01.addParm('OBNvFile',OBNvFile,PARM_STR);  
  obcs_parm01.addParm('OBEvFile',OBEvFile,PARM_STR);  
%   obcs_parm01.addParm('OBWvFile',OBWvFile,PARM_STR);

%   obcs_parm01.addParm('OBNetaFile',OBNetaFile,PARM_STR);  
%   obcs_parm01.addParm('OBEetaFile',OBEetaFile,PARM_STR);  
%   obcs_parm01.addParm('OBWetaFile',OBWetaFile,PARM_STR);

  obcs_parm01.addParm('OBNaFile',OBNaFile,PARM_STR);  
  obcs_parm01.addParm('OBEaFile',OBEaFile,PARM_STR);  
%   obcs_parm01.addParm('OBWaFile',OBWaFile,PARM_STR);

  obcs_parm01.addParm('OBNhFile',OBNhFile,PARM_STR);  
  obcs_parm01.addParm('OBEhFile',OBEhFile,PARM_STR);  
%   obcs_parm01.addParm('OBWhFile',OBWhFile,PARM_STR);

  obcs_parm01.addParm('OBNsnFile',OBNsnFile,PARM_STR);  
  obcs_parm01.addParm('OBEsnFile',OBEsnFile,PARM_STR);  
%   obcs_parm01.addParm('OBWsnFile',OBWsnFile,PARM_STR);
% 
  obcs_parm01.addParm('OBNuiceFile',OBNuiceFile,PARM_STR);  
  obcs_parm01.addParm('OBEuiceFile',OBEuiceFile,PARM_STR);  
%   obcs_parm01.addParm('OBWuiceFile',OBWuiceFile,PARM_STR);

  obcs_parm01.addParm('OBNviceFile',OBNviceFile,PARM_STR);  
  obcs_parm01.addParm('OBEviceFile',OBEviceFile,PARM_STR);  
%   obcs_parm01.addParm('OBWviceFile',OBWviceFile,PARM_STR);


   
  %%% Enforces mass conservation across the northern boundary by adding a
  %%% barotropic inflow/outflow
  useOBCSbalance = true;
  OBCS_balanceFacN = 1; 
  OBCS_balanceFacE = 0;
%   OBCS_balanceFacS = 0;
%   OBCS_balanceFacW = -1;
  obcs_parm01.addParm('useOBCSbalance',useOBCSbalance,PARM_BOOL);  
  obcs_parm01.addParm('OBCS_balanceFacN',OBCS_balanceFacN,PARM_REAL);  
  obcs_parm01.addParm('OBCS_balanceFacE',OBCS_balanceFacE,PARM_REAL);  
%   obcs_parm01.addParm('OBCS_balanceFacS',OBCS_balanceFacS,PARM_REAL);  
%   obcs_parm01.addParm('OBCS_balanceFacW',OBCS_balanceFacW,PARM_REAL);  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ORLANSKI OPTIONS (OBCS_PARM02) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  





  
  %%% Velocity averaging time scale - must be larger than deltaT.
  %%% The Orlanski radiation condition computes the characteristic velocity
  %%% at the boundary by averaging the spatial derivative normal to the 
  %%% boundary divided by the time step over this period.
  %%% At the moment we're using the magic engineering factor of 3.
%   cvelTimeScale = 3*deltaT;




  %%% Max dimensionless CFL for Adams-Basthforth 2nd-order method
%   CMAX = 0.45; 
%   
%   obcs_parm02.addParm('cvelTimeScale',cvelTimeScale,PARM_REAL);
%   obcs_parm02.addParm('CMAX',CMAX,PARM_REAL);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Sponge Layer Parms (OBCS_PARM03) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Values taken from SOSE.
%%% Urelaxobcsinner = relaxation time scale at the innermost sponge layer point of a meridional OB
%%% Vrelaxobcsinner = relaxation time scale at the innermost sponge layer point of a zonal OB
%%% Urelaxobcsbound = relaxation time scale at the outermost sponge layer point of a meridional OB
%%% Vrelaxobcsbound = relaxation time scale at the outermost sponge layer point of a zonal OB



    Urelaxobcsinner = 864000;  %%% 10 days
    Urelaxobcsbound = 43200;  %%% half a day
    Vrelaxobcsinner = 864000;
    Vrelaxobcsbound = 43200;
    
%%%%%% sponge thickness - finer in high-res simulation due to placement of
%%%%%% eastern boundary, increased number of gridpoints
  if (res_fac == 24)
    spongethickness = 48;
  else
    spongethickness = round(10*res_fac/3);
%     spongethickness = 5;
  end

 
  obcs_parm03.addParm('Urelaxobcsinner',Urelaxobcsinner,PARM_REAL);
  obcs_parm03.addParm('Urelaxobcsbound',Urelaxobcsbound,PARM_REAL);
  obcs_parm03.addParm('Vrelaxobcsinner',Vrelaxobcsinner,PARM_REAL);
  obcs_parm03.addParm('Vrelaxobcsbound',Vrelaxobcsbound,PARM_REAL);

  obcs_parm03.addParm('spongethickness',spongethickness,PARM_INT);
  


  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Sea ice Sponge Parms (OBCS_PARM05) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  seaiceSpongeThickness = spongethickness;
  Arelaxobcsinner = Urelaxobcsinner;
  Arelaxobcsbound = Urelaxobcsbound;
  Hrelaxobcsinner = Urelaxobcsinner;
  Hrelaxobcsbound = Urelaxobcsbound;
  SLrelaxobcsinner = Urelaxobcsinner;
  SLrelaxobcsbound = Urelaxobcsbound;
  SNrelaxobcsinner = Urelaxobcsinner;
  SNrelaxobcsbound = Urelaxobcsbound;  
  
  obcs_parm05.addParm('Arelaxobcsinner',Arelaxobcsinner,PARM_REAL);
  obcs_parm05.addParm('Arelaxobcsbound',Arelaxobcsbound,PARM_REAL); 
  obcs_parm05.addParm('Hrelaxobcsinner',Hrelaxobcsinner,PARM_REAL);
  obcs_parm05.addParm('Hrelaxobcsbound',Hrelaxobcsbound,PARM_REAL); 
  obcs_parm05.addParm('SLrelaxobcsinner',SLrelaxobcsinner,PARM_REAL);
  obcs_parm05.addParm('SLrelaxobcsbound',SLrelaxobcsbound,PARM_REAL); 
  obcs_parm05.addParm('SNrelaxobcsinner',SNrelaxobcsinner,PARM_REAL);
  obcs_parm05.addParm('SNrelaxobcsbound',SNrelaxobcsbound,PARM_REAL); 
  obcs_parm05.addParm('seaiceSpongeThickness',seaiceSpongeThickness,PARM_INT);

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.obcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data.obcs' file
  write_data_obcs(inputpath,OBCS_PARM,listterm,realfmt);
  
  
  
  
  
  
  
  
  
    
%   %%%%%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%
%   %%%%% GMREDI %%%%%
%   %%%%%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%
%   
%   
%   
%   
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%%% GMREDI SET-UP %%%%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   
%   %%% To store parameter names and values
%   gmredi_parm01 = parmlist;
%   GMREDI_PARM = {gmredi_parm01};
%   
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%%% GMREDI PARAMETERS %%%%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   
%   %%% Define parameters for gmredi package %%%
% 
%   %%% Isopycnal diffusivity
%   GM_isopycK          = 100;
% 
%   %%% Thickness diffusivity
%   GM_background_K     = 100;
% 
%   %%% Maximum isopycnal slope (I think this is only used by certain
%   %%% tapering schemes)
%   GM_maxSlope         = 0.025;
%   
%   %%% Tapering scheme
%   GM_taper_scheme     = 'dm95';
% 
%   %%% DM95 critical slope
%   GM_Scrit            = 0.025;
% 
%   %%% DM95 tapering width
%   GM_Sd               = 0.0025;
% 
%   %%% Add parameters
%   gmredi_parm01.addParm('GM_isopycK',GM_isopycK,PARM_REAL);
%   gmredi_parm01.addParm('GM_background_K',GM_background_K,PARM_REAL);
%   gmredi_parm01.addParm('GM_maxSlope',GM_maxSlope,PARM_REAL);
%   gmredi_parm01.addParm('GM_taper_scheme',GM_taper_scheme,PARM_STR);
%   gmredi_parm01.addParm('GM_Scrit',GM_Scrit,PARM_REAL);
%   gmredi_parm01.addParm('GM_Sd',GM_Sd,PARM_REAL);
%   
%   %%z% Create the data.gmredi file
%   write_data_gmredi(inputpath,GMREDI_PARM,listterm,realfmt);
  
















    %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%
    %%%%%%EXF PKG%%%%%
    %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%
    
    
  %%% To store parameter names and values
    
     EXF_NML_01 = parmlist;
     EXF_NML_02 = parmlist;
     EXF_NML_03 = parmlist;
     EXF_NML_04 = parmlist;
     EXF_NML_OBCS = parmlist;

     EXF_PARM = {EXF_NML_01,EXF_NML_02,EXF_NML_03,EXF_NML_04,EXF_NML_OBCS};  
    
    


% %     EXF_NML_01

  	exf_albedo        = 0.15;
 	exf_scal_BulkCdn  = 1.015;
 	exf_iprec         = 64;  
 	useExfYearlyFields= false;
 	useExfCheckRange  = false;
 	useRelativeWind   = true;
%  	useRelativeWind   = false;
    repeatPeriod      = 283996800;
%     repeatPeriod = 31536000;
    exf_offset_atemp =  273.16;
    
    
%%%runoff from ERA is in hours, need to convert to seconds
%     exf_inscal_runoff = 1.14e-04;
    
    
    apressurefile     = 'apressurefile.bin';
    atempfile         = 'atempfile.bin';
    aqhfile           = 'aqhfile.bin';
    uwindfile         = 'uwindfile.bin';
    vwindfile         = 'vwindfile.bin';
    precipfile        = 'precipfile.bin';
    swdownfile        = 'swdownfile.bin';
    lwdownfile        = 'lwdownfile.bin';
%     runofffile        = 'runofffile.bin';
    
    
% "*period=-12" specifies monthly-mean forcing
   apressurestartdate1 = str2num([num2str(start_year),'0101']);
   apressurestartdate2 = 000000;
   apressureperiod     = 86400.0;
%    apressureperiod     = 10800.0;
    
    aqhstartdate1 = str2num([num2str(start_year),'0101']);
    aqhstartdate2 = 000000;
    aqhperiod = 86400.0;
%     aqhperiod           = 10800.0;

    atempstartdate1 = str2num([num2str(start_year),'0101']);
    atempstartdate2 = 000000;
    atempperiod = 86400.0;
%     atempperiod         = 10800.0;
 

    uwindstartdate1 = str2num([num2str(start_year),'0101']);
    uwindstartdate2 = 000000;
   uwindperiod = 86400.0;
%     uwindperiod   = 10800.0;
 
    vwindstartdate1 = str2num([num2str(start_year),'0101']);
    vwindstartdate2 = 000000;
    vwindperiod = 86400.0;
%     vwindperiod         = 10800.0; 
 
    precipstartdate1 = str2num([num2str(start_year),'0101']);
    precipstartdate2 = 000000;
    precipperiod = 86400.0;
%     precipperiod        = 10800.0; 
 

    swdownstartdate1 = str2num([num2str(start_year),'0101']);
    swdownstartdate2 = 000000;
    swdownperiod = 86400.0;
%     swdownperiod        = 10800.0;
% 

    lwdownstartdate1 = str2num([num2str(start_year),'0101']);
    lwdownstartdate2 = 000000;
    lwdownperiod = 86400.0;
%     lwdownperiod        = 10800.0;


%     runoffstartdate1 = str2num([num2str(start_year),'0101']);
%     runoffstartdate2 = 000000;
%     runoffperiod = 2592000.0;
%    runoffperiod        = 10800.0;

   EXF_dmxg = EXF_dmxg*(sum(EXF_dmxg)/sum(EXF_dmxg(1:end-1)));
   EXF_dmyg = EXF_dmyg*(sum(EXF_dmyg)/sum(EXF_dmyg(1:end-1)));
   
 
   precip_lon0 = EXF_xmin;
   precip_lon_inc = EXF_dmxg(1);
   precip_lat0 = EXF_ymin;
   precip_lat_inc = EXF_dmyg(1:end-1);
   precip_nlon = EXF_Nx;
   precip_nlat = EXF_Ny;
   
   atemp_lon0 = EXF_xmin;
   atemp_lon_inc = EXF_dmxg(1);
   atemp_lat0 = EXF_ymin;
   atemp_lat_inc = EXF_dmyg(1:end-1);
   atemp_nlon = EXF_Nx;
   atemp_nlat = EXF_Ny;
   
   apressure_lon0 = EXF_xmin;
   apressure_lon_inc = EXF_dmxg(1);
   apressure_lat0 = EXF_ymin;
   apressure_lat_inc = EXF_dmyg(1:end-1);
   apressure_nlon = EXF_Nx;
   apressure_nlat = EXF_Ny;
    
   aqh_lon0 = EXF_xmin;
   aqh_lon_inc = EXF_dmxg(1);
   aqh_lat0 = EXF_ymin;
   aqh_lat_inc = EXF_dmyg(1:end-1);
   aqh_nlon = EXF_Nx;
   aqh_nlat = EXF_Ny;
   
   uwind_lon0 = EXF_xmin;
   uwind_lon_inc = EXF_dmxg(1);
   uwind_lat0 = EXF_ymin;
   uwind_lat_inc = EXF_dmyg(1:end-1);
   uwind_nlon = EXF_Nx;
   uwind_nlat = EXF_Ny;
   
   
   vwind_lon0 = EXF_xmin;
   vwind_lon_inc = EXF_dmxg(1);
   vwind_lat0 = EXF_ymin;
   vwind_lat_inc = EXF_dmyg(1:end-1);
   vwind_nlon = EXF_Nx;
   vwind_nlat = EXF_Ny;
   
   swdown_lon0 = EXF_xmin;
   swdown_lon_inc = EXF_dmxg(1);
   swdown_lat0 = EXF_ymin;
   swdown_lat_inc = EXF_dmyg(1:end-1);
   swdown_nlon = EXF_Nx;
   swdown_nlat = EXF_Ny;
   
   lwdown_lon0 = EXF_xmin;
   lwdown_lon_inc = EXF_dmxg(1);
   lwdown_lat0 = EXF_ymin;
   lwdown_lat_inc = EXF_dmyg(1:end-1);
   lwdown_nlon = EXF_Nx;
   lwdown_nlat = EXF_Ny;

%    runoff_lon0 = EXF_xmin;
%    runoff_lon_inc = dmxg(1);
%    runoff_lat0 = EXF_ymin;
%    runoff_lat_inc = dmyg;
%    runoff_nlon = Nx;
%    runoff_nlat = Ny;
  
 
  %%%%%%%%%%%%%%%%%%% EXF_NML_OBCS %%%%%%%%%%%%%%%%%%%%
  
  %%% High-res configuration uses monthly means from lower-res simulations
  if (res_fac == 24)
        
    obcsPeriod = round(t1month);
    obcsstartdate =  datenum([num2str(start_year),'-01-01'])-t1month/2/t1day;
    obcsstartdate1 = datestr(obcsstartdate,'yyyymmdd');
    obcsstartdate2 = datestr(obcsstartdate,'HHMMSS');

    EXF_NML_OBCS.addParm('obcsNstartdate1',obcsstartdate1,PARM_MISC);
    EXF_NML_OBCS.addParm('obcsNstartdate2',obcsstartdate2,PARM_MISC);
    
    EXF_NML_OBCS.addParm('obcsEstartdate1',obcsstartdate1,PARM_MISC);
    EXF_NML_OBCS.addParm('obcsEstartdate2',obcsstartdate2,PARM_MISC);
   
    EXF_NML_OBCS.addParm('siobEstartdate1',obcsstartdate1,PARM_MISC);
    EXF_NML_OBCS.addParm('siobEstartdate2',obcsstartdate2,PARM_MISC);
    
    EXF_NML_OBCS.addParm('siobNstartdate1',obcsstartdate1,PARM_MISC);
    EXF_NML_OBCS.addParm('siobNstartdate2',obcsstartdate2,PARM_MISC);

    EXF_NML_OBCS.addParm('obcsNperiod',obcsPeriod,PARM_REAL);
    EXF_NML_OBCS.addParm('obcsEperiod',obcsPeriod,PARM_REAL);

    EXF_NML_OBCS.addParm('siobNperiod',obcsPeriod,PARM_REAL);
    EXF_NML_OBCS.addParm('siobEperiod',obcsPeriod,PARM_REAL);

    
  else

  %  obcsNstartdate1   = str2num([num2str(start_year),'0101']);
  %  obcsNstartdate2   = 000000;
   obcsNperiod       = -12;

  %  obcsEstartdate1   = str2num([num2str(start_year),'0101']);
  %  obcsEstartdate2   = 000000;
   obcsEperiod       = -12;


  %  siobNstartdate1   = str2num([num2str(start_year),'0101']);
  %  siobNstartdate2   = 000000;
   siobNperiod       = -12;

  %  siobEstartdate1   = str2num([num2str(start_year),'0101']);
  %  siobEstartdate2   = 000000;
   siobEperiod       = -12;



  %   EXF_NML_OBCS.addParm('obcsNstartdate1',obcsNstartdate1,PARM_INT);
  %   EXF_NML_OBCS.addParm('obcsNstartdate2',obcsNstartdate2,PARM_INT);
  %   
  %   EXF_NML_OBCS.addParm('obcsEstartdate1',obcsEstartdate1,PARM_INT);
  %   EXF_NML_OBCS.addParm('obcsEstartdate2',obcsEstartdate2,PARM_INT);
  %  
  %   EXF_NML_OBCS.addParm('siobEstartdate1',siobEstartdate1,PARM_INT);
  %   EXF_NML_OBCS.addParm('siobEstartdate2',siobEstartdate2,PARM_INT);
  %   
  %   EXF_NML_OBCS.addParm('siobNstartdate1',siobNstartdate1,PARM_INT);
  %   EXF_NML_OBCS.addParm('siobNstartdate2',siobNstartdate2,PARM_INT);

    EXF_NML_OBCS.addParm('obcsNperiod',obcsNperiod,PARM_REAL);
    EXF_NML_OBCS.addParm('obcsEperiod',obcsEperiod,PARM_REAL);
  %   EXF_NML_OBCS.addParm('obcsWperiod',obcsWperiod,PARM_REAL);

    EXF_NML_OBCS.addParm('siobNperiod',siobNperiod,PARM_REAL);
    EXF_NML_OBCS.addParm('siobEperiod',siobEperiod,PARM_REAL);
  %   EXF_NML_OBCS.addParm('siobWperiod',siobWperiod,PARM_REAL);

  end

     
%%%%%%%%%%%%%%%%%%%% ADD PARAMETERS

  EXF_NML_01.addParm('exf_albedo',exf_albedo,PARM_INT);
  EXF_NML_01.addParm('exf_scal_BulkCdn',exf_scal_BulkCdn,PARM_REAL);
  EXF_NML_01.addParm('exf_iprec',exf_iprec,PARM_INT);
  EXF_NML_01.addParm('useExfYearlyFields',useExfYearlyFields,PARM_BOOL);
  EXF_NML_01.addParm('useExfCheckRange',useExfCheckRange,PARM_BOOL);
  EXF_NML_01.addParm('useRelativeWind',useRelativeWind,PARM_BOOL);
  EXF_NML_01.addParm('repeatPeriod',repeatPeriod,PARM_REAL);
  EXF_NML_03.addParm('exf_offset_atemp',exf_offset_atemp,PARM_REAL);
%   EXF_NML_03.addParm('exf_inscal_runoff',exf_inscal_runoff,PARM_REAL);
  
  
  EXF_NML_02.addParm('apressurefile',apressurefile,PARM_STR);
  EXF_NML_02.addParm('atempfile',atempfile,PARM_STR);
  EXF_NML_02.addParm('aqhfile',aqhfile,PARM_STR);
  EXF_NML_02.addParm('uwindfile',uwindfile,PARM_STR);
  EXF_NML_02.addParm('vwindfile',vwindfile,PARM_STR);
  EXF_NML_02.addParm('precipfile',precipfile,PARM_STR);
  EXF_NML_02.addParm('swdownfile',swdownfile,PARM_STR);
  EXF_NML_02.addParm('lwdownfile',lwdownfile,PARM_STR);
%   EXF_NML_02.addParm('runofffile',runofffile,PARM_STR);


  EXF_NML_02.addParm('apressurestartdate1',apressurestartdate1,PARM_INT);
  EXF_NML_02.addParm('apressurestartdate2',apressurestartdate2,PARM_INT);
  EXF_NML_02.addParm('apressureperiod',apressureperiod,PARM_REAL);

  EXF_NML_02.addParm('atempstartdate1',atempstartdate1,PARM_INT);
  EXF_NML_02.addParm('atempstartdate2',atempstartdate2,PARM_INT);
  EXF_NML_02.addParm('atempperiod',atempperiod,PARM_REAL);
  
  EXF_NML_02.addParm('aqhstartdate1',aqhstartdate1,PARM_INT);
  EXF_NML_02.addParm('aqhstartdate2',aqhstartdate2,PARM_INT);
  EXF_NML_02.addParm('aqhperiod',aqhperiod,PARM_REAL);
  
  EXF_NML_02.addParm('uwindstartdate1',uwindstartdate1,PARM_INT);
  EXF_NML_02.addParm('uwindstartdate2',uwindstartdate2,PARM_INT);
  EXF_NML_02.addParm('uwindperiod',uwindperiod,PARM_REAL);
  
  EXF_NML_02.addParm('vwindperiod',vwindperiod,PARM_REAL);
  EXF_NML_02.addParm('vwindstartdate1',vwindstartdate1,PARM_INT);
  EXF_NML_02.addParm('vwindstartdate2',vwindstartdate2,PARM_INT);
  
  EXF_NML_02.addParm('precipperiod',precipperiod,PARM_REAL);
  EXF_NML_02.addParm('precipstartdate1',precipstartdate1,PARM_INT);
  EXF_NML_02.addParm('precipstartdate2',precipstartdate2,PARM_INT);
  
  EXF_NML_02.addParm('swdownperiod',swdownperiod,PARM_REAL);
  EXF_NML_02.addParm('swdownstartdate1',swdownstartdate1,PARM_INT);
  EXF_NML_02.addParm('swdownstartdate2',swdownstartdate2,PARM_INT);
  
  EXF_NML_02.addParm('lwdownperiod',lwdownperiod,PARM_REAL);
  EXF_NML_02.addParm('lwdownstartdate1',lwdownstartdate1,PARM_INT);
  EXF_NML_02.addParm('lwdownstartdate2',lwdownstartdate2,PARM_INT);
  
%   EXF_NML_02.addParm('runoffperiod',runoffperiod,PARM_REAL);

  EXF_NML_04.addParm('precip_lon0',precip_lon0,PARM_REAL);
  EXF_NML_04.addParm('precip_lon_inc',precip_lon_inc,PARM_REAL);
  EXF_NML_04.addParm('precip_lat0',precip_lat0,PARM_REAL);
  EXF_NML_04.addParm('precip_lat_inc',precip_lat_inc,PARM_REALS);
  EXF_NML_04.addParm('precip_nlon',precip_nlon,PARM_INT);
  EXF_NML_04.addParm('precip_nlat',precip_nlat,PARM_INT);

  EXF_NML_04.addParm('atemp_lon0',atemp_lon0,PARM_REAL);
  EXF_NML_04.addParm('atemp_lon_inc',atemp_lon_inc,PARM_REAL);
  EXF_NML_04.addParm('atemp_lat0',atemp_lat0,PARM_REAL);
  EXF_NML_04.addParm('atemp_lat_inc',atemp_lat_inc,PARM_REALS);
  EXF_NML_04.addParm('atemp_nlon',atemp_nlon,PARM_INT);
  EXF_NML_04.addParm('atemp_nlat',atemp_nlat,PARM_INT);

  EXF_NML_04.addParm('apressure_lon0',apressure_lon0,PARM_REAL);
  EXF_NML_04.addParm('apressure_lon_inc',apressure_lon_inc,PARM_REAL);
  EXF_NML_04.addParm('apressure_lat0',apressure_lat0,PARM_REAL);
  EXF_NML_04.addParm('apressure_lat_inc',apressure_lat_inc,PARM_REALS);
  EXF_NML_04.addParm('apressure_nlon',apressure_nlon,PARM_INT);
  EXF_NML_04.addParm('apressure_nlat',apressure_nlat,PARM_INT);

  EXF_NML_04.addParm('aqh_lon0',aqh_lon0,PARM_REAL);
  EXF_NML_04.addParm('aqh_lon_inc',aqh_lon_inc,PARM_REAL);
  EXF_NML_04.addParm('aqh_lat0',aqh_lat0,PARM_REAL);
  EXF_NML_04.addParm('aqh_lat_inc',aqh_lat_inc,PARM_REALS);
  EXF_NML_04.addParm('aqh_nlon',aqh_nlon,PARM_INT);
  EXF_NML_04.addParm('aqh_nlat',aqh_nlat,PARM_INT);

  EXF_NML_04.addParm('uwind_lon0',uwind_lon0,PARM_REAL);
  EXF_NML_04.addParm('uwind_lon_inc',uwind_lon_inc,PARM_REAL);
  EXF_NML_04.addParm('uwind_lat0',uwind_lat0,PARM_REAL);
  EXF_NML_04.addParm('uwind_lat_inc',uwind_lat_inc,PARM_REALS);
  EXF_NML_04.addParm('uwind_nlon',uwind_nlon,PARM_INT);
  EXF_NML_04.addParm('uwind_nlat',uwind_nlat,PARM_INT);

  EXF_NML_04.addParm('vwind_lon0',vwind_lon0,PARM_REAL);
  EXF_NML_04.addParm('vwind_lon_inc',vwind_lon_inc,PARM_REAL);
  EXF_NML_04.addParm('vwind_lat0',vwind_lat0,PARM_REAL);
  EXF_NML_04.addParm('vwind_lat_inc',vwind_lat_inc,PARM_REALS);
  EXF_NML_04.addParm('vwind_nlon',vwind_nlon,PARM_INT);
  EXF_NML_04.addParm('vwind_nlat',vwind_nlat,PARM_INT);

  EXF_NML_04.addParm('swdown_lon0',swdown_lon0,PARM_REAL);
  EXF_NML_04.addParm('swdown_lon_inc',swdown_lon_inc,PARM_REAL);
  EXF_NML_04.addParm('swdown_lat0',swdown_lat0,PARM_REAL);
  EXF_NML_04.addParm('swdown_lat_inc',swdown_lat_inc,PARM_REALS);
  EXF_NML_04.addParm('swdown_nlon',swdown_nlon,PARM_INT);
  EXF_NML_04.addParm('swdown_nlat',swdown_nlat,PARM_INT);

  EXF_NML_04.addParm('lwdown_lon0',lwdown_lon0,PARM_REAL);
  EXF_NML_04.addParm('lwdown_lon_inc',lwdown_lon_inc,PARM_REAL);
  EXF_NML_04.addParm('lwdown_lat0',lwdown_lat0,PARM_REAL);
  EXF_NML_04.addParm('lwdown_lat_inc',lwdown_lat_inc,PARM_REALS);
  EXF_NML_04.addParm('lwdown_nlon',lwdown_nlon,PARM_INT);
  EXF_NML_04.addParm('lwdown_nlat',lwdown_nlat,PARM_INT);

%   EXF_NML_04.addParm('runoff_lon0',runoff_lon0,PARM_REAL);
%   EXF_NML_04.addParm('runoff_lon_inc',runoff_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('runoff_lat0',runoff_lat0,PARM_REAL);
%   EXF_NML_04.addParm('runoff_lat_inc',runoff_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('runoff_nlon',runoff_nlon,PARM_INT);
%   EXF_NML_04.addParm('runoff_nlat',runoff_nlat,PARM_INT);


  %%z% Create the data.exf file
  write_data_exf(inputpath,EXF_PARM,listterm,realfmt);


  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%% CAL%%% %%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CAL    SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  
  cal_parm01 = parmlist;
  CAL_PARM = {cal_parm01};
  
  cal_parm01.addParm('TheCalendar','gregorian',PARM_STR);
  cal_parm01.addParm('startDate_1',str2num([num2str(start_year),'0101']),PARM_INT);
  cal_parm01.addParm('startDate_2',000000,PARM_INT);
  cal_parm01.addParm('calendarDumps',false,PARM_BOOL);
  
  %%z% Create the data.cal file
  write_data_cal(inputpath,CAL_PARM,listterm,realfmt);
  
  
  
  
  
  



  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%% KPP%%% %%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% KPP    SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%
  


%%%%%%%% WRITTEN IN DEFAULTS/INPUT


  
  %%% To store parameter names and values
  
  kpp_parm01 = parmlist;
  KPP_PARM = {kpp_parm01};
  KPPmixingMaps   = false;
  KPPwriteState   = true;
  
  kpp_parm01.addParm('KPPmixingMaps',KPPmixingMaps,PARM_BOOL);  
  kpp_parm01.addParm('KPPwriteState',KPPwriteState,PARM_BOOL)
  
  %%z% Create the data.kpp file
  write_data_kpp(inputpath,KPP_PARM,listterm,realfmt);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% RBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  useRBCS = false;
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  rbcs_parm01 = parmlist;
  rbcs_parm02 = parmlist;
  RBCS_PARM = {rbcs_parm01,rbcs_parm02};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  useRBCtemp = useRBCS;
  useRBCsalt = useRBCS;
  useRBCuVel = false;
  useRBCvVel = false;
%   tauRelaxT = 60*t1day;
%   tauRelaxS = 60*t1day;  
  tauRelaxT = 30*t1day;
  tauRelaxS = 30*t1day; 
  rbcs_parm01.addParm('useRBCtemp',useRBCtemp,PARM_BOOL);
  rbcs_parm01.addParm('useRBCsalt',useRBCsalt,PARM_BOOL);
  rbcs_parm01.addParm('useRBCuVel',useRBCuVel,PARM_BOOL);
  rbcs_parm01.addParm('useRBCvVel',useRBCvVel,PARM_BOOL);
  rbcs_parm01.addParm('tauRelaxT',tauRelaxT,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxS',tauRelaxS,PARM_REAL);  
  
  if (useRBCS)
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% RELAXATION TEMPERATURE/SALINITY %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

    %%% Restore salinity everywhere under ice shelves to ISW salinity
    salt_ISW = 34.5;
  %   salt_ISW = 34;
    salt_relax = salt_ISW*ones(Nx,Ny,Nr);

    %%% Restore temperature to local freezing point under ice shelves
  %   Pressure = repmat(reshape(-rho0*g*zz/Pa1dbar,[1 1 Nr]),[Nx Ny 1]);
  %   temp_relax = .0901 - .0575*salt_relax - (7.61e-4 * Pressure);
    %%%% Restore temperature to surface freezing point 

    Pressure = repmat(reshape(-rho0*g*zz(1)/Pa1dbar,[1 1 1]),[Nx Ny 1]);
    temp_relax = .0901 - .0575*salt_relax - (7.61e-4 * Pressure);


    %%% Save as parameters
    writeDataset(temp_relax,fullfile(inputpath,'sponge_temp.bin'),ieee,prec); 
    rbcs_parm01.addParm('relaxTFile','sponge_temp.bin',PARM_STR);
    writeDataset(salt_relax,fullfile(inputpath,'sponge_salt.bin'),ieee,prec); 
    rbcs_parm01.addParm('relaxSFile','sponge_salt.bin',PARM_STR);  


    %%%%%%%%%%%%%%%%%%%%%  
    %%%%% RBCS MASK %%%%%
    %%%%%%%%%%%%%%%%%%%%%  

    %%% Mask is zero everywhere by default, i.e. no relaxation
    mskT=zeros(Nx,Ny,Nr);
    mskS=zeros(Nx,Ny,Nr);     

    %%% If a sponge BC is required, gradually relax in the gridpoints 
    %%% approaching the wall (no relaxation at the wall)   
    for i=1:Nx
      for j=1:Ny      
        if (SHELFICEtopo(i,j) < 0)
          mskT(i,j,:) = 1;
          mskS(i,j,:) = 1;
        end                  
      end       
    end   

    %%% Save as parameters
    writeDataset(mskT,fullfile(inputpath,'rbcs_temp_mask.bin'),ieee,prec); 
    rbcs_parm01.addParm('relaxMaskFile(1)','rbcs_temp_mask.bin',PARM_STR);
    writeDataset(mskS,fullfile(inputpath,'rbcs_salt_mask.bin'),ieee,prec); 
    rbcs_parm01.addParm('relaxMaskFile(2)','rbcs_salt_mask.bin',PARM_STR);


  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.rbcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
  %%% Creates the 'data.rbcs' file
  write_data_rbcs(inputpath,RBCS_PARM,listterm,realfmt);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  layers_parm01 = parmlist;
  LAYERS_PARM = {layers_parm01};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Define parameters for layers package %%%
  
  %%% Number of fields for which to calculate layer fluxes
  layers_maxNum = 1;
  
  %%% Specify potential temperature
  layers_name = 'RHO';  
  
  %%% Potential temperature bounds for layers  
  layers_bounds = [30.5:1:36.5 36.6:0.1:36.8 36.9:0.01:37.4 37.42:0.02:37.6] - 9.38; 
  
  %%% Reference level for calculation of potential density
  layers_krho = 1;    
  
  %%% If set true, the GM bolus velocity is added to the calculation
  layers_bolus = false;  
   
  %%% Layers
  layers_parm01.addParm('layers_bounds',layers_bounds,PARM_REALS); 
  layers_parm01.addParm('layers_krho',layers_krho,PARM_INT); 
  layers_parm01.addParm('layers_name',layers_name,PARM_STR); 
  layers_parm01.addParm('layers_bolus',layers_bolus,PARM_BOOL); 

  %%z% Create the data.layers file
  write_data_layers(inputpath,LAYERS_PARM,listterm,realfmt);
  
  %%% Create the LAYERS_SIZE.h file
  createLAYERSSIZEh(codepath,length(layers_bounds)-1,layers_maxNum); 
  
  
  
  
  
  
  
  
  
     
  
  
  
  
  
  
 
 
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
    
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  diag_parm01 = parmlist;
  diag_parm02 = parmlist;
  DIAG_PARM = {diag_parm01,diag_parm02};
  diag_matlab_parm01 = parmlist;
  DIAG_MATLAB_PARM = {diag_matlab_parm01}; %%% Matlab parameters need to be different
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%% Stores total number of diagnostic quantities
  ndiags = 0;
  
  
  diag_fields_avg = ...
  { ...
    'UVEL','VVEL','WVEL'...      
    'THETA' ... 
    'SALT','ETAN' ... 
    'SIheff','SIarea','SIhsnow','SItices','SIhsalt' ... %%ice diagnostics
    'SIuice','SIvice' ...
    'oceFWflx','oceSflux','oceQnet','oceTAUX','oceTAUY','TFLUX','SFLUX'...
    'SHIfwFlx','SHIhtFlx','SHIUDrag','SHIVDrag' ...
    'EXFtaux','EXFtauy','EXFlwnet','EXFswnet','EXFlwdn','EXFswdn','EXFqnet','EXFhs','EXFhl','EXFevap','EXFpreci','EXFatemp' ...
    'SIqnet','SIqsw','SIatmQnt','SItflux','SIaaflux','SIhl','SIqneto','SIqneti','SIempmr','SIatmFW','SIsnPrcp','SIactLHF','SIacSubl'
  };

  if (use_extended_diags)
    diag_fields_avg = {diag_fields_avg{:}, ...
      'UVELSLT','VVELSLT','WVELSLT', ... %%% Salt fluxes
      'UVELTH','VVELTH','WVELTH', ... %%% Temperature fluxes  
      'UVELSQ','VVELSQ','WVELSQ', ... %%% For Kinetic Energy
      'UV_VEL_Z','WU_VEL','WV_VEL', ... %%% Momentum fluxes
      'SALTSQ','THETASQ','THSLT' %%% Thermodynamic variances
    };
  end     
  
  if (use_layers)
    diag_fields_avg = {diag_fields_avg{:}, ...
      'LaUH1RHO', 'LaVH1RHO', ... %%% LAYERS fluxes
      'LaHw1RHO', 'LaHs1RHO' ... %%% LAYERS thicknesses
    };
  end  
  
  numdiags_avg = length(diag_fields_avg);  
  diag_freq_avg = t1month; %%% Approximately monthly output
%   diag_freq_avg = 10*deltaT;
%   diag_freq_avg = 30*t1day;
%   diag_freq_avg = .0417*t1day;
%  diag_freq_avg = 5*t1year;
  diag_phase_avg = 0;    
     
  diag_parm01.addParm('diag_mnc',true,PARM_BOOL);  
  for n=1:numdiags_avg    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields{1,',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames{',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL);  
    
  end
  
  diag_fields_inst = ...
  {...      
  %'ETAN'...
%     'UVEL','VVEL','WVEL'...      
%     'THETA' ... 
%     'SALT','ETAN' ... 
%     'SIheff','SIarea' ... %%ice diagnostics
%     'SIuice','SIvice' ...
%     'PHIBOT','oceFWflx','oceSflux','oceQnet'...
%     'SHIfwFlx','SHIhtFlx','SHIUDrag','SHIVDrag' ...
  }; 
  
  numdiags_inst = length(diag_fields_inst);    
%   diag_freq_inst = 0.1*t1year;
%    diag_freq_inst = .0417*t1day;
%   diag_freq_inst = 30*t1day;
%   diag_freq_inst = 10*deltaT;
  diag_freq_inst = t1day/24;
  diag_phase_inst = 0;
  
  for n=1:numdiags_inst    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);  
    
  end
  
  %%% Create the data.diagnostics file
  write_data_diagnostics(inputpath,DIAG_PARM,listterm,realfmt);
  
  %%% Create the DIAGNOSTICS_SIZE.h file
  createDIAGSIZEh(codepath,ndiags,Nr);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%% PACKAGES %%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  
  packages = parmlist;
  PACKAGE_PARM = {packages};  
  
  packages.addParm('useDiagnostics',true,PARM_BOOL);    
  packages.addParm('useKPP',true,PARM_BOOL);
  packages.addParm('useOBCS',true,PARM_BOOL);  
  packages.addParm('useRBCS',false,PARM_BOOL);  
  packages.addParm('useEXF',true,PARM_BOOL);        
  packages.addParm('useCAL',true,PARM_BOOL);         
  packages.addParm('useSEAICE',true,PARM_BOOL);      
  packages.addParm('useSHELFICE',true,PARM_BOOL);      
  packages.addParm('useFRAZIL',false,PARM_BOOL);        
  packages.addParm('useLayers',use_layers,PARM_BOOL); 
  packages.addParm('useBBL',false,PARM_BOOL); 
  packages.addParm('useGMREDI',false,PARM_BOOL);          
  
  %%% Create the data.pkg file
  write_data_pkg(inputpath,PACKAGE_PARM,listterm,realfmt);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE PARAMETERS TO A MATLAB FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates a matlab file defining all input parameters
%   write_matlab_params(inputpath,[PARM DIAG_MATLAB_PARM],realfmt);
%   write_matlab_params(inputpath,[PARM SHELFICE_PARM EXF_PARM DIAG_MATLAB_PARM],realfmt);
  write_matlab_params(inputpath,[PARM RBCS_PARM KPP_PARM OBCS_PARM SHELFICE_PARM SEAICE_PARM EXF_PARM CAL_PARM LAYERS_PARM DIAG_MATLAB_PARM PACKAGE_PARM],realfmt);
  
end
