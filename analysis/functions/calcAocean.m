% %%%
% %%% calcAocean.m
% %%%
% %%% Calculates ocean area above density surfaces using the MITgcm 'layers' package.
% %%%
% 
% %%% Options
% expdir = '../experiments';
% % expname = 'hires_seq_onethird_RTOPO2';
% % tmin = 18.05;
% % tmax = 27.05;
% % expname = 'hires_seq_onesixth_RTOPO2';
% % tmin = 9.05;
% % tmax = 18.05;
% % expname = 'hires_seq_onetwelfth_RTOPO2';
% % tmin = 1.05;
% % tmax = 9.05;
% expname = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% tmin = 1.01;
% tmax = 7.01;
% 
% %%% Load experiment
% loadexp;
% 
% %%% Set true to deform coordinates in the cavity
% deform_cavity = false;
% 
% %%% Set true to use barotropic streamfunction as the coordinate system
% use_PsiBT = false;
% 
% 
% %%% Select density variable in which to compute isopycnal fluxes
% densvar = 'PD0';
% % densvar = 'ND1';
% % densvar = 'ND2';
% % densvar = 'PT';
% 
% %%% Density bins for MOC calculation  
% densvar = 'PD0';
% dens_levs = layers_bounds;
% Nd = length(dens_levs)-1;
% p_ref = -rhoConst*gravity*RC(1)/1e4; %%% Reference pressure for surface-referenced potential density
% 
% %%% Define coordinate system for integrating to compute streamfunction
% if (use_PsiBT)
% 
%   infname = [expname,'_TSfluxes'];
%   load(fullfile('products',infname),'uvel_tavg');
% 
%   %%% Calculate depth-averaged zonal velocity
%   UU = sum(uvel_tavg.*repmat(DRF,[Nx Ny 1]).*hFacW,3);
%   clear('uvel_tavg');
%   
%   %%% Calculate barotropic streamfunction
%   Psi = zeros(Nx+1,Ny+1);
%   Psi(2:Nx+1,2:Ny+1) = -cumsum(UU.*DYG,2);
%   Psi = Psi(1:Nx,1:Ny);
%   
%   %%% Interpolate to cell centers
%   ETA = 0.25*(Psi(1:Nx,1:Ny)+Psi([2:Nx 1],1:Ny)+Psi(1:Nx,[2:Ny 1])+Psi([2:Nx 1],[2:Ny 1]))/1e6;
%   
%   %%% Streamunction grid for flux calculation
%   eta = -2:.1:10;
%   Neta = length(eta);
% 
% else
% 
%   ETA = defineMOCgrid(XC,YC,SHELFICEtopo,bathy,deform_cavity);
%   eta = -9:.1:11;
%   Neta = length(eta);
% 
% end
% 
% %%% Frequency of diagnostic output - should match that specified in
% %%% data.diagnostics.
% dumpFreq = abs(diag_frequency(1));
% nDumps = round(endTime/dumpFreq);
% dumpIters = round((1:nDumps)*dumpFreq/deltaT);
% dumpIters = dumpIters(dumpIters > nIter0);
% nDumps = length(dumpIters);
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%
% %%% GRIDS %%%
% %%%%%%%%%%%%%
% 
% %%% 3D horizontal grid spacing matrices
% DXG_3D = repmat(DXG,[1 1 Nd]);
% DYG_3D = repmat(DYG,[1 1 Nd]);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% ISOPYCNAL AREA CALCULATION %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Load time-mean isopycnal thicknesses
% hisop_u = readIters(exppath,'LaHw1RHO',dumpIters,deltaT,tmin*86400*365,tmax*86400*365,Nx,Ny,Nd);
% hisop_v = readIters(exppath,'LaHs1RHO',dumpIters,deltaT,tmin*86400*365,tmax*86400*365,Nx,Ny,Nd);
% 
% %%% Compute area below each isopycnal
% Aisop = -calcIsopStreamfunction(...
%   hisop_u,hisop_v, ...
%   Nx,Ny,Neta,Nd, ...  
%   DXG_3D,DYG_3D,ETA,eta);

%%% Area of ocean below each z-level (for mapping streamfunction back to
%%% real space
Ax = hFacW.*repmat(DRF,[Nx Ny 1]);
Ay = hFacS.*repmat(DRF,[Nx Ny 1]);
Aocean = -calcIsopStreamfunction(...
  Ax,Ay, ...
  Nx,Ny,Neta,Nr, ...  
  repmat(DXG,[1 1 Nr]),repmat(DXG,[1 1 Nr]),ETA,eta);

%%% TODO I think we need to rethink this for the case of area under the ice
%%% shelf


%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

%%% Construct output file name
outfname = [expname,'_Aocean_',densvar];
outfname = [outfname,'.mat'];

%%% Store computed data for later
save(fullfile('products',outfname),'-v7.3', ...
  'eta','ETA','dens_levs', ...
  'Aisop','Aocean');
























function psi = calcIsopStreamfunction(...
  uflux,vflux, ...
  Nx,Ny,Neta,Nd, ...  
  DXG_3D,DYG_3D,ETA,eta)

  %%% Compute horizontal divergence of isopycnal fluxes
  fluxdiv = zeros(Nx,Ny,Nd);
  fluxdiv(1:Nx-1,1:Ny-1,:) = uflux(2:Nx,1:Ny-1,:) .* DYG_3D(2:Nx,1:Ny-1,:) ...
                              - uflux(1:Nx-1,1:Ny-1,:) .* DYG_3D(1:Nx-1,1:Ny-1,:) ...
                              + vflux(1:Nx-1,2:Ny,:) .* DXG_3D(1:Nx-1,2:Ny,:) ...
                              - vflux(1:Nx-1,1:Ny-1,:) .* DXG_3D(1:Nx-1,1:Ny-1,:);
                       
  %%% Integrate flux divergence across lines of constant eta (parallel to FRIS face)
  eflux = zeros(Neta,Nd);
  for m = 1:Neta
    msk = repmat(ETA<eta(m),[1 1 Nd]);
    eflux(m,:) = squeeze(sum(sum(fluxdiv.*msk,1),2));
  end

  %%% Sum fluxes to obtain streamfunction
  psi = zeros(Neta,Nd+1);
  for m=1:Nd  
    psi(:,m) = -sum(eflux(:,m:Nd),2);     
  end

end















