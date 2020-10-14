%%%
%%% createND2.m
%%%
%%% Creates a neutral density (of the second kind) output file from 
%%% potential temperature and salinity output files.
%%%
function createND2 (expdir,expname,iter)




  %%%%%%%%%%%%%
  %%% Setup %%%
  %%%%%%%%%%%%%

  addpath ../NeutDens;
  addpath ../NeutDens/matlab-interface;
  addpath ../GSW;
  addpath ../GSW/library/
 
  %%% Experiment directories
  exppath = fullfile(expdir,expname);
  inputpath = fullfile(exppath,'input');
  resultspath = fullfile(exppath,'results');
  
  %%% Check whether the file already exists: if so, don't create it
  ncname = fullfile(resultspath,['ND2.',num2str(iter,'%.10d'),'.dat']);
  if (exist(ncname))
    ['WARNING: ',ncname,' not generated because the file already exists']
    return;
  end
  
  %%% Read temperature and salinity data files
  pt = rdmdsWrapper(fullfile(resultspath,'THETA'),iter);         
  ss = rdmdsWrapper(fullfile(resultspath,'SALT'),iter);         
   
  %%% Load parameters used for this experiment
  run(fullfile(inputpath,'params.m'));
  
  %%% Grids
  XC = rdmds(fullfile(resultspath,'XC'));
  YC = rdmds(fullfile(resultspath,'YC'));  
  RC = rdmds(fullfile(resultspath,'RC'));
  hFacC = rdmds(fullfile(resultspath,'hFacC'));
  
  %%% Pressure is just Boussinesq hydrostatic reference pressu
  pp_ref = -gravity*rhoConst*squeeze(RC)/1e4;

  %%% Grid dimensions (not specified explicitly in params.m)
  Nx = length(delX);
  Ny = length(delY);
  Nr = length(delR);
 
  %%% To store neutral density
  gg = zeros(Nx,Ny,Nr);
  
  
  
  %%% TODO remove
  tstart = tic();
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Compute neutral densities %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Neutral density computed directly from mean state via JM97 algorithm
  for i=1:Nx
    
    for j=1:Ny
      
      %%% Min and max indices of wet cells in this water column
      kmin = find(hFacC(i,j,:)>0,1,'first');
      kmax = find(hFacC(i,j,:)>0,1,'last');
      cast_len = kmax-kmin+1;
      
      if (isempty(kmin))
        continue;
      end
        
      ss_vec = squeeze(ss(i,j,kmin:kmax));
      pt_vec = squeeze(pt(i,j,kmin:kmax));
      pp_vec = squeeze(pp_ref(kmin:kmax));
      lon = XC(i,j);
      lat = YC(i,j);
      if (lat <= -80)
        lat = -79.9;
      end
      if (lon <= -64)
        lon = -63.9;
      end
      ssa_vec = gsw_SA_from_SP(ss_vec,pp_vec,lon,lat);  
      ttc_vec = gsw_CT_from_pt(ssa_vec,pt_vec);
      ttis_vec = gsw_t_from_CT(ssa_vec,ttc_vec,pp_vec);    
      [gg_vec,dg_lo,dg_hi,wts,sc,tc,pc] = gamma_n(ss_vec,ttis_vec,pp_vec,lon,lat);  
      gg(i,j,kmin:kmax) = reshape(gg_vec,[1 1 cast_len]);

    end
  end
  gg(gg<0) = NaN; %%% Flag missing values
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Loop through casts to fill in missing densities %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for i=1:Nx
    i
    for j=1:Ny
      
      %%% Min and max indices of wet cells in this water column
      kmin = find(hFacC(i,j,:)>0,1,'first');
      kmax = find(hFacC(i,j,:)>0,1,'last');
    
      %%% Search for NaNs
      k_firstnan = [];
      k_lastnan = [];
      for k=kmin:kmax

        %%% NaNs indicate that the search for a matching density level in
        %%% this parcel's reference cast failed
        if (isnan(gg(i,j,k)))
          if (isempty(k_firstnan))     
            k_firstnan = k;
          end
          if (k~=kmax)
            continue;
          else
            k_lastnan = k;
          end
        else         
          if (isempty(k_firstnan))
            continue;
          else     
            k_lastnan = k-1;
          end             
        end

        %%% Just fix b-factor in this case because it doesn't matter too much
        %%% here
        bfac = 2;
        
        %%% Whole water column is NaNs  
        if ((k_firstnan == kmin) && (k_lastnan == kmax))
          
          %%% Do nothing
          break;

        %%% Block of NaNs at top of water column
        elseif (k_firstnan == kmin)

          %%% Extrapolate recursively upward
          for l=k_lastnan:-1:k_firstnan
            pp_mid = 0.5*(pp_ref(l+1)+pp_ref(l));
            drho = densjmd95(ss(i,j,l),pt(i,j,l),pp_mid) ...
                - densjmd95(ss(i,j,l+1),pt(i,j,l+1),pp_mid);
              gg(i,j,l) = gg(i,j,l+1) + drho .* bfac;
          end

        %%% Block of NaNs at bottom of water column
        elseif (k_lastnan == kmax)

          %%% Extrapolate recursively downward
          for l=k_firstnan:1:k_lastnan
            pp_mid = 0.5*(pp_ref(l-1)+pp_ref(l));
            drho = densjmd95(ss(i,j,l),pt(i,j,l),pp_mid) ...
                - densjmd95(ss(i,j,l-1),pt(i,j,l-1),pp_mid);
              gg(i,j,l) = gg(i,j,l-1) + drho .* bfac;
          end

        %%% Block of NaNs in middle of water column
        else            

          %%% Extrapolate upward and downward and take the mean
          for l=k_firstnan:1:k_lastnan
            pp_mid = 0.5*(pp_ref(k_firstnan-1)+pp_ref(l));
            drho1 = densjmd95(ss(i,j,l),pt(i,j,l),pp_mid) ...
                - densjmd95(ss(i,j,k_firstnan-1),pt(i,j,k_firstnan-1),pp_mid);
            pp_mid = 0.5*(pp_ref(k_lastnan+1)+pp_ref(l));
            drho2 = densjmd95(ss(i,j,l),pt(i,j,l),pp_mid) ...
                - densjmd95(ss(i,j,k_lastnan+1),pt(i,j,k_lastnan+1),pp_mid);
            gg(i,j,l) = 0.5 * ( gg(i,j,k_firstnan-1) + drho1 .* bfac ...
                                  + gg(i,j,k_lastnan+1) + drho2 .* bfac );
          end

        end 
        
        k_firstnan = [];
        k_lastnan = [];

      end
      
    end
  end
  
  
  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Inpaint to fill in any remaining NaNs %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  gg(gg<20) = NaN;
  for k = 1:Nr
    gg(:,:,k) = inpaint_nans(gg(:,:,k));
  end
  gg(hFacC==0) = NaN;
  
  
  
  
  
  %%% TODO remove
  toc(tstart)
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Write neutral density to a file %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fid = fopen(ncname,'w','b');
  fwrite(fid,gg,'real*8');
  fclose(fid);
  
  
  
  
  

end
