%%%
%%% createND1.m
%%%
%%% Creates a neutral density (of the first kind) output file from 
%%% potential temperature and salinity output files.
%%%
function createND1 (expdir,expname,iter)
 


  %%%%%%%%%%%%%
  %%% Setup %%%
  %%%%%%%%%%%%%
  
  %%% Experiment directories
  exppath = fullfile(expdir,expname);
  inputpath = fullfile(exppath,'input');
  resultspath = fullfile(exppath,'results');
  
  %%% Check whether the file already exists: if so, don't create it
  ncname = fullfile(resultspath,['ND1.',num2str(iter,'%.10d'),'.dat']);
  if (exist(ncname))
    ['WARNING: ',ncname,' not generated because the file already exists']
    return;
  end
  
  %%% Read temperature and salinity data files
  pt = rdmdsWrapper(fullfile(resultspath,'THETA'),iter);         
  ss = rdmdsWrapper(fullfile(resultspath,'SALT'),iter);         
  
  %%% Load ND1 reference data
  load([expname,'_ND1.mat']);
  
  %%% Load parameters used for this experiment
  run(fullfile(inputpath,'params.m'));
  
  %%% Grids
  XC = rdmds(fullfile(resultspath,'XC'));
  YC = rdmds(fullfile(resultspath,'YC'));  
  RC = rdmds(fullfile(resultspath,'RC'));  
  hFacC = rdmds(fullfile(resultspath,'hFacC'));

  %%% Grid dimensions (not specified explicitly in params.m)
  Nx = length(delX);
  Ny = length(delY);
  Nr = length(delR);
  
  %%% To store neutral density
  gg = zeros(Nx,Ny,Nr);
  
  %%% 'b-factor' of Jackett and McDougall 1997
  bfac = 2;
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Compute neutral densities %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Neutral density computed directly from mean state via JM97 algorithm
  for i=1:Nx
    i
    for j=1:Ny
      
      %%% Min and max indices of wet cells in this water column
      kmin = find(hFacC(i,j,:)>0,1,'first');
      kmax = find(hFacC(i,j,:)>0,1,'last');     
      
      if (isempty(kmin))
        continue;
      end
    
      %%% Extract "cast" properties
      ss_cast = ss_ref(i,j,:);
      pt_cast = pt_ref(i,j,:);
      pp_cast = pp_ref(i,j,:);
      gg_cast = gg_ref(i,j,:);
      cast_len = find(~isnan(ss_cast),1,'last');             
            
      %%% Vertical loop
      for k=kmin:kmax 
                    
        %%% "Parcel" properties (just the properties of this grid cell)
        ss_parcel = ss(i,j,k);
        pt_parcel = pt(i,j,k);
        pp_parcel = pp_ref(i,j,k); %%% Reference data are on the same grid as instantaneous data
                         
        %%% Density differences between the "parcel" and the lightest and
        %%% densest points on the reference "cast"
        topdiff = -densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,kmin,kmax,pp_cast(kmin));
        botdiff = densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,kmin,kmax,pp_cast(kmax));

        %%% If the "parcel" falls within the density range of the "cast",
        %%% assign a neutral density label via linear interpolation
        if ((topdiff < 0) && (botdiff < 0))                   
    
          %%% Initial guess for neutral pressure on cast
          pp_neut = pp_ref(i,j,k);
                              
          pp_mean = 0.5*(pp_parcel+pp_neut);
          E = densjmd95(ss_parcel,pt_parcel,pp_mean) - densjmd95(ss_cast(k),pt_cast(k),pp_mean);
          if (E > 0)
            Enext = E;            
            knext = k;
            while (Enext > 0)
              knext = knext + 1;
              pp_neut = pp_cast(knext);
              Eprev = Enext;
%               Enext = densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast(knext),pt_cast(knext),pp_cast(knext),kmin,kmax,pp_neut);
              pp_mean = 0.5*(pp_parcel+pp_neut);
              Enext = densjmd95(ss_parcel,pt_parcel,pp_mean) - densjmd95(ss_cast(knext),pt_cast(knext),pp_mean);
            end
            kprev = knext - 1;
          else
            Eprev = E;            
            kprev = k;
            while (Eprev < 0)
              kprev = kprev - 1;
              pp_neut = pp_cast(kprev);
              Enext = Eprev;
%               Eprev = densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast(kprev),pt_cast(kprev),pp_cast(kprev),cast_len,pp_neut);
              pp_mean = 0.5*(pp_parcel+pp_neut);
              Eprev = densjmd95(ss_parcel,pt_parcel,pp_mean) - densjmd95(ss_cast(kprev),pt_cast(kprev),pp_mean);
            end
            knext = kprev + 1;
          end
          
          %%% Linearly interpolate to find the zero
          pp_neut = (pp_cast(kprev)*Enext - pp_cast(knext)*Eprev) / (Enext - Eprev);

          %%% Find the pressure on the "cast" corresponding to the point that is
          %%% neutral to the "parcel"          
          %options = optimset('TolX',1e-16,'Display','iter');
%           pp_neut = fzero(@(p) densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,cast_len,p), pp_neut);%, options);
          
          %%% Label this point with neutral density. Currently using linear
          %%% interpolation vertically on "cast" data to assign label, which is
          %%% somewhat crude, though our very tightly spaced gridpoints help
          %%% here.
          gg(i,j,k) = interpCast(gg_cast,pp_cast,kmin,kmax,pp_neut);
    
        %%% If the parcel is lightest than the lightest point on the
        %%% cast, or densest than the densest point, assign neutral density
        %%% via linear interpolation.
        else
                    
          if (topdiff > 0)
            ss_diff = ss_cast(kmin);
            pt_diff = pt_cast(kmin);
            pp_diff = pp_cast(kmin);
            gg_diff = gg_cast(kmin);
          else
            ss_diff = ss_cast(kmax);
            pt_diff = pt_cast(kmax);
            pp_diff = pp_cast(kmax);
            gg_diff = gg_cast(kmax);
          end
            
          %%% Calculate T,S,P,dT,dS,alpha,beta,drho all mid-way between this grid cell and the
          %%% reference grid cell                    
          pp_mid = 0.5*(pp_parcel + pp_diff);              
          drho = densjmd95(ss_parcel,pt_parcel,pp_mid) - densjmd95(ss_diff,pt_diff,pp_mid);

          %%% Assign density label via linear extrapolation    
          gg(i,j,k) = gg_diff + drho * bfac;  
          
        end     
        
      end            
      
    end
  
  end
  
    
    
    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Write neutral density to a file %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Write neutral density to a file
  fid = fopen(ncname,'w','b');
  fwrite(fid,gg,'real*8');
  fclose(fid);
  

end
