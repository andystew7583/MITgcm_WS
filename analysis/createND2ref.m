%%%
%%% createND2ref.m
%%%
%%% Computes JM97 neutral density using the time-mean model state, 
%%% for comparison with the reference ND1 variable.
%%%


addpath ../NeutDens;
addpath ../NeutDens/matlab-interface;
addpath ../GSW;
addpath ../GSW/library/



%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
expname = 'hires_seq_onethird_RTOPO2';
loadexp;






%%%%%%%%%%%%%%%%%%%
%%%%% OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%

%%% Time frame over which to average thermodynamic variables to create
%%% climatology
tmin = 18.05*86400*365;
tmax = 27.05*86400*365;








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% COMPUTE REFERENCE STATE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(endTime/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Time-average temperature and salinity
pt_ref = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
ss_ref = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);

%%% Pressure is just Boussinesq hydrostatic reference pressure
pp_ref = -gravity*rhoConst*repmat(RC,[Nx Ny 1])/1e4;

%%% Remove topography
pt_ref(hFacC==0) = NaN;
ss_ref(hFacC==0) = NaN;
 
%%% To store neutral density
gg_ref = zeros(Nx,Ny,Nr);



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
    cast_len = kmax-kmin+1;

    if (isempty(kmin))
      continue;
    end

    ss_vec = squeeze(ss_ref(i,j,kmin:kmax));
    pt_vec = squeeze(pt_ref(i,j,kmin:kmax));
    pp_vec = squeeze(pp_ref(i,j,kmin:kmax));
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
    gg_ref(i,j,kmin:kmax) = reshape(gg_vec,[1 1 cast_len]);

  end
end
gg_ref(gg_ref<0) = NaN; %%% Flag missing values






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
      if (isnan(gg_ref(i,j,k)))
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
          drho = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
              - densjmd95(ss_ref(i,j,l+1),pt_ref(i,j,l+1),pp_mid);
            gg_ref(i,j,l) = gg_ref(i,j,l+1) + drho .* bfac;
        end

      %%% Block of NaNs at bottom of water column
      elseif (k_lastnan == kmax)

        %%% Extrapolate recursively downward
        for l=k_firstnan:1:k_lastnan
          pp_mid = 0.5*(pp_ref(l-1)+pp_ref(l));
          drho = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
              - densjmd95(ss_ref(i,j,l-1),pt_ref(i,j,l-1),pp_mid);
            gg_ref(i,j,l) = gg_ref(i,j,l-1) + drho .* bfac;
        end

      %%% Block of NaNs in middle of water column
      else            

        %%% Extrapolate upward and downward and take the mean
        for l=k_firstnan:1:k_lastnan
          pp_mid = 0.5*(pp_ref(k_firstnan-1)+pp_ref(l));
          drho1 = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
              - densjmd95(ss_ref(i,j,k_firstnan-1),pt_ref(i,j,k_firstnan-1),pp_mid);
          pp_mid = 0.5*(pp_ref(k_lastnan+1)+pp_ref(l));
          drho2 = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
              - densjmd95(ss_ref(i,j,k_lastnan+1),pt_ref(i,j,k_lastnan+1),pp_mid);
          gg_ref(i,j,l) = 0.5 * ( gg_ref(i,j,k_firstnan-1) + drho1 .* bfac ...
                                + gg_ref(i,j,k_lastnan+1) + drho2 .* bfac );
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

gg_ref(gg_ref<20) = NaN;
for k = 1:Nr
  gg_ref(:,:,k) = inpaint_nans(gg_ref(:,:,k));
end
gg_ref(hFacC==0) = NaN;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write neutral density to a file %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Finally, save the new reference dataset
save(fullfile('products',[expname,'_ND2.mat']),'ss_ref','pt_ref','pp_ref','gg_ref');
