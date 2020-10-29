
%%%
%%% createND1ref.m
%%%
%%% Creates a Neutral Density of the first kind using the time-mean
%%% simulation output.
%%%2

addpath ../NeutDens;
addpath ../NeutDens/matlab-interface;
addpath ../GSW;
addpath ../GSW/library/



%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = '../experiments';
% expname = 'hires_seq_onethird_RTOPO2';
% expname = 'hires_seq_onesixth_RTOPO2';
expname = 'hires_seq_onetwelfth_notides_RTOPO2';
loadexp;






%%%%%%%%%%%%%%%%%%%
%%%%% OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%

%%% Time frame over which to average thermodynamic variables to create
%%% climatology
% tmin = 18.05*86400*365;
% tmax = 27.05*86400*365;
% tmin = 9.05*86400*365;
% tmax = 18.05*86400*365;
tmin = 1.05*86400*365;
tmax = 9.05*86400*365;

%%% b-factor of Jackett and McDougall 1997
bfac = 2;

%%% Informational flags 
FLAG_UNASSIGNED = -1;
FLAG_ASSIGNED = 0;
FLAG_UNASSIGNABLE = 1;
FLAG_UNSTABLE_STRAT = 2;


















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












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZE NEUTRAL DENSITY %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% To store reference density
gg_ref = 0*pt_ref;

%%% Calculate initial neutral density profile in north-east corner
ss_init = squeeze(ss_ref(Nx,Ny,:));
pt_init = squeeze(pt_ref(Nx,Ny,:));
pp_init = squeeze(pp_ref(Nx,Ny,:));
lon_init = XC(Nx,Ny);
lat_init = YC(Nx,Ny);
ssa_init = gsw_SA_from_SP(ss_init,pp_init,lon_init,lat_init);  
ttc_init = gsw_CT_from_pt(ssa_init,pt_init);
ttis_init = gsw_t_from_CT(ssa_init,ttc_init,pp_init);    
[gg_init,dg_lo,dg_hi,wts,sc,tc,pc] = gamma_n(ss_init,ttis_init,pp_init,lon_init,lat_init);  
gg_ref(Nx,Ny,:) = reshape(gg_init,[1 1 Nr]);












%%% TODO consider restricting the comparison with assigned casts to a fixed
%%% distance

%%% TODO consider calculating true distance (at least pythagorean distance,
%%% if not true on a sphere)

%%% TODO compare resulting reference ND product against that calculated by
%%% gamma-n






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ITERATIVELY GENERATE REFERENCE DENSITY %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% To store indices of assigned casts
ncasts = length(find(sum(hFacC,3)>0));
i_assigned = zeros(ncasts,1);
j_assigned = zeros(ncasts,1);
n_assigned = 1;
cast_flags = FLAG_UNASSIGNED*ones(ncasts,1);
cast_flags(1) = FLAG_ASSIGNED;
i_assigned(1) = Nx;
j_assigned(1) = Ny;

%%% These confusing loops step across the domain, looping through
%%% gridpoints in diagnonal lines. Each line starts from the northern
%%% boundary or the eastern boundary, and moves southeastward to the
%%% eastern or southern boundary.
for i0 = Nx-1:-1:-(Ny-2)
% for i0 = Nx-1
  for i = max(i0,1):1:min(Nx,Ny-1+i0)
%   for i = Nx-1    
    j = Ny - (i-i0);
    
    %%% Track progrss
    disp(['i=',num2str(i),'j=',num2str(j)]);    
    tstart = tic();
 
    %%% Min and max indices of wet cells in this water column
    kmin = find(hFacC(i,j,:)>0,1,'first');
    kmax = find(hFacC(i,j,:)>0,1,'last');
    
    %%% Ignore dry columns
    if (isempty(kmin))
      continue;
    end

    %%% Stores the indices of the column that will be used as a reference
    %%% "cast" for each point in the current water column. This allows
    %%% distant water columns to be used as a reference if the adjacent water
    %%% column does not contain waters that are sufficiently dense/light to
    %%% fully label the cast.
    i_refcast = zeros(1,Nr);
    j_refcast = zeros(1,Nr);

    %%% In cases where no suitable water column can be found to serve as a
    %%% reference "cast", these indices record the water columns that have densities
    %%% closest to the density of the point requiring a neutral density
    %%% label.
    topdiff_min = 10000*ones(1,Nr);
    botdiff_min = 10000*ones(1,Nr);
    i_top_mindiff = zeros(1,Nr);
    i_bot_mindiff = zeros(1,Nr);
    j_top_mindiff = zeros(1,Nr);
    j_bot_mindiff = zeros(1,Nr);
    k_top_mindiff = zeros(1,Nr);
    k_bot_mindiff = zeros(1,Nr);

    %%% Stores a list of grid points that cannot be assigned neutral density
    %%% labels via vertical interpolation along any previously-assigned cast
%     k_flag_top = [];
%     k_flag_bot = [];
%     k_unassigned_top = kmax;
%     k_unassigned_bot = kmin;
    k_unassigned = ones(1,Nr);
    k_unassigned(1:kmin-1) = NaN;
    k_unassigned(kmax+1:Nr) = NaN;
    
    %%% Sort all previously-assigned but unflagged casts by distance from the current
    %%% cast
    searchidx = find(cast_flags==FLAG_ASSIGNED);
    Nsearch = length(searchidx);
    dist = sqrt((i-i_assigned(searchidx)).^2 + (j-j_assigned(searchidx)).^2);
    [unused,sortidx]=sort(dist);
    
    %%% Loop over previously-assigned casts to find reference casts within
    %%% whose density range each parcel in the current cast falls
    for n=1:Nsearch
      
      %%% Condition to stop: we have found casts from which we can assign 
      %%% neutral density labels to all points in this water column   
      if (nansum(k_unassigned)==0)
        break;
      end
      
      %%% Indices of the next-nearest cast
      i_cast = i_assigned(sortidx(n));
      j_cast = j_assigned(sortidx(n));

      %%% Vertical limits of cast data
      kmin_cast = find(hFacC(i_cast,j_cast,:)>0,1,'first');
      kmax_cast = find(hFacC(i_cast,j_cast,:)>0,1,'last');
      
      %%% Properties at the top and bottom of the candidate reference cast
      ss_cast_top = ss_ref(i_cast,j_cast,kmin_cast);
      pt_cast_top = pt_ref(i_cast,j_cast,kmin_cast);
      pp_cast_top = pp_ref(i_cast,j_cast,kmin_cast);
      ss_cast_bot = ss_ref(i_cast,j_cast,kmax_cast);
      pt_cast_bot = pt_ref(i_cast,j_cast,kmax_cast);
      pp_cast_bot = pp_ref(i_cast,j_cast,kmax_cast);
      
      %%% Max and min indices of the parcels that can be assigned a density
      %%% from the candidate reference cast
      k_max_assignable = kmin - 1;
      k_min_assignable = kmax + 1;
            
      %%% STEP 1:
      %%% Find the first parcel that is denser than the lightest point in
      %%% the candidate reference cast
      for k=kmin:kmax
        
        %%% "Parcel" properties (just the properties of this grid cell)
        ss_parcel = ss_ref(i,j,k);
        pt_parcel = pt_ref(i,j,k);
        pp_parcel = pp_ref(i,j,k);
        
        %%% Density difference between parcel and top of the candidate reference cast
        pp_mid = 0.5*(pp_parcel+pp_cast_top);
        topdiff = - densjmd95(ss_parcel,pt_parcel,pp_mid) ...
                + densjmd95(ss_cast_top,pt_cast_top,pp_mid);
           
        %%% Parcel is denser than the top of the candidate reference cast,
        %%% so record this index and move to next step.
        if (topdiff < 0)
          k_min_assignable = k;
          break;
        else
          if (topdiff < topdiff_min(k))
            topdiff_min(k) = topdiff;
            i_top_mindiff(k) = i_cast;
            j_top_mindiff(k) = j_cast;
            k_top_mindiff(k) = kmin_cast;
          end
        end
        
      end
      
      %%% STEP 2:
      %%% Find the first parcel that is lighter than the densest point in
      %%% the candidate reference cast
      for k=kmax:-1:kmin
        
        %%% "Parcel" properties (just the properties of this grid cell)
        ss_parcel = ss_ref(i,j,k);
        pt_parcel = pt_ref(i,j,k);
        pp_parcel = pp_ref(i,j,k);
        
        %%% Density difference between parcel and bottom of the candidate reference cast
        pp_mid = 0.5*(pp_parcel+pp_cast_bot);
        botdiff = densjmd95(ss_parcel,pt_parcel,pp_mid) ...
                - densjmd95(ss_cast_bot,pt_cast_bot,pp_mid);
                            
        %%% Parcel is lighter than the bottom of the candidate reference cast,
        %%% so record this index and move to next step.
        if (botdiff < 0)
          k_max_assignable = k;
          break;
        else
          if (botdiff < botdiff_min(k))
            botdiff_min(k) = botdiff;
            i_bot_mindiff(k) = i_cast;
            j_bot_mindiff(k) = j_cast;
            k_bot_mindiff(k) = kmax_cast;
          end
        end
        
        
      end      
      
      %%% STEP 3:
      %%% If there are no parcels that can be assigned using this candidate reference
      %%% cast then just move on
      if (k_min_assignable > k_max_assignable)
        continue;
      end

      %%% There are unassigned parcels that fall into the range that can be
      %%% assigned using the current candidate reference cast, so assign
      %%% them this cast as their reference cast
      kk = 1:Nr;
      i_refcast((k_unassigned==1) & (kk>=k_min_assignable) & (kk<=k_max_assignable)) = i_cast;
      j_refcast((k_unassigned==1) & (kk>=k_min_assignable) & (kk<=k_max_assignable)) = j_cast;
      k_unassigned((k_unassigned==1) & (kk>=k_min_assignable) & (kk<=k_max_assignable)) = 0;
%       if (k_min_assignable <= k_unassigned_top)
%         i_refcast(k_min_assignable:min(k_max_assignable,k_unassigned_top)) = i_cast;
%         j_refcast(k_min_assignable:min(k_max_assignable,k_unassigned_top)) = j_cast;
%         k_unassigned_top = k_min_assignable - 1;
%       end
%       if (k_max_assignable >= k_unassigned_bot)        
%         i_refcast(max(k_min_assignable,k_unassigned_bot):k_max_assignable) = i_cast;
%         j_refcast(max(k_min_assignable,k_unassigned_bot):k_max_assignable) = j_cast;
%         k_unassigned_bot = k_max_assignable + 1;          
%       end
      
    end
              
      

%     %%% Determine i_refcast and j_refcast: indices of the "cast" of data to use 
%     %%% to label the kth point in this water column. We search for the closest cast within
%     %%% whose density range the point falls.       
%     for k=kmin:kmax
% 
%       %%% Track progress
% %       disp(['k=',num2str(k)]);
%       
%       %%% "Parcel" properties (just the properties of this grid cell)
%       ss_parcel = ss_ref(i,j,k);
%       pt_parcel = pt_ref(i,j,k);
%       pp_parcel = pp_ref(i,j,k);
% 
%       %%% Search to find the closest cast whose densities span the density 
%       %%% of the kth point on the current cast
%       topdiff_max = 10000;
%       botdiff_max = 10000;
%       for n=1:n_assigned
%         
%         %%% Indices of the next-nearest cast
%         i_cast = i_assigned(sortidx(n));
%         j_cast = j_assigned(sortidx(n));
%         
%         %%% Vertical limits of cast data
%         kmin_cast = find(hFacC(i_cast,j_cast,:)>0,1,'first');
%         kmax_cast = find(hFacC(i_cast,j_cast,:)>0,1,'last');
% 
% %         %%% Version 1: using densDiff
% %         ss_cast = squeeze(ss_ref(i_cast,j_cast,:));
% %         pt_cast = squeeze(pt_ref(i_cast,j_cast,:));
% %         pp_cast = squeeze(pp_ref(i_cast,j_cast,:));               
% %         topdiff = -densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,kmin_cast,kmax_cast,pp_ref(i_cast,j_cast,kmin_cast));
% %         botdiff = densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,kmin_cast,kmax_cast,pp_ref(i_cast,j_cast,kmax_cast));
%         
%         %%% Version 2: using dens function directly        
%         pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i_cast,j_cast,kmin_cast));
%         topdiff = - densjmd95(ss_parcel,pt_parcel,pp_mid) ...
%                 + densjmd95(ss_ref(i_cast,j_cast,kmin_cast),pt_ref(i_cast,j_cast,kmin_cast),pp_mid);
%         pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i_cast,j_cast,kmax_cast));
%         botdiff = densjmd95(ss_parcel,pt_parcel,pp_mid) ...
%                 - densjmd95(ss_ref(i_cast,j_cast,kmax_cast),pt_ref(i_cast,j_cast,kmax_cast),pp_mid);
%         
%         %%% If we find a cast in which this parcel's density exceeds that of
%         %%% the lightest point and is lower than that of the densest point,
%         %%% then we define this j-index to be j_cast, and use it as a
%         %%% reference "cast" for labeling this point with neutral density.
%         if ( (topdiff < 0) && (botdiff < 0) )
%           i_refcast(k) = i_cast;
%           j_refcast(k) = j_cast;
%           i_top_mindiff(k) = i_cast;
%           j_top_mindiff(k) = j_cast;
%           k_top_mindiff(k) = kmin_cast;
%           i_bot_mindiff(k) = i_cast;
%           j_bot_mindiff(k) = j_cast;          
%           k_bot_mindiff(k) = kmax_cast;
%           break;
%         end
% 
%         %%% Finds the "cast" with the smallest density difference between
%         %%% this point and its densest/lightest point
%         if (topdiff < topdiff_max)
%           topdiff_max = topdiff;
%           i_top_mindiff(k) = i_cast;
%           j_top_mindiff(k) = j_cast;
%           k_top_mindiff(k) = kmin_cast;
%         end
%         if (botdiff < botdiff_max)
%           botdiff_max = botdiff;          
%           i_bot_mindiff(k) = i_cast;
%           j_bot_mindiff(k) = j_cast;
%           k_bot_mindiff(k) = kmax_cast;
%         end
%         
%       end %%% for n=1:n_assigned
% 
%       %%% If we couldn't find a cast within whose density range the
%       %%% parcel's density lies, we flag this point as needing
%       %%% extrapolation
%       if (j_refcast(k) == 0)    
%         
%         if (topdiff_max > 0)
%           k_flag_top = [k k_flag_top];
%         end
%         if (botdiff_max > 0)
%           k_flag_bot = [k k_flag_bot];
%         end
% 
%       end      
% 
%     end

%     %%% Loop through vertical grid cells and assign density labels
%     for k=kmin:kmax %%% Only loop over wet grid cells          
% 
%       %%% Exclude points for which no appropriate reference cast could be
%       %%% found
%       if (j_refcast(k) == 0 )
%         continue;
%       end

    %%% Loop through vertical grid cells and assign density labels
    %%% Exclude points for which no appropriate reference cast could be
    %%% found
%     for k=k_unassigned_top+1:k_unassigned_bot-1
    for k = find(k_unassigned==0)
     
      %%% "Parcel" properties (just the properties of this grid cell)
      ss_parcel = ss_ref(i,j,k);
      pt_parcel = pt_ref(i,j,k);
      pp_parcel = pp_ref(i,j,k);

      %%% "Cast" properties (just the properties of the adjacent column of
      %%% grid cells, which already have neutral density labels)
      ss_cast = ss_ref(i_refcast(k),j_refcast(k),:);
      pt_cast = pt_ref(i_refcast(k),j_refcast(k),:);
      pp_cast = pp_ref(i_refcast(k),j_refcast(k),:);      
      kmin_cast = find(hFacC(i_refcast(k),j_refcast(k),:)>0,1,'first');
      kmax_cast = find(hFacC(i_refcast(k),j_refcast(k),:)>0,1,'last');

      %%% Initial guess for neutral pressure on cast
      pp_neut = pp_ref(i,j,k);

      %%% Find the pressure on the "cast" corresponding to the point that is
      %%% neutral to the "parcel"
      [pp_neut,fval,exitflag] = fzero(@(p) densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,kmin_cast,kmax_cast,p), pp_neut);
      
      %%% If we fail to assign a neutral density then flag this parcel as
      %%% not having had a reference cast assigned to it
      if (exitflag ~= 1)
        gg_ref(i,j,k)= NaN;
        cast_flags(n_assigned+1) = FLAG_UNASSIGNABLE;
        continue;
      end
      
      %%% Label this point with neutral density. Currently using linear
      %%% interpolation vertically on "cast" data to assign label, which is
      %%% somewhat crude, though our very tightly spaced gridpoints help
      %%% here.
      gg_ref(i,j,k) = interpCast(gg_ref(i_refcast(k),j_refcast(k),:),pp_cast,kmin_cast,kmax_cast,pp_neut); 

    end

    %%% Fill in points in the water column that could not be
    %%% assign a neutral density label via linear extrapolation from the
    %%% closest-density point that was previously assigned a neutral density label
    for k = find(k_unassigned==1)
    
      %%% Calculate T,S,P,dT,dS,alpha,beta,drho all mid-way between this grid cell and the
      %%% reference grid cell
      if (topdiff_min(k) < botdiff_min(k))
        i_mindiff = i_top_mindiff(k);
        j_mindiff = j_top_mindiff(k);
        k_mindiff = k_top_mindiff(k);      
      else
        i_mindiff = i_bot_mindiff(k);
        j_mindiff = j_bot_mindiff(k);
        k_mindiff = k_bot_mindiff(k);
      end
        
      %%% Compute density difference at midpoint pressure
      pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i_mindiff,j_mindiff,k_mindiff));
      drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
           - densjmd95(ss_ref(i_mindiff,j_mindiff,k_mindiff),pt_ref(i_mindiff,j_mindiff,k_mindiff),pp_mid);
      
      %%% Assign density label via linear extrapolation    
      gg_ref(i,j,k) = gg_ref(i_mindiff,j_mindiff,k_mindiff) + drho .* bfac;
      
    end
    
    %%% Fill in any points that could not be assigned a density value
    if (cast_flags == FLAG_UNASSIGNABLE)   

      %%% Loop through cast to find missing densities
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

        %%% Block of NaNs at top of water column
        if (k_firstnan == kmin)

          %%% Extrapolate recursively upward
          for l=k_lastnan:-1:k_firstnan
            pp_mid = 0.5*(pp_ref(i,j,l+1)+pp_ref(i,j,l));
            drho = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
                - densjmd95(ss_ref(i,j,l+1),pt_ref(i,j,l+1),pp_mid);
              gg_ref(i,j,l) = gg_ref(i,j,l+1) + drho .* bfac;
          end

        %%% Block of NaNs at bottom of water column
        elseif (k_lastnan == kmax)

          %%% Extrapolate recursively downward
          for l=k_firstnan:1:k_lastnan
            pp_mid = 0.5*(pp_ref(i,j,l-1)+pp_ref(i,j,l));
            drho = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
                - densjmd95(ss_ref(i,j,l-1),pt_ref(i,j,l-1),pp_mid);
              gg_ref(i,j,l) = gg_ref(i,j,l-1) + drho .* bfac;
          end

        %%% Block of NaNs in middle of water column
        else

          %%% Extrapolate upward and downward and take the mean
          for l=k_firstnan:1:k_lastnan
            pp_mid = 0.5*(pp_ref(i,j,k_firstnan-1)+pp_ref(i,j,l));
            drho1 = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
                - densjmd95(ss_ref(i,j,k_firstnan-1),pt_ref(i,j,k_firstnan-1),pp_mid);
            pp_mid = 0.5*(pp_ref(i,j,k_lastnan+1)+pp_ref(i,j,l));
            drho2 = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
                - densjmd95(ss_ref(i,j,k_lastnan+1),pt_ref(i,j,k_lastnan+1),pp_mid);
            gg_ref(i,j,l) = 0.5 * ( gg_ref(i,j,k_firstnan-1) + drho1 .* bfac ...
                                  + gg_ref(i,j,k_lastnan+1) + drho2 .* bfac );
          end

        end 
        
        k_firstnan = [];
        k_lastnan = [];

      end
      
      cast_flags(n_assigned+1) = FLAG_UNASSIGNED;

    end

         
    
    %%% Flag casts with unstable stratification
    %%% Finally, fill in any values that could not be assigned density
    %%% values or which have produced unstable stratification
    for k=kmin+1:kmax
      
      dg = gg_ref(i,j,k)-gg_ref(i,j,k-1); 
      if (dg < 0)
         pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i,j,k-1));
         drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
            - densjmd95(ss_ref(i,j,k-1),pt_ref(i,j,k-1),pp_mid);
        if ((drho >= 0) || (dg/drho > 1.5*bfac))
          cast_flags(n_assigned+1) = FLAG_UNSTABLE_STRAT;    
        end        
      end
      
    end
    
    %%% Fix up any casts with inappropriately unstable stratification
    if (cast_flags(n_assigned+1) == FLAG_UNSTABLE_STRAT)
      
      %%% Keep searching until there are no instances of inappropriately
      %%% unstable stratification remaining
      cast_is_stable = false;
      while (~cast_is_stable)

        %%% Loop through to find unstable stratification
        k_unstable = [];
        for k=kmin+1:kmax

          %%% If the stratification is unstable, and should really be
          %%% stable, then replace this density
          dg = gg_ref(i,j,k)-gg_ref(i,j,k-1);
          if (dg < 0)
            pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i,j,k-1));
            drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
                - densjmd95(ss_ref(i,j,k-1),pt_ref(i,j,k-1),pp_mid);
            if ( (drho >= 0) || (dg/drho > 1.5*bfac) )         
              k_unstable = k;
              break;
            end
          end

        end
       
        %%% No more instances of unstable stratification -> we can all go home
        if (isempty(k_unstable))
          cast_is_stable = true;
          continue;
        end        

        %%% Otherwise, recursively replace densities to make the cast
        %%% stably stratified. We create two test casts: one by extrapolating
        %%% density upward, and one by extrapolating density downward. We favor
        %%% the test cast that incurs fewer replacements of the parcel
        %%% densities.
        gg_test_up = gg_ref(i,j,:);
        gg_test_dn = gg_ref(i,j,:);
        n_mods_up = 0;
        n_mods_dn = 0;

        %%% Loop up toward top of cast
        k = k_unstable;    
        while (k > kmin)

          k = k - 1;

          %%% Compute true density difference
          pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i,j,k+1));
          drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
                - densjmd95(ss_ref(i,j,k+1),pt_ref(i,j,k+1),pp_mid);

          %%% If stratification is inappropriately unstable then extrapolate
          dg = gg_test_up(k) - gg_test_up(k+1);
          if ((dg > 0) && ((drho <= 0) || (dg/drho > 1.5*bfac)))
            gg_test_up(k) = gg_test_up(k+1) + drho * bfac;
            n_mods_up = n_mods_up + 1;
          else
            %%% Stop iterating if stratification is no longer inappropriately
            %%% unstable
            break;
          end

        end 
        kstop_up = k;
        diff_up = gg_test_up(k_unstable) - mean(gg_test_up(kstop_up:k_unstable-1));    

        %%% Loop down toward bottom of cast
        k = k_unstable-1;    
        while (k < kmax)      

          k = k + 1;

          %%% Compute true density difference
          pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i,j,k-1));
          drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
                - densjmd95(ss_ref(i,j,k-1),pt_ref(i,j,k-1),pp_mid);

          %%% If stratification is inappropriately unstable then extrapolate
          dg = gg_test_dn(k) - gg_test_dn(k-1);
          if ((dg < 0) && ((drho >= 0) || (dg/drho > 1.5*bfac)))
            gg_test_dn(k) = gg_test_dn(k-1) + drho * bfac;
            n_mods_dn = n_mods_dn + 1;
          else
            %%% Stop iterating if stratification is no longer inappropriately
            %%% unstable
            break;
          end

        end    
        kstop_dn = k;
        diff_dn = mean(gg_test_dn(k_unstable:kstop_dn)) - gg_test_dn(k_unstable-1);
        
        %%% Replace this entire cast with whichever modified cast incurs fewer
        %%% density changes
        if ( ((kstop_up == kmin) && (kstop_dn == kmax)) ...
          || ((kstop_up > kmin) && (kstop_dn < kmax)))
           if (n_mods_up < n_mods_dn)    
            gg_ref(i,j,:) = gg_test_up;
          else
            gg_ref(i,j,:) = gg_test_dn;
          end
        else

    %       if (nanstd(gg_test_up)<nanstd(gg_test_dn))
      %     if (diff_up < diff_dn)
          if (kstop_dn == kmax)
            gg_ref(i,j,:) = gg_test_up;
          else
            gg_ref(i,j,:) = gg_test_dn;
          end
        end

      end
      
      cast_flags(n_assigned+1) = FLAG_UNASSIGNED;
  
    end
    
    
    
    
    
    
    %%% Flag this cast as having been assigned
    if (cast_flags(n_assigned+1) == FLAG_UNASSIGNED)
      cast_flags(n_assigned+1) = FLAG_ASSIGNED;      
    end
      
      
    
      
    
%     %%% Fill in points at the bottom of the water column that were too dense to
%     %%% assign a neutral density label via linear extrapolation from the
%     %%% densest point that was previously assigned a neutral density label
%     for k = kmin:k_unassigned_top
% 
%       %%% Calculate T,S,P,dT,dS,alpha,beta,drho all mid-way between this grid cell and the
%       %%% reference grid cell
%       i_mindiff = i_top_mindiff(k);
%       j_mindiff = j_top_mindiff(k);
%       k_mindiff = k_top_mindiff(k);      
% 
%       %%% Compute density difference at midpoint pressure
%       pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i_mindiff,j_mindiff,k_mindiff));
%       drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
%            - densjmd95(ss_ref(i_mindiff,j_mindiff,k_mindiff),pt_ref(i_mindiff,j_mindiff,k_mindiff),pp_mid);
% 
%       %%% Assign density label via linear extrapolation    
%       gg_ref(i,j,k) = gg_ref(i_mindiff,j_mindiff,k_mindiff) + drho .* bfac;
% 
%     end 
% 
%     %%% Fill in points at the bottom of the water column that were too dense to
%     %%% assign a neutral density label via linear extrapolation from the
%     %%% densest point that was previously assigned a neutral density label
%     for k = k_unassigned_bot:kmax
% 
%       %%% Calculate T,S,P,dT,dS,alpha,beta,drho all mid-way between this grid cell and the
%       %%% reference grid cell
%       i_mindiff = i_bot_mindiff(k);
%       j_mindiff = j_bot_mindiff(k);
%       k_mindiff = k_bot_mindiff(k);      
% 
%       %%% Compute density difference at midpoint pressure
%       pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i_mindiff,j_mindiff,k_mindiff));
%       drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
%            - densjmd95(ss_ref(i_mindiff,j_mindiff,k_mindiff),pt_ref(i_mindiff,j_mindiff,k_mindiff),pp_mid);
% 
%       %%% Assign density label via linear extrapolation    
%       gg_ref(i,j,k) = gg_ref(i_mindiff,j_mindiff,k_mindiff) + drho .* bfac;        
% 
%     end  
    
%     %%% Fill in points at the bottom of the water column that were too dense to
%     %%% assign a neutral density label via linear extrapolation from the
%     %%% densest point that was previously assigned a neutral density label
%     for k = k_flag_top
% 
%       %%% Calculate T,S,P,dT,dS,alpha,beta,drho all mid-way between this grid cell and the
%       %%% reference grid cell
%       i_mindiff = i_top_mindiff(k);
%       j_mindiff = j_top_mindiff(k);
%       k_mindiff = k_top_mindiff(k);
%       
% %       %%% Version 1: following Stewart and Thompson (2015)
% %       pt_mid = 0.5*(pt_ref(i,j,k)+pt_ref(i_mindiff,j_mindiff,k_mindiff));
% %       ss_mid = 0.5*(ss_ref(i,j,k)+ss_ref(i_mindiff,j_mindiff,k_mindiff));
% %       pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i_mindiff,j_mindiff,k_mindiff));    
% %       ds_mid = ss_ref(i,j,k)-ss_ref(i_mindiff,j_mindiff,k_mindiff);
% %       dt_mid = pt_ref(i,j,k)-pt_ref(i_mindiff,j_mindiff,k_mindiff);
% %       [alpha_m,beta_m] = calcAlphaBeta(ss_mid,pt_mid,pp_mid);
% %       drho_mid = (beta_m.*ds_mid - alpha_m.*dt_mid);
% %       dens_mid = densjmd95(ss_mid,pt_mid,pp_mid);
% %       drho = dens_mid.*drho_mid;     
%       
%       %%% Version 2: following Jackett and McDougall (1997)
%       pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i_mindiff,j_mindiff,k_mindiff));
%       drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
%            - densjmd95(ss_ref(i_mindiff,j_mindiff,k_mindiff),pt_ref(i_mindiff,j_mindiff,k_mindiff),pp_mid);
%       
%       %%% Assign density label via linear extrapolation    
%       gg_ref(i,j,k) = gg_ref(i_mindiff,j_mindiff,k_mindiff) + drho .* bfac;
% 
%     end 
% 
%     %%% Fill in points at the bottom of the water column that were too dense to
%     %%% assign a neutral density label via linear extrapolation from the
%     %%% densest point that was previously assigned a neutral density label
%     for k = k_flag_bot
% 
%       %%% Calculate T,S,P,dT,dS,alpha,beta,drho all mid-way between this grid cell and the
%       %%% reference grid cell
%       i_mindiff = i_bot_mindiff(k);
%       j_mindiff = j_bot_mindiff(k);
%       k_mindiff = k_bot_mindiff(k);      
%       
% %       %%% Version 1: following Stewart and Thompson (2015)
% %       pt_mid = 0.5*(pt_ref(i,j,k)+pt_ref(i_mindiff,j_mindiff,k_mindiff));
% %       ss_mid = 0.5*(ss_ref(i,j,k)+ss_ref(i_mindiff,j_mindiff,k_mindiff));
% %       pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i_mindiff,j_mindiff,k_mindiff));    
% %       ds_mid = ss_ref(i,j,k)-ss_ref(i_mindiff,j_mindiff,k_mindiff);
% %       dt_mid = pt_ref(i,j,k)-pt_ref(i_mindiff,j_mindiff,k_mindiff);
% %       [alpha_m,beta_m,dalpha_dT_m,dalpha_dS_m,dalpha_dz_m, ...
% %               dbeta_dT_m,dbeta_dS_m,dbeta_dz_m] = calcAlphaBeta(ss_mid,pt_mid,pp_mid);
% %       drho_mid = (beta_m.*ds_mid - alpha_m.*dt_mid);
% %       dens_mid = densjmd95(ss_mid,pt_mid,pp_mid);
% %       drho = drho_mid .* dens_mid;
% 
%       %%% Version 2: following Jackett and McDougall (1997)
%       pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i_mindiff,j_mindiff,k_mindiff));
%       drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
%            - densjmd95(ss_ref(i_mindiff,j_mindiff,k_mindiff),pt_ref(i_mindiff,j_mindiff,k_mindiff),pp_mid);
%          
%       %%% Assign density label via linear extrapolation    
%       gg_ref(i,j,k) = gg_ref(i_mindiff,j_mindiff,k_mindiff) + drho .* bfac;        
% 
%     end      
    
    %%% Increment list of assigned casts
    n_assigned = n_assigned + 1;
    i_assigned(n_assigned) = i;
    j_assigned(n_assigned) = j;
    
    %%% Track progress
    disp(toc(tstart));
   
  end

end








% %%% Fill in any points that could not be assigned a density value
% searchidx = find(cast_flags == FLAG_UNASSIGNABLE);
% for n = 1:length(searchidx)
%   
%   %%% Horizontal grid indices
%   i = i_assigned(searchidx(n));
%   j = j_assigned(searchidx(n));
%   
%   %%% Min and max indices of wet cells in this water column
%   kmin = find(hFacC(i,j,:)>0,1,'first');
%   kmax = find(hFacC(i,j,:)>0,1,'last');
% 
%   %%% Loop through cast to find missing densities
%   k_firstnan = [];
%   k_lastnan = [];
%   for k=kmin:kmax
% 
%     %%% NaNs indicate that the search for a matching density level in
%     %%% this parcel's reference cast failed
%     if (isnan(gg_ref(i,j,k)))
%       if (isempty(k_firstnan))     
%         k_firstnan = k;
%       end
%       if (k~=kmax)
%         continue;
%       else
%         k_lastnan = k;
%       end
%     else         
%       if (isempty(k_firstnan))
%         continue;
%       else     
%         k_lastnan = k-1;
%       end             
%     end
%     
%     %%% Block of NaNs at top of water column
%     if (k_firstnan == kmin)
%       
%       %%% Extrapolate recursively upward
%       for l=k_lastnan:-1:k_firstnan
%         pp_mid = 0.5*(pp_ref(i,j,l+1)+pp_ref(i,j,l));
%         drho = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
%             - densjmd95(ss_ref(i,j,l+1),pt_ref(i,j,l+1),pp_mid);
%           gg_ref(i,j,l) = gg_ref(i,j,l+1) + drho .* bfac;
%       end
%       
%     %%% Block of NaNs at bottom of water column
%     elseif (k_lastnan == kmax)
%      
%       %%% Extrapolate recursively downward
%       for l=k_firstnan:1:k_lastnan
%         pp_mid = 0.5*(pp_ref(i,j,l-1)+pp_ref(i,j,l));
%         drho = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
%             - densjmd95(ss_ref(i,j,l-1),pt_ref(i,j,l-1),pp_mid);
%           gg_ref(i,j,l) = gg_ref(i,j,l-1) + drho .* bfac;
%       end
%     
%     %%% Block of NaNs in middle of water column
%     else
%       
%       %%% Extrapolate upward and downward and take the mean
%       for l=k_firstnan:1:k_lastnan
%         pp_mid = 0.5*(pp_ref(i,j,k_firstnan-1)+pp_ref(i,j,l));
%         drho1 = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
%             - densjmd95(ss_ref(i,j,k_firstnan-1),pt_ref(i,j,k_firstnan-1),pp_mid);
%         pp_mid = 0.5*(pp_ref(i,j,k_lastnan+1)+pp_ref(i,j,l));
%         drho2 = densjmd95(ss_ref(i,j,l),pt_ref(i,j,l),pp_mid) ...
%             - densjmd95(ss_ref(i,j,k_lastnan+1),pt_ref(i,j,k_lastnan+1),pp_mid);
%         gg_ref(i,j,l) = 0.5 * ( gg_ref(i,j,firstnan-1) + drho1 .* bfac ...
%                               + gg_ref(i,j,lastnan+1) + drho2 .* bfac );
%       end
%       
%     end     
% 
%   end
% 
% end






% %%% Fix up any casts with inappropriately unstable stratification
% searchidx = find(cast_flags == FLAG_UNSTABLE_STRAT);
% for n = 1:length(searchidx)
%   n
%   %%% Horizontal grid indices
%   i = i_assigned(searchidx(n));
%   j = j_assigned(searchidx(n));
%   
%   %%% Min and max indices of wet cells in this water column
%   kmin = find(hFacC(i,j,:)>0,1,'first');
%   kmax = find(hFacC(i,j,:)>0,1,'last');
% 
%   %%% Keep searching until there are no instances of inappropriately
%   %%% unstable stratification remaining
%   cast_is_stable = false;
%   while (~cast_is_stable)
%     
%     %%% Loop through to find unstable stratification
%     k_unstable = [];
%     for k=kmin+1:kmax
% 
%       %%% If the stratification is unstable, and should really be
%       %%% stable, then replace this density
%       if (gg_ref(i,j,k)<gg_ref(i,j,k-1))
%         pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i,j,k-1));
%         drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
%             - densjmd95(ss_ref(i,j,k-1),pt_ref(i,j,k-1),pp_mid);
%         if (drho > 0)         
%           k_unstable = k;
%           break;
%         end
%       end
%       
%     end
%     
%     %%% No more instances of unstable stratification -> we can all go home
%     if (isempty(k_unstable))
%       cast_is_stable = true;
%       continue;
%     end        
%     
%     %%% Otherwise, recursively replace densities to make the cast
%     %%% stably stratified. We create two test casts: one by extrapolating
%     %%% density upward, and one by extrapolating density downward. We favor
%     %%% the test cast that incurs fewer replacements of the parcel
%     %%% densities.
%     gg_test_up = gg_ref(i,j,:);
%     gg_test_dn = gg_ref(i,j,:);
%     n_mods_up = 0;
%     n_mods_dn = 0;
%     
%     %%% Loop up toward top of cast
%     k = k_unstable;    
%     while (k > kmin)
%             
%       k = k - 1;
%       
%       %%% Compute true density difference
%       pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i,j,k+1));
%       drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
%             - densjmd95(ss_ref(i,j,k+1),pt_ref(i,j,k+1),pp_mid);
%           
%       %%% If stratification is inappropriately unstable then extrapolate
%       if ((gg_test_up(k) > gg_test_up(k+1)) && (drho < 0))
%         gg_test_up(k) = gg_test_up(k+1) + drho * bfac;
%         n_mods_up = n_mods_up + 1;
%       else
%         %%% Stop iterating if stratification is no longer inappropriately
%         %%% unstable
%         break;
%       end
%            
%     end 
%     kstop_up = k;
%     diff_up = gg_test_up(k_unstable) - mean(gg_test_up(kstop_up:k_unstable-1));    
%     
%     %%% Loop down toward bottom of cast
%     k = k_unstable-1;    
%     while (k < kmax)      
%       
%       k = k + 1;
%       
%       %%% Compute true density difference
%       pp_mid = 0.5*(pp_ref(i,j,k)+pp_ref(i,j,k-1));
%       drho = densjmd95(ss_ref(i,j,k),pt_ref(i,j,k),pp_mid) ...
%             - densjmd95(ss_ref(i,j,k-1),pt_ref(i,j,k-1),pp_mid);
%           
%       %%% If stratification is inappropriately unstable then extrapolate
%       if ((gg_test_dn(k) < gg_test_dn(k-1)) && (drho > 0))
%         gg_test_dn(k) = gg_test_dn(k-1) + drho * bfac;
%         n_mods_dn = n_mods_dn + 1;
%       else
%         %%% Stop iterating if stratification is no longer inappropriately
%         %%% unstable
%         break;
%       end
%             
%     end    
%     kstop_dn = k;
%     diff_dn = mean(gg_test_dn(k_unstable:kstop_dn)) - gg_test_dn(k_unstable-1);
%     
%     %%% Replace this entire cast with whichever modified cast incurs fewer
%     %%% density changes
%     if ( ((kstop_up == kmin) && (kstop_dn == kmax)) ...
%       || ((kstop_up > kmin) && (kstop_dn < kmax)))
%        if (n_mods_up < n_mods_dn)    
%         gg_ref(i,j,:) = gg_test_up;
%       else
%         gg_ref(i,j,:) = gg_test_dn;
%       end
%     else
%       
% %       if (nanstd(gg_test_up)<nanstd(gg_test_dn))
%   %     if (diff_up < diff_dn)
%       if (kstop_dn == kmax)
%         gg_ref(i,j,:) = gg_test_up;
%       else
%         gg_ref(i,j,:) = gg_test_dn;
%       end
%     end
%   
%   end
% 
% end






%%% Finally, save the new reference dataset
save(fullfile('products',[expname,'_ND1.mat']),'ss_ref','pt_ref','pp_ref','gg_ref');






