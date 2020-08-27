%%%
%%% topog.m
%%%
%%% Defines topography and ice shelf draft for MITgcm_WS.
%%% 



%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%% SET UP %%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%    

% Getting our grid
run defineGrid.m

%%% Gaussian kernel smoothing width (in grid points)
do_smooth = true;
smooth_ice_face = false;
gauss_width = 4;

%%% Minimum water column thickness
min_thickness = 50;

%%% Minimum ice draft to qualify as an ice shelf
min_ice_thickness = 10;

%%% Set true to remove coastal embayment in SW Weddell Sea, which is prone
%%% to spurious build-up of ice
remove_embayment = false;
embayment_lat = -74; %%% Max latitude to search for the SW embayment. Original was -74.2.




%%%% Setting Directories

datadir = fullfile(gendir,'/MITgcm_WS/data/RTOPO/');

addpath ../newexp_utils




%%% Reading in grid from RTOPO dataset
% RTOPO_latitude = ncread(fullfile(datadir,'RTOPO.nc'),'lat');
% RTOPO_longitude = ncread(fullfile(datadir,'RTOPO.nc'),'lon');
% RTOPO_bathymetry = ncread(fullfile(datadir,'RTOPO.nc'),'bathy');
% RTOPO_ice = ncread(fullfile(datadir,'RTOPO.nc'),'draft');
% RTOPO_elev = ncread(fullfile(datadir,'RTOPO.nc'),'height');
RTOPO_latitude = ncread(fullfile(datadir,'RTOPO2.nc'),'lat');
RTOPO_longitude = ncread(fullfile(datadir,'RTOPO2.nc'),'lon');
RTOPO_bathymetry = ncread(fullfile(datadir,'RTOPO2.nc'),'bedrock_topography');
RTOPO_ice = ncread(fullfile(datadir,'RTOPO2.nc'),'ice_base_topography');
RTOPO_elev = ncread(fullfile(datadir,'RTOPO2.nc'),'surface_elevation');

%%% Remove superfluous additional longitude value
RTOPO_longitude = RTOPO_longitude(1:end-1);
RTOPO_bathymetry = RTOPO_bathymetry(1:end-1,:);
RTOPO_ice = RTOPO_ice(1:end-1,:);
RTOPO_elev = RTOPO_elev(1:end-1,:);

%%% Do the twist
% RTOPO_bathymetry = Tweddell(RTOPO_longitude,RTOPO_latitude,RTOPO_bathymetry);
% RTOPO_ice = Tweddell(RTOPO_longitude,RTOPO_latitude,RTOPO_ice);
% RTOPO_elev = Tweddell(RTOPO_longitude,RTOPO_latitude,RTOPO_elev);


%Scaling bathymetry to fit area of interest as defined by scaled lat/lon
lonidx = find(RTOPO_longitude >= xmin & RTOPO_longitude <= xmax);
latidx = find(RTOPO_latitude <= ymax & RTOPO_latitude >= ymin);
newbathy = RTOPO_bathymetry(lonidx,latidx);
Model_draft = RTOPO_ice(lonidx,latidx);
Model_elev = RTOPO_elev(lonidx,latidx);


%Scaling lat/lon to fit area of interest

RTOPO_latitude=RTOPO_latitude(RTOPO_latitude <= ymax & RTOPO_latitude >= ymin);
RTOPO_longitude=RTOPO_longitude(RTOPO_longitude >= xmin & RTOPO_longitude <= xmax); 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MESHGRID FOR RTOPO DATA %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LO,LA] = meshgrid(double(RTOPO_longitude),double(RTOPO_latitude));






%%%
%%% Interpolate to model grid
%%%

%%% NOTE: SPECIFIC TO CURRENT MODEL CONFIG

%%% Nearest-neighbor interpolation in order to avoid creating spurious
%%% water column thicknesses due to interpolation

MDepth=interp2(LO,LA,double(newbathy)',XMC,YMC,'nearest');
Mice_draft = interp2(LO,LA,double(Model_draft)',XMC,YMC,'nearest');
Mice_elev = interp2(LO,LA,double(Model_elev)',XMC,YMC,'nearest');
Mthickness = Mice_draft-MDepth;
















%%%
%%% Topographic smoothing (new version)
%%%

%%% Grids for smoothing
% XMC_smooth = repmat(xmc,[1 Ny]);
% YMC_smooth = repmat(ymc,[Nx 1]);
% DXG_smooth = repmat(dmxg,[1 Ny]);
% DYG_smooth = repmat(dmyg,[Nx 1]);






%%% Smooth sea floor
if (do_smooth)
  
  
  %%% Grids for smoothing
  XMC_smooth = repmat((1:Nx)',[1 Ny]);
  YMC_smooth = repmat((1:Ny),[Nx 1]);
  DXG_smooth = ones(Nx,Ny);
  DYG_smooth = ones(Nx,Ny);

  Depth_smooth = smooth2D(XMC_smooth,YMC_smooth,DXG_smooth,DYG_smooth,MDepth',gauss_width);

  %%% Smooth water column thickness where ice is present and thickness is
  %%% nonzero
  % % Mthickness_tmp = Mthickness';
  % Mthickness_tmp = Mice_draft' - Depth_smooth;
  % Mthickness_tmp(Mthickness'==0) = NaN;
  % Mthickness_tmp(Mice_draft'==0) = NaN;
  % Mthickness_smooth = smooth2D(XMC_smooth,YMC_smooth,DXG_smooth,DYG_smooth,Mthickness_tmp,gauss_width);
  % Mthickness_smooth(Mthickness'==0) = 0;
  % Mthickness_smooth(Mice_draft'==0) = -Depth_smooth(Mice_draft'==0);
  % 
  % %%% Reconstruct ice draft using smooth thickness and bathymetry
  % Smooth_Mice = Depth_smooth + Mthickness_smooth;

  %%% Smooth ice draft
  Mice_draft_tmp = Mice_draft';
  if (~smooth_ice_face)
    Mice_draft_tmp(Mice_draft'==0) = NaN;
    Mice_draft_tmp(Mthickness'==0) = NaN;
  end
  Smooth_Mice = smooth2D(XMC_smooth,YMC_smooth,DXG_smooth,DYG_smooth,Mice_draft_tmp,gauss_width);
  Smooth_Mice(Mthickness'==0) = Depth_smooth(Mthickness'==0);
  if (smooth_ice_face)
    Smooth_Mice(Smooth_Mice>-min_ice_thickness) = 0;
  else
    Smooth_Mice(Mice_draft'==0) = Mice_draft(Mice_draft==0)';
  end

else
  
  Smooth_Mice = Mice_draft';
  Depth_smooth = MDepth';
  
end
























%%%%%%%%%%%%%%%%
%%% Remove WAP%%
%%%%%%%%%%%%%%%%

WAP_lat_min = -74;
WAP_lat_mid = -67.5;
WAP_lon_mid = -66;
WAP_lon_top = -60;
WAP_lon = zeros(1,Ny);
for j=1:Ny
  
  if (ymc(j) < WAP_lat_min)
    WAP_lon(j) = xmin;
    continue;
  end
  
  if (ymc(j) < WAP_lat_mid)
    WAP_lon(j) = WAP_lon_mid;
    continue;
  end  
  
  WAP_lon(j) = WAP_lon_mid + (ymc(j)-WAP_lat_mid)*(WAP_lon_top-WAP_lon_mid)/(ymax-WAP_lat_mid);
  
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% REMOVE WAP
for i=1:Nx
  for j=1:Ny
    if (XMC(j,i)<WAP_lon(j))
      Depth_smooth(i,j) = 0;
      Smooth_Mice(i,j) = 0;
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOPOGRAPHIC SMOOTHING (OBSOLETE CODE) %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% XMC_smooth = repmat((1:Nx),[Ny 1]);
% YMC_smooth = repmat((1:Ny)',[1 Nx]);
% 
% 
% if (smooth_depth_only)
%   
%   %%% Apply smoothing filter
%   if (gauss_width > 0)
%     Depth_smooth = smooth2D(XMC_smooth',YMC_smooth',MDepth',gauss_width,gauss_width)';    
%   else
%     Depth_smooth = MDepth;
%   end
%     
%   %%% No smoothing under ice shelves  
%   Depth_smooth(Mice_draft<0 | MDepth>=0) = MDepth(Mice_draft<0 | MDepth>=0);  
%   Depth_smooth = Depth_smooth';  
%   Smooth_Mice = Mice_draft';
%   
% else
% 
%   %%% VERSION 1 - SMOOTH ICE AND BATHYMETRY SEPARATELY %%%
%   
% %   %%% Smoothing bathymetry
% %   if (gauss_width > 0)
% %     Depth_smooth = smooth2D(XMC_smooth',YMC_smooth',MDepth',gauss_width,gauss_width)';
% %     Depth_smooth = Depth_smooth';
% %   else
% %     Depth_smooth = MDepth';
% %   end
% %   
% %   %%%%%%% Smoothing Ice Draft
% %   if (gauss_width_shelfice > 0) 
% %     Smooth_Mice = smooth2D(XMC_smooth',YMC_smooth',Mice_draft',gauss_width_shelfice,gauss_width_shelfice)';
% %     Smooth_Mice = Smooth_Mice';
% %   else
% %     Smooth_Mice = Mice_draft';
% %   end
% 
% 
% %%% VERSION 2 - SMOOTH BATHYMETRY AND WATER COLUMN THICKNESS %%%
% 
% %   %%% Smoothing bathymetry
% %   if (gauss_width > 0)
% %     Depth_smooth = smooth2D(XMC_smooth',YMC_smooth',MDepth',gauss_width,gauss_width)';
% %     Depth_smooth = Depth_smooth';
% %   else
% %     Depth_smooth = MDepth';
% %   end
% %   
% %   %%% Thickness of water column between smoothed sea floor and unsmoothed
% %   %%% ice draft
% %   thickness = Mice_draft - Depth_smooth';
% %   
% %   %%%%%%% Smoothing water column thickness
% %   if (gauss_width_shelfice > 0) 
% %     thickness_smooth = smooth2D(XMC_smooth',YMC_smooth',thickness',gauss_width_shelfice,gauss_width_shelfice)';
% %     thickness_smooth = thickness_smooth';
% %   else
% %     thickness_smooth = thickness_smooth';
% %   end
% %   
% %   %%% Now reconstruct ice shelf draft from smoothed water column thickness
% %   Smooth_Mice = Depth_smooth + thickness_smooth;
% 
%   
%   %%% VERSION 3 - SMOOTH ICE AND WATER COLUMN THICKNESS %%%
%   
%   %%% First, remove all ice draft above sea level to remove errors where
%   %%% ocean meets land with ice, but there is no ice shelf
%   Mice_draft_noland = Mice_draft;
%   Mice_draft_noland(Mice_draft_noland > 0) = 0;
% 
%   %%%%%%% Smoothing Ice Draft
%   if (gauss_width_shelfice > 0) 
%     Smooth_Mice = smooth2D(XMC_smooth',YMC_smooth',Mice_draft_noland',gauss_width_shelfice,gauss_width_shelfice)';
%     Smooth_Mice = Smooth_Mice';
%   else
%     Smooth_Mice = Mice_draft_noland';
%   end
%     
%   %%% Eliminate points with very thin ice shelves
%   Smooth_Mice(Smooth_Mice>=-min_ice_thickness) = 0;
%   
%   %%% Thickness of water column between smoothed ice shelf and unsmoothed
%   %%% sea floor
%   thickness = Smooth_Mice' - MDepth;
%   
%   %%%%%%% Smoothing water column thickness
%   if (gauss_width_shelfice > 0) 
%     thickness_smooth = smooth2D(XMC_smooth',YMC_smooth',thickness',gauss_width_shelfice,gauss_width_shelfice)';
%     thickness_smooth = thickness_smooth';
%   else
%     thickness_smooth = thickness';
%   end  
%   thickness_smooth(thickness_smooth<0) = 0; %%% Can't have negative thickness
%   
%   %%% Now reconstruct bathymetry from smoothed water column thickness
%   Depth_smooth = Smooth_Mice - thickness_smooth;
% 
% end
% 
% 
% %%%%%% Anywhere that the original ice draft was above land, set equal to
% %%%%%% zero in smoothed version.
% % Smooth_Mice(Mice_draft'>=0) = 0;
%   
% %%%%%% Remove points where depth or ice draft lie above the surface
% Depth_smooth(Depth_smooth>0) = 0;
% Smooth_Mice(Smooth_Mice>0) = 0;
% 
% %%%%%% Smoothed Ice draft and Smoothed bathymetry match
% %%%%%% where the original draft minus the bathymetry equals zero.
% Smooth_Mice((Mice_draft-MDepth)'==0) = Depth_smooth((Mice_draft-MDepth)'==0);
% 
% %%% Avoid negative water column thicknesses
% Smooth_Mice(Smooth_Mice-Depth_smooth < 0) = Depth_smooth(Smooth_Mice-Depth_smooth < 0);






%%%
%%% Code to remove embayment in SW corner
%%%
if (remove_embayment)

  %%% Here we find the easternmost extent of 
  %%% the shelf ice at this latitude. We'll only add ice to points southwest
  %%% of this location.
  jmax = min(find(ymc>embayment_lat));
  imax = min(find((Depth_smooth(:,jmax)<0) & (Smooth_Mice(:,jmax)==0)))-1;

  %%% Make a copy (Smooth_Mice_tmp), then we'll alter Smooth_Mice itself
  Smooth_Mice_tmp = Smooth_Mice;
  Smooth_Mice_tmp(Smooth_Mice==0) = NaN;

  %%% Search over points southwest of (imax,jmax)
  for i=1:imax
    for j=1:jmax

      %%% Infill points with no shelf ice and nonzero ocean depth
      if ((Depth_smooth(i,j)<0) && (Smooth_Mice(i,j)==0))
                
        %%% Linear interpolation in latitude
        jnext = j-1+min(find(Smooth_Mice_tmp(i,j:end)<0));
        jprev = max(find(Smooth_Mice_tmp(i,1:j)<0));
        wnext = (ymc(j)-ymc(jprev))/(ymc(jnext)-ymc(jprev));
        wprev = 1-wnext;
        Smooth_Mice(i,j) = wnext*Smooth_Mice_tmp(i,jnext) + wprev*Smooth_Mice_tmp(i,jprev);

      end
    end
  end

end
      
      
 
      
%%%% Topographic mask: used to mask out topography. Only really useful for
%%%% making plots, as MITgcm doesn't care if dry grid cells are assigned
% %%%% non-zero values
% topog_msk = ones(Ny,Nx,Nz);
% for i=1:Nx
%     for j=1:Ny        
%         for k=1:Nz
%            if ( ( (zz(k) - dz(k)/2 > Smooth_Mice(j,i)) || ...
%                   (zz(k) + dz(k)/2 < Depth_smooth(j,i)) ) ...
%                || (Depth_smooth(j,i) == Smooth_Mice(j,i)) ) 
%                topog_msk(j,i,k) = 0;
%            end
%         end
%     end
% end


%%% If we have created negative depths then remove them

Smooth_Mice(Smooth_Mice-Depth_smooth < 0) = Depth_smooth(Smooth_Mice-Depth_smooth < 0);






%%% Calculate hFacC based on MITgcm algorithm
hFacC = zeros(Nx,Ny,Nz);
for k=1:Nz  
  hFacMnSz = max( hFacMin, min(hFacMinDr/dz(k),1) );
  for j=1:Ny        
    for i=1:Nx
 
%   o Non-dimensional distance between grid bound. and domain lower_R bound.
      hFacCtmp = (min(zzf(k),Smooth_Mice(i,j))-max(zzf(k+1),Depth_smooth(i,j))) / dz(k);
%       hFacCtmp = (max(zzf(k+1),thickness_smooth(i,j))) / dz(k);

%   o Select between, closed, open or partial (0,1,0-1)
      hFacCtmp = min( max( hFacCtmp, 0 ) , 1 );
%   o Impose minimum fraction and/or size (dimensional)
      if ( hFacCtmp < hFacMnSz ) 
        if ( hFacCtmp < hFacMnSz*0.5) 
          hFacC(i,j,k) = 0;
        else
          hFacC(i,j,k) = hFacMnSz;
        end
      else
        hFacC(i,j,k) = hFacCtmp;
      end
    end
  end
end






%%%%% Remove isolated wet grid cells at the ocean bed

Depth_smooth_old = Depth_smooth;
Smooth_Mice_old = Smooth_Mice;

found_cells = true;
while (found_cells)
  found_cells = false;
  for i=2:Nx-1
    for j=2:Ny-1
      idx = find(hFacC(i,j,:)>0);      
      kidx = max(idx);
      
      if (isempty(kidx))
        continue;
      end
            
      wet_ip1 = (hFacC(i+1,j,kidx)>0);
      wet_im1 = (hFacC(i-1,j,kidx)>0);
      wet_jp1 = (hFacC(i,j+1,kidx)>0);
      wet_jm1 = (hFacC(i,j-1,kidx)>0);
      make_dry = true;
      if (wet_ip1+wet_im1+wet_jp1+wet_jm1 >= 2)        
        if (wet_ip1 && wet_jp1 && (hFacC(i+1,j+1,kidx)>0))
          make_dry = false;
        end
        if (wet_im1 && wet_jp1 && (hFacC(i-1,j+1,kidx)>0))
          make_dry = false;
        end
        if (wet_ip1 && wet_jm1 && (hFacC(i+1,j-1,kidx)>0))
          make_dry = false;
        end
        if (wet_im1 && wet_jm1 && (hFacC(i-1,j-1,kidx)>0))
          make_dry = false;
        end
      end        
      
      if (make_dry)
        Depth_smooth(i,j) = zzf(kidx);
        hFacC(i,j,kidx) = 0;
        found_cells = true;
      end
      
    end
  end
  
end





%%%%% Remove isolated wet grid cells at ice shelf base

found_cells = true;
while (found_cells)
  found_cells = false;
  for i=2:Nx-1
    for j=2:Ny-1
      idx = find(hFacC(i,j,:)>0);      
      kidx = min(idx);
      
      if (isempty(kidx))
        continue;
      end
      %%%% defining cells (i+1,i-1,j+1,j-1, which are 'neighbors')      
      wet_ip1 = (hFacC(i+1,j,kidx)>0);
      wet_im1 = (hFacC(i-1,j,kidx)>0);
      wet_jp1 = (hFacC(i,j+1,kidx)>0);
      wet_jm1 = (hFacC(i,j-1,kidx)>0);
      make_dry = true;
      if (wet_ip1+wet_im1+wet_jp1+wet_jm1 >= 2)    %%%sum of neighboring HFaCcs    
        if (wet_ip1 && wet_jp1 && (hFacC(i+1,j+1,kidx)>0))
          make_dry = false;
        end
        if (wet_im1 && wet_jp1 && (hFacC(i-1,j+1,kidx)>0))
          make_dry = false;
        end
        if (wet_ip1 && wet_jm1 && (hFacC(i+1,j-1,kidx)>0))
          make_dry = false;
        end
        if (wet_im1 && wet_jm1 && (hFacC(i-1,j-1,kidx)>0))
          make_dry = false;
        end
      end        
      
      if (make_dry)        
        Smooth_Mice(i,j) = zzf(kidx+1);
        hFacC(i,j,kidx) = 0;
        found_cells = true;
      end
      
    end
  end
  
end

%%% If we have created negative depths then remove them

Smooth_Mice(Smooth_Mice-Depth_smooth < 0) = Depth_smooth(Smooth_Mice-Depth_smooth < 0);
















%%%
%%% Eliminate areas below minimum water column thickness
%%%

thickness_smooth = Smooth_Mice - Depth_smooth;
Depth_smooth(thickness_smooth < min_thickness) = Smooth_Mice(thickness_smooth < min_thickness);











% 
% %%%%%%%%%%%  ShElFiCe - DepTh %%%%%%%%%%%%%%
% %%%% can get a handle of true depth
% 
% thickness_smooth = Smooth_Mice - Depth_smooth;
% 
% thickness = Mice_draft - MDepth;
% 
% DZ = repmat(reshape(dz,[1 1 Nz]),[Nx,Ny,1]);
% thickness_true = sum(hFacC(:,:,:).*DZ,3);





%%% Write to files

addpath ../newexp_utils

data = Depth_smooth;
writeDataset(data,fullfile(inputconfigdir,bathyFile),ieee,prec);
clear data

data = Smooth_Mice;
writeDataset(data,fullfile(inputconfigdir,SHELFICEtopoFile),ieee,prec);
clear data

data = hFacC;
writeDataset(data,fullfile(inputconfigdir,'hFacC.bin'),ieee,prec);
clear data;



















