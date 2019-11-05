%%% Defining topography 
%%% Also setting intial conditions for ice shelf 



%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%% SET UP %%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%    

% Getting our grid

run defineGrid.m

%%% Gaussian kernel smoothing width (in grid points)
gauss_width = 0; %%% Original width 1.5
gauss_width_shelfice = 0;
min_ice_thickness = 10; %%% Minimum shelf ice thickness in m
% gauss_width = 4;

%%% Set true to smooth water column thickness
smooth_depth_only = false;

%%% Set true to remove coastal embayment in SW Weddell Sea, which is prone
%%% to spurious build-up of ice
remove_embayment = true;
embayment_lat = -74; %%% Max latitude to search for the SW embayment. Original was -74.2.




%%%% Setting Directories

datadir = fullfile(gendir,'/MITgcm_WS/data/RTOPO/');

addpath ../newexp_utils




%Reading in grid from RTOPO dataset

RTOPO_latitude = ncread(fullfile(datadir,'RTOPO.nc'),'lat');
RTOPO_longitude = ncread(fullfile(datadir,'RTOPO.nc'),'lon');
RTOPO_bathymetry = ncread(fullfile(datadir,'RTOPO.nc'),'bathy');
RTOPO_ice = ncread(fullfile(datadir,'RTOPO.nc'),'draft');
RTOPO_elev = ncread(fullfile(datadir,'RTOPO.nc'),'height');
%%% Remove superfluous additional longitude value
RTOPO_longitude = RTOPO_longitude(1:end-1);
RTOPO_bathymetry = RTOPO_bathymetry(1:end-1,:);
RTOPO_ice = RTOPO_ice(1:end-1,:);
RTOPO_elev = RTOPO_elev(1:end-1,:);

%%% Do the twist
RTOPO_bathymetry = Tweddell(RTOPO_longitude,RTOPO_latitude,RTOPO_bathymetry);
RTOPO_ice = Tweddell(RTOPO_longitude,RTOPO_latitude,RTOPO_ice);
RTOPO_elev = Tweddell(RTOPO_longitude,RTOPO_latitude,RTOPO_elev);


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
      MDepth(j,i) = 0;
      Mice_draft(j,i) = 0;
%       Mice_draft(j,i) = MDepth(j,i);
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAUSSIAN KERNAL REGRESSION %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


XMC_smooth = repmat((1:Nx),[Ny 1]);
YMC_smooth = repmat((1:Ny)',[1 Nx]);


if (smooth_depth_only)
  
  %%% Apply smoothing filter
  if (gauss_width > 0)
    Depth_smooth = smooth2D(XMC_smooth',YMC_smooth',MDepth',gauss_width,gauss_width)';    
  else
    Depth_smooth = MDepth;
  end
    
  %%% No smoothing under ice shelves  
  Depth_smooth(Mice_draft<0 | MDepth>=0) = MDepth(Mice_draft<0 | MDepth>=0);  
  Depth_smooth = Depth_smooth';  
  Smooth_Mice = Mice_draft';
  
else

  %%% VERSION 1 - SMOOTH ICE AND BATHYMETRY SEPARATELY %%%
  
%   %%% Smoothing bathymetry
%   if (gauss_width > 0)
%     Depth_smooth = smooth2D(XMC_smooth',YMC_smooth',MDepth',gauss_width,gauss_width)';
%     Depth_smooth = Depth_smooth';
%   else
%     Depth_smooth = MDepth';
%   end
%   
%   %%%%%%% Smoothing Ice Draft
%   if (gauss_width_shelfice > 0) 
%     Smooth_Mice = smooth2D(XMC_smooth',YMC_smooth',Mice_draft',gauss_width_shelfice,gauss_width_shelfice)';
%     Smooth_Mice = Smooth_Mice';
%   else
%     Smooth_Mice = Mice_draft';
%   end


%%% VERSION 2 - SMOOTH BATHYMETRY AND WATER COLUMN THICKNESS %%%

%   %%% Smoothing bathymetry
%   if (gauss_width > 0)
%     Depth_smooth = smooth2D(XMC_smooth',YMC_smooth',MDepth',gauss_width,gauss_width)';
%     Depth_smooth = Depth_smooth';
%   else
%     Depth_smooth = MDepth';
%   end
%   
%   %%% Thickness of water column between smoothed sea floor and unsmoothed
%   %%% ice draft
%   thickness = Mice_draft - Depth_smooth';
%   
%   %%%%%%% Smoothing water column thickness
%   if (gauss_width_shelfice > 0) 
%     thickness_smooth = smooth2D(XMC_smooth',YMC_smooth',thickness',gauss_width_shelfice,gauss_width_shelfice)';
%     thickness_smooth = thickness_smooth';
%   else
%     thickness_smooth = thickness_smooth';
%   end
%   
%   %%% Now reconstruct ice shelf draft from smoothed water column thickness
%   Smooth_Mice = Depth_smooth + thickness_smooth;

  
  %%% VERSION 3 - SMOOTH ICE AND WATER COLUMN THICKNESS %%%
  
  %%% First, remove all ice draft above sea level to remove errors where
  %%% ocean meets land with ice, but there is no ice shelf
  Mice_draft_noland = Mice_draft;
  Mice_draft_noland(Mice_draft_noland > 0) = 0;

  %%%%%%% Smoothing Ice Draft
  if (gauss_width_shelfice > 0) 
    Smooth_Mice = smooth2D(XMC_smooth',YMC_smooth',Mice_draft_noland',gauss_width_shelfice,gauss_width_shelfice)';
    Smooth_Mice = Smooth_Mice';
  else
    Smooth_Mice = Mice_draft_noland';
  end
    
  %%% Eliminate points with very thin ice shelves
  Smooth_Mice(Smooth_Mice>=-min_ice_thickness) = 0;
  
  %%% Thickness of water column between smoothed ice shelf and unsmoothed
  %%% sea floor
  thickness = Smooth_Mice' - MDepth;
  
  %%%%%%% Smoothing water column thickness
  if (gauss_width_shelfice > 0) 
    thickness_smooth = smooth2D(XMC_smooth',YMC_smooth',thickness',gauss_width_shelfice,gauss_width_shelfice)';
    thickness_smooth = thickness_smooth';
  else
    thickness_smooth = thickness';
  end  
  thickness_smooth(thickness_smooth<0) = 0; %%% Can't have negative thickness
  
  %%% Now reconstruct bathymetry from smoothed water column thickness
  Depth_smooth = Smooth_Mice - thickness_smooth;

end


%%%%%% Anywhere that the original ice draft was above land, set equal to
%%%%%% zero in smoothed version.
% Smooth_Mice(Mice_draft'>=0) = 0;
  
%%%%%% Remove points where depth or ice draft lie above the surface
Depth_smooth(Depth_smooth>0) = 0;
Smooth_Mice(Smooth_Mice>0) = 0;

%%%%%% Smoothed Ice draft and Smoothed bathymetry match
%%%%%% where the original draft minus the bathymetry equals zero.
Smooth_Mice((Mice_draft-MDepth)'==0) = Depth_smooth((Mice_draft-MDepth)'==0);

%%% Avoid negative water column thicknesses
Smooth_Mice(Smooth_Mice-Depth_smooth < 0) = Depth_smooth(Smooth_Mice-Depth_smooth < 0);






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

data = Depth_smooth;
writeDataset(data,fullfile(inputfolder,bathyFile),ieee,prec);
clear data


data = Smooth_Mice;
writeDataset(data,fullfile(inputfolder,SHELFICEtopoFile),ieee,prec);
clear data

tides =1;
if tides ==1
    addpath ../newexp_utils
    data = hFacC;
    writeDataset(data,fullfile(inputfolder,'hFacC.bin'),ieee,prec);
    clear data;
end


















