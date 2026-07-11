%%%
%%% regridPickup.m
%%%
%%% Regrids a pickup file horizontall and vertically onto a grid of a different size.
%%%

%%% For file I/O
addpath ../newexp_utils/
addpath ../utils/matlab

%%% Configuration
defineGrid;

%%% Load low-resolution grid
expdir_lo = '../experiments/';
expname_lo = 'hires_seq_onetwentyfourth_notides_RTOPO2';
% expname_lo = 'hires_seq_onethird_notides_RTOPO2';
inputdir = fullfile(expdir_lo,expname_lo,'results');
XC_lo = rdmds(fullfile(inputdir,'XC'));
YC_lo = rdmds(fullfile(inputdir,'YC'));
RC_lo = rdmds(fullfile(inputdir,'RC'));
hFacC_lo = rdmds(fullfile(inputdir,'hFacC'));
hFacW_lo = rdmds(fullfile(inputdir,'hFacW'));
hFacS_lo = rdmds(fullfile(inputdir,'hFacS'));

outputdir = '/data/data3/MITgcm_WS/experiments/hires_nest_oneseventysecond_notides_RTOPO2/results';

expiter_in = 1578034; %%% 3 years into 1/24 run, i.e. start of 2011
% expiter_in = 1512020; %%% Arbitrarily-chosen test iteration in 1/3 run
expiter_out = 1;
Nx_out = Nx;
Ny_out = Ny;
Nr_out = Nz;
Nx_in = size(XC_lo,1);
Ny_in = size(YC_lo,2);
Nr_in = length(RC_lo);
pickupfiles = {'pickup','pickup_seaice'};
pickupIs3D = {[1 1 1 1 1 1 0 0 0],... % Flags whether fields are 2D or 3D
  [0 0 0 0 0 0 0 0 0 0 0 0 0]}; % N.B. sea ice pickkup has extra levels for different ice categories
pickupFace = {{'W','S','C','C','W','S','C','C','C'},...
{'C','C','C','C','C','C','C','C','C','C','C','W','S'}};                          
pickupIsNonZeroInCavities = [true,false]; % Flags whether 2D fields can be ignored in ice shelf cavities
write_init_files = false;

%%% Grids for interpolation, assuming regular lat/lon grids
xx_out = xmc;
yy_out = ymc;
zz_out = zz;
xx_in = XC_lo(:,1);
yy_in = YC_lo(1,:);
zz_in = squeeze(RC_lo);

%%% Formatting
ieee='b';
prec='real*8';

%%% Loop through pickup files for the specified iteration and regrid each
for m = 1:length(pickupfiles)
  
  disp(['m = ',num2str(m)])

  %%% Open I/O streams
  fid_in = fopen(fullfile(inputdir,[pickupfiles{m},'.',num2str(expiter_in,'%.10d'),'.data']),'r','b');
  fid_out = fopen(fullfile(outputdir,[pickupfiles{m},'.',num2str(expiter_out,'%.10d'),'.data']),'w','b');
  
  %%% Loop over fields in the pickup file
  is_3D = pickupIs3D{m};
  for n = 1:length(is_3D)

    disp(['n = ',num2str(n)])

    %%% Select hFac matrix based on placement of current field on model
    %%% C-grid
    switch (pickupFace{m}{n})
      case 'C'
        hFac = hFacC_lo;
      case 'W'
        hFac = hFacW_lo;
      case 'S'
        hFac = hFacS_lo;
    end

    %%% 3D fields
    if (is_3D(n))

      %%% Load 3D field
      pdata_in = zeros(Nx_in,Ny_in,Nr_in);
      for k = 1:Nr_in
        pdata_in(:,:,k) = fread(fid_in,[Nx_in Ny_in],'real*8');
      end

      %%% Remove land points
      pdata_in(hFac==0) = NaN;

      %%% Do 3D interpolation
      pdata_out = interpPickup3_lowmem (pdata_in,xx_in,yy_in,zz_in,xx_out,yy_out,zz_out);      

    %%% 2D fields
    else

      %%% Read 2D field as the next slice in the pickup file
      pdata_in = fread(fid_in,[Nx_in Ny_in],'real*8');   

      %%% Mask land points with NaNs. For some fields (e.g. SSH), ice shelf 
      %%% cavities don't count as land points, whereas they do for others
      %%% (e.g. sea ice area)
      if (pickupIsNonZeroInCavities(m))
        msk = sum(hFac,3) == 0;
      else
        msk = hFac(:,:,1) == 0;
      end
      pdata_in(msk) = NaN;
      
      %%% Do 2D interpolation
      pdata_out = interpPickup2 (pdata_in,xx_in,yy_in,xx_out,yy_out,true);

    end

    %%% Write pickup.****.data file
    fwrite(fid_out,pdata_out,'real*8');   
    
  end
  
  %%% Close I/O
  fclose(fid_in);
  fclose(fid_out);
  
  %%% Create .meta file
  ifid = fopen(fullfile(inputdir,[pickupfiles{m},'.',num2str(expiter_in,'%.10d'),'.meta']),'r');
  ofid = fopen(fullfile(outputdir,[pickupfiles{m},'.',num2str(expiter_out,'%.10d'),'.meta']),'w');
  iline = fgetl(ifid);
  while (ischar(iline))
    if (startsWith(iline,' dimList')) %%% Update dimension list
      fprintf(ofid,'%s\n',iline);
      fprintf(ofid,'%s\n',['   ',num2str(Nx_out),',    1,  ',num2str(Nx_out),',']);
      fprintf(ofid,'%s\n',['   ',num2str(Ny_out),',    1,  ',num2str(Ny_out)]);
      fgetl(ifid); %%% Discard next two lines because we are replacing them
      fgetl(ifid);
    elseif (startsWith(iline,' timeStepNumber'))
      fprintf(ofid,'%s\n',[' timeStepNumber = [     ',num2str(expiter_out),' ];']);
    else
      fprintf(ofid,'%s\n',iline);
    end
    iline = fgetl(ifid);
  end
  fclose(ifid);
  fclose(ofid);
  
end




%%%
%%% interpPickup2
%%%
%%% Convenience function to interpolate two-dimensional pickup fields.
%%%
function data_interp = interpPickup2 (data,xx_in,yy_in,xx_out,yy_out,extrap)

  %%% Grids for linear interpolation (this is rather crude)      
  [XD,YD] = meshgrid(xx_in,yy_in);
  [XI,YI] = meshgrid(xx_out,yy_out);

  %%% To store interpolated data
  data_interp = zeros(size(XI,2),size(XI,1));
  
  %%% Do 2D interpolation
  data_interp(:,:) = interp2(XD,YD,data(:,:)',XI,YI,'linear')';

  %%% Crude extrapolation
  if (extrap)
    data_interp(:,:) = inpaint_nans(data_interp(:,:),4); 
  end

end


%%%
%%% interpPickup3
%%%
%%% Convenience function to interpolate three-dimensional pickup fields.
%%%
function data_interp = interpPickup3 (data,xx_in,yy_in,zz_in,xx_out,yy_out,zz_out)

  %%% 3D meshgrids for interpolation
  [XD,YD,ZD] = meshgrid(xx_in,yy_in,zz_in);
  [XI,YI,ZI] = meshgrid(xx_out,yy_out,zz_out);

  %%% Prepare data matrix for interpolation - needs to match meshgrid
  %%% orientation
  data = transpose3D(data);

  %%% To store interpolated data
  data_interp = interp3(XD,YD,ZD,data,XI,YI,ZI,'linear');
  
  %%% Crude extrapolation
  for k = 1:size(data_interp,3)
    data_interp(:,:,k) = inpaint_nans(data_interp(:,:,k),4); 
  end

  %%% Revert back to x/y/z matrix format
  data_interp = transpose3D(data_interp);

end

%%%
%%% interpPickup3_lowmem
%%%
%%% Convenience function to interpolate three-dimensional pickup fields.
%%%
function data_interp = interpPickup3_lowmem (data,xx_in,yy_in,zz_in,xx_out,yy_out,zz_out)

  %%% Extract grid dimensions
  Nx_in = length(xx_in);
  Ny_in = length(yy_in);
  Nr_in = length(zz_in);
  Nx_out = length(xx_out);
  Ny_out = length(yy_out);
  Nr_out = length(zz_out);

  %%% First step: interpolate horizontally
  data_2Dinterp = zeros(Nx_out,Ny_out,Nr_in);
  for k = 1:Nr_in
    disp(['interpPickup3_lowmem: k = ',num2str(k)])
    data_2Dinterp(:,:,k) = interpPickup2(data(:,:,k),xx_in,yy_in,xx_out,yy_out,false);
  end

   %%% Second step: interpolate vertically
  data_interp = zeros(Nx_out,Ny_out,Nr_out);
  for k = 1:Nr_out
    if (zz_out(k) > zz_in(1))
      data_interp(:,:,k) = data_2Dinterp(:,:,1);
    elseif (zz_out(k) < zz_in(Nr_in))
      data_interp(:,:,k) = data_2Dinterp(:,:,Nr_in);
    else
      kprev = find(zz_in>zz_out(k),1,'last');
      knext = kprev + 1;
      wtnext = ones(Nx_out,Ny_out)*(zz_in(kprev)-zz_out(k)) / (zz_in(kprev) - zz_in(knext));
      wtprev = ones(Nx_out,Ny_out)*(zz_out(k)-zz_in(knext)) / (zz_in(kprev) - zz_in(knext));
      wtnext(isnan(data_2Dinterp(:,:,knext))) = 0;
      wtprev(isnan(data_2Dinterp(:,:,knext))) = 1;
      wtnext(isnan(data_2Dinterp(:,:,kprev))) = 1;
      wtprev(isnan(data_2Dinterp(:,:,kprev))) = 0;
      data_interp(:,:,k) = wtnext.*data_2Dinterp(:,:,knext) + wtprev.*data_2Dinterp(:,:,kprev);
    end
  end

  %%% Nearest neighbor extrapolation vertically
  for i=1:Nx_out
    for j=1:Ny_out
      kmin = find(~isnan(squeeze(data_interp(i,j,:))),1,'first');
      if (isempty(kmin))
        continue;
      end
      kmax = find(~isnan(squeeze(data_interp(i,j,:))),1,'last');
      if (kmin > 1)
        data_interp(i,j,1:kmin-1) = repmat(data_interp(i,j,kmin),[1 1 kmin-1]);
      end
      if (kmax < Nr_out)
        data_interp(i,j,kmax+1:Nr_out) = repmat(data_interp(i,j,kmax),[1 1 Nr_out-kmax]);
      end
    end
  end 
  
  %%% Crude extrapolation horizontally
  for k = 1:size(data_interp,3)
    data_interp(:,:,k) = inpaint_nans(data_interp(:,:,k),4); 
  end

end
