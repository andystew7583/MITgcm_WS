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
inputdir = fullfile(expdir_lo,expname_lo,'results');
XC_lo = rdmds(fullfile(inputdir,'XC'));
YC_lo = rdmds(fullfile(inputdir,'YC'));
RC_lo = rdmds(fullfile(inputdir,'RC'));
hFacC_lo = rdmds(fullfile(inputdir,'hFacC'));
hFacW_lo = rdmds(fullfile(inputdir,'hFacW'));
hFacS_lo = rdmds(fullfile(inputdir,'hFacS'));

outputdir = '/data/data3/MITgcm_WS/experiments/hires_nest_oneseventysecond_notides_RTOPO2/results';

expiter_in = 1578034; %%% 3 years into 1/24 run, i.e. start of 2011
expiter_out = 1;
Nx_out = Nx;
Ny_out = Ny;
Nr_out = Nz;
Nx_in = size(XC_lo,1);
Ny_in = size(YC_lo,2);
Nr_in = length(RC_lo);
pickupfiles = {'pickup','pickup_seaice'};
pickupIs3D = {[1 1 1 1 1 1 0 0 0],[0 0 0 0 0 0 0]}; % Flags whether fields are 2D or 3D
pickupFace = {{'W','S','C','C','W','S','C','C','C'},{'C','C','C','C','C','W','S'}};                          
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
  

  %%% Open I/O streams
  fid_in = fopen(fullfile(inputdir,[pickupfiles{m},'.',num2str(expiter_in,'%.10d'),'.data']),'r','b');
  fid_out = fopen(fullfile(outputdir,[pickupfiles{m},'.',num2str(expiter_out,'%.10d'),'.data']),'w','b');
  
  %%% Loop over fields in the pickup file
  is_3D = pickupIs3D{m};
  for n = 1:length(is_3D)

    %%% Select hFac matrix based on placement of current field on model
    %%% C-grid
    switch (pickupFace{m}{n})
      case 'C'
        hFac = hFacC;
      case 'W'
        hFac = hFacW;
      case 'S'
        hFac = hFacS;
    end

    %%% 3D fields
    if (is_3D(n))

      %%% Load 3D field
      pdata_in = zeros(Nx_in,Ny_in,Nr_in);
      for k = 1:Nr_in
        pdata_in(:,:,k) = fread(fid_in,[Nx_in Ny_in],'real*8');
      end

      %%% 3D meshgrids for interpolation
      [XX_in,YY_in,ZZ_in] = meshgrid(xx_in,yy_in,zz_in);
      [XX_out,YY_out,ZZ_out] = meshgrid(xx_out,yy_out,zz_in);

      %%% Remove land points
      pdata_in(hFac==0) = NaN;

      %%% Do 3D interpolation
      pdata_out = interpPickup3 (pdata_in,XX_in,YY_in,ZZ_in,XX_out,YY_out,ZZ_out);

    %%% 2D fields
    else

      %%% Read 2D field as the next slice in the pickup file
      pdata_in = fread(fid_in,[Nx_in Ny_in],'real*8');

      %%% Grids for linear interpolation (this is rather crude)      
      [XX_in,YY_in] = meshgrid(xx_in,yy_in);
      [XX_out,YY_out] = meshgrid(xx_out,yy_out);

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
      pdata_out = interpPickup2 (pdata_in,XX_in,YY_in,XX_out,YY_out);

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
function data_interp = interpPickup2 (data,XD,YD,XI,YI)
   
  %%% To store interpolated data
  data_interp = zeros(size(XI,2),size(XI,1));
  
  %%% Do 2D interpolation
  data_interp(:,:) = interp2(XD,YD,data(:,:)',XI,YI,'linear')';

  %%% Crude extrapolation
  data_interp(:,:) = inpaint_nans(data_interp(:,:),4); 

end


%%%
%%% interpPickup3
%%%
%%% Convenience function to interpolate three-dimensional pickup fields.
%%%
function data_interp = interpPickup3 (data,XD,YD,ZD,XI,YI,ZI)

  %%% Prepare data matrix for interpolation - needs to match meshgrid
  %%% orientation
  data = transpose3D(data);

  %%% To store interpolated data
  data_interp = interp3(XD,YD,ZD,data,XI,YI,ZI,'linear')';
  
  %%% Crude extrapolation
  for k = 1:size(data_interp,3)
    data_interp(:,:,k) = inpaint_nans(data_interp(:,:,k),4); 
  end

end