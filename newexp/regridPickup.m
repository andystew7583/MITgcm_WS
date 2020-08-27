%%%
%%% regridPickup.m
%%%
%%% Regrids a pickup file HORIZONTALLY onto a grid of a different size.
%%%

%%% For file I/O
addpath ../newexp_utils/
addpath ../utils/matlab

%%% Configuration
% inputdir = '/data3/MITgcm_WS/experiments/hires_seq_onethird/results';
% inputdir = '/data3/MITgcm_WS/experiments/hires_seq_onesixth_modbdrystrat/results';
inputdir = '/data3/MITgcm_WS/experiments/hires_seq_onethird_RTOPO2/results';
% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onesixth_nobdrymods/results';
% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onetwelfth/results';
% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onethird_RTOPO2/results';
outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onesixth_RTOPO2/results';
% expiter_in = 1183320;
% expiter_in = 2366640;
expiter_in = 591660;
expiter_out = 1;
% Nx_out = 1224;
% Ny_out = 912;
% Nx_out = 304;
% Ny_out = 232;
Nx_out = 612;
Ny_out = 432;
pickupfiles = {'pickup','pickup_seaice'};
nan_zeros = [true,false]; %%% Set true to replace zeros with NaNs before interpolation. 
                          %%% Usually a good idea e.g. to prevent tiny salinities when 
                          %%% interpolating near land points, but may be undesireable 
                          %%% for sea ice fields, which can be zero in non-land points.                          
write_init_files = false;

%%% Formatting
ieee='b';
prec='real*8';

%%% Loop through pickup files for the specified iteration and regrid each
for m = 1:length(pickupfiles)
  
  %%% Load all data from pickup file
  pdata_in = rdmds(fullfile(inputdir,pickupfiles{m}),expiter_in);
  Nx_in = size(pdata_in,1);
  Ny_in = size(pdata_in,2);
  Nr = size(pdata_in,3);
  pdata_out = NaN*ones(Nx_out,Ny_out,Nr);
  
  %%% Grids for linear interpolation (this is rather crude)
  dx_in = 1.0/Nx_in;
  dy_in = 1.0/Ny_in;
  xx_in = 0.5*dx_in:dx_in:1-0.5*dx_in;
  yy_in = 0.5*dy_in:dy_in:1-0.5*dy_in;
  dx_out = 1.0/Nx_out;
  dy_out = 1.0/Ny_out;
  xx_out = 0.5*dx_out:dx_out:1-0.5*dx_out;
  yy_out = 0.5*dy_out:dy_out:1-0.5*dy_out;
  [XX_in,YY_in] = meshgrid(xx_in,yy_in);
  [XX_out,YY_out] = meshgrid(xx_out,yy_out);
  
  %%% Mask land points with NaNs
  if (nan_zeros(m))
    pdata_in(pdata_in==0) = NaN;
  end
  
  %%% Loop through all columns of pickup data (i.e. through all vertical
  %%% levels of all variables) and interpolate
  for k = 1:size(pdata_in,3)
    k
    pdata_out(:,:,k) = interp2(XX_in,YY_in,pdata_in(:,:,k)',XX_out,YY_out,'linear')'; %%% Interpolation
    pdata_out(:,:,k) = inpaint_nans(pdata_out(:,:,k),4); %%% Extrapolation
  end
  
  %%% Write initialization .bin files
  if (write_init_files)
    %%% TODO
  %%% Write pickup.****.data file
  else
    writeDataset(pdata_out,fullfile(outputdir,[pickupfiles{m},'.',num2str(expiter_out,'%.10d'),'.data']),ieee,prec);
  end
  
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
