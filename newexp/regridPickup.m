%%%
%%% regridPickup.m
%%%
%%% Regrids a pickup file HORIZONTALLY onto a grid of a different size.
%%% (I.e. this script assumes that the vertical grids are identical)
%%%

%%% For file I/O
addpath ../newexp_utils/
addpath ../utils/matlab

%%% Configuration
defineGrid;

%%% Load low-resolution grid
expdir_lo = '../experiments/';
expname_lo = 'hires_seq_onetwelfth_notides_RTOPO2';
inputdir = fullfile(expdir_lo,expname_lo,'results');
XC_lo = rdmds(fullfile(inputdir,'XC'));
YC_lo = rdmds(fullfile(inputdir,'YC'));
RC_lo = rdmds(fullfile(inputdir,'RC'));
hFacC_lo = rdmds(fullfile(inputdir,'hFacC'));

% inputdir = '/data3/MITgcm_WS/experiments/hires_seq_onethird/results';
% inputdir = '/data3/MITgcm_WS/experiments/hires_seq_onesixth_modbdrystrat/results';
% inputdir = '/data3/MITgcm_WS/experiments/hires_seq_onethird_RTOPO2/results';
% inputdir = '/data3/MITgcm_WS/experiments/hires_seq_onesixth_RTOPO2/results';
% inputdir = '/data3/MITgcm_WS/experiments/hires_seq_onetwelfth_RTOPO2/results';
% inputdir = '/data3/MITgcm_WS/experiments/hires_seq_onesixth_notides_RTOPO2/results';

% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onesixth_nobdrymods/results';
% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onetwelfth/results';
% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onethird_RTOPO2/results';
% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onetwelfth_RTOPO2/results';
% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onetwentyfourth_notides_RTOPO2/results';
% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onetwentyfourth_RTOPO2/results';
% outputdir = '/data3/MITgcm_WS/experiments/hires_seq_onetwelfth_notides_RTOPO2/results'
% outputdir = '/data3/MITgcm_WS/experiments/hires_nest_onethirtieth_notides_RTOPO2/results';
outputdir = '/data3/MITgcm_WS/experiments/hires_nest_onethirtysecond_notides_RTOPO2/results';

% expiter_in = 1183320;
% expiter_in = 2366640;
% expiter_in = 591660;
expiter_in = 262960; %%% 1 years into 1/12 run, i.e. start of 2008
expiter_out = 1;
Nx_out = Nx;
Ny_out = Ny;
Nx_in = size(XC_lo,1);
Ny_in = size(YC_lo,2);
pickupfiles = {'pickup','pickup_seaice'};
% pickupfiles = {'pickup_seaice'};
nan_zeros = [true,false]; %%% Set true to replace zeros with NaNs before interpolation. 
                          %%% Usually a good idea e.g. to prevent tiny salinities when 
                          %%% interpolating near land points, but may be undesireable 
                          %%% for sea ice fields, which can be zero in non-land points.                          
% nan_zeros = [false];
write_init_files = false;

%%% Formatting
ieee='b';
prec='real*8';

%%% Loop through pickup files for the specified iteration and regrid each
for m = 1:length(pickupfiles)
  
%   %%% Load all data from pickup file
%   pdata_in = rdmds(fullfile(inputdir,pickupfiles{m}),expiter_in);
%   Nx_in = size(pdata_in,1);
%   Ny_in = size(pdata_in,2);
%   Nr = size(pdata_in,3);
%   pdata_out = NaN*ones(Nx_out,Ny_out,Nr);
  
  %%% Grids for linear interpolation (this is rather crude)
  xx_out = xmc;
  yy_out = ymc;
  xx_in = XC_lo(:,1);
  yy_in = YC_lo(1,:);
  [XX_in,YY_in] = meshgrid(xx_in,yy_in);
  [XX_out,YY_out] = meshgrid(xx_out,yy_out);
  
  %%% Open I/O streams
  fid_in = fopen(fullfile(inputdir,[pickupfiles{m},'.',num2str(expiter_in,'%.10d'),'.data']),'r','b');
  fid_out = fopen(fullfile(outputdir,[pickupfiles{m},'.',num2str(expiter_out,'%.10d'),'.data']),'w','b');
  
  %%% Loop through all columns of pickup data (i.e. through all vertical
  %%% levels of all variables) and interpolate
  k = 0;
  pdata_in = fread(fid_in,[Nx_in Ny_in],'real*8');
  while (size(pdata_in) == [Nx_in Ny_in])
    k = k+1

    %%% Mask land points with NaNs
    if (nan_zeros(m))
      pdata_in(pdata_in==0) = NaN;
    end

    %%% Interpolate
    pdata_out = interp2(XX_in,YY_in,pdata_in',XX_out,YY_out,'linear')'; %%% Interpolation
    pdata_out = inpaint_nans(pdata_out,4); %%% Extrapolation

    %%% Write initialization .bin files
    if (write_init_files)
      %%% TODO
    %%% Write pickup.****.data file
    else
      fwrite(fid_out,pdata_out,'real*8');
    end
    
    %%% Read next level of input file+1
    pdata_in = fread(fid_in,[Nx_in Ny_in],'real*8');
    
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
