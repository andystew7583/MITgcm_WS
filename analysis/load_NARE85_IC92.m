%%%
%%% load_NARE85_IC92.m
%%%
%%% Loads cast data from the NARE85 and IC92 cruises.
%%%

%%% Pull data CTD data files
NARE85_casts = readCTDdata('../data/CTD/NARE85_phys_oce.txt');
IC92_casts = readCTDdata('../data/CTD/IC92_phys_oce.txt');

%%% Extract specific zonal sections that we need
Nz = 200;
NARE85_plotidx = 4:15;
[NARE85_HH,NARE85_ZZ,NARE85_pt,NARE85_ss,NARE85_pd] = extractSection(NARE85_casts,NARE85_plotidx,Nz,false);
IC92_plotidx = [32    34    36    38    40    42    44    48    50    52    54    56];
[IC92_HH,IC92_ZZ,IC92_pt,IC92_ss,IC92_pd] = extractSection(IC92_casts,IC92_plotidx,Nz,true);






%%%
%%% Convenience function to load cast data from tab-delimited text file.
%%%
function casts = readCTDdata (fname)

  %%% Open data file
  fid = fopen(fname,'r');

  %%% To store cast data
  casts = cell(0);
  ncasts = 0;

  %%% Loop over lines of data and assign to casts
  dline = fgetl(fid);
  firstLine = true;
  while (dline~=-1)

    %%% Split line by whitespace
    dline_split = strsplit(dline);

    %%% Check whether we have moved to the next cast yet
    castname = dline_split{1};  
    if (firstLine || (~strcmp(castname,casts{ncasts}.name)))
      newcast.name = castname;
      newcast.lat = [];
      newcast.lon = [];
      newcast.temp = [];
      newcast.salt = [];
      newcast.maxdepth = [];
      newcast.depth = [];
      newcast.pres = [];
      newcast.pottemp = [];
      newcast.sigma = [];
      casts = {casts{:},newcast}; %%% Create a new cast if need be
      ncasts = ncasts + 1;
    end
    
    if (length(dline_split)<10)
      dline
      dline_split
    end

    %%% Load data values from this line
    casts{ncasts}.lat = [casts{ncasts}.lat str2num(dline_split{3})];
    casts{ncasts}.lon = [casts{ncasts}.lon str2num(dline_split{4})];
    casts{ncasts}.temp = [casts{ncasts}.temp str2num(dline_split{8})];
    casts{ncasts}.salt = [casts{ncasts}.salt str2num(dline_split{10})];  
    casts{ncasts}.maxdepth = [casts{ncasts}.maxdepth str2num(dline_split{5})];  
    casts{ncasts}.depth = [casts{ncasts}.depth str2num(dline_split{6})];
    casts{ncasts}.pres = [casts{ncasts}.pres str2num(dline_split{7})];  
    casts{ncasts}.pottemp = [casts{ncasts}.pottemp str2num(dline_split{11})];
    casts{ncasts}.sigma = [casts{ncasts}.sigma str2num(dline_split{12})];

    %%% Read next line
    dline = fgetl(fid);
    firstLine = false;

  end

  %%% Close the file
  fclose(fid);

end







%%%
%%% Convenience function to extract a section from cast data.
%%%
function [HH,ZZ,pt,ss,pd] = extractSection (casts,plotidx,Nz,do_sort)

  smoothwidth = 60;
  smoothmethod = 'moving';

  %%% Extract data from casts, smoothing each cast vertically and
  %%% interpolating linearly onto a grid with a uniform number of vertical
  %%% points
  Nh = length(plotidx);
  HH = zeros(Nh,Nz);
  ZZ = zeros(Nh,Nz);
  pt = zeros(Nh,Nz);
  ss = zeros(Nh,Nz);
  pd = zeros(Nh,Nz);
  for n = 1:length(plotidx)
    ncast = plotidx(n);
    zmax = casts{ncast}.depth(end);
%     zmax = 3500;
    dz = zmax/Nz;  
    ZZ(n,:) = dz:dz:zmax;
    HH(n,:) = interp1(casts{ncast}.depth,casts{ncast}.lon,ZZ(n,:),'linear');
    pt(n,:) = interp1(casts{ncast}.depth,smooth(casts{ncast}.pottemp,smoothwidth,smoothmethod),ZZ(n,:),'linear');    
    ss(n,:) = interp1(casts{ncast}.depth,smooth(casts{ncast}.salt,smoothwidth,smoothmethod),ZZ(n,:),'linear');   
    pd(n,:) = interp1(casts{ncast}.depth,smooth(casts{ncast}.sigma,smoothwidth,smoothmethod),ZZ(n,:),'linear');
  end  
  
  %%% Sort in longitudinal order
  if (do_sort)
    [unused,sortidx] = sort(HH(:,1));
    HH = HH(sortidx,:);
    ZZ = ZZ(sortidx,:);
    pt = pt(sortidx,:);
    ss = ss(sortidx,:);
    pd = pd(sortidx,:);
  end
  
end