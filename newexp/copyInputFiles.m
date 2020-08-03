%%%
%%% copyInputFiles.m
%%% 
%%% Copies input files configured for the current grid configuration into
%%% DEFAULTS/input so that they will be included in experiments generated
%%% via 'newexp'.
%%%

%%% Grid parameters
defineGrid;

%%% Loop through all files in input configuration directory and copy each
%%% to the default input directory
filelist = dir(inputconfigdir);
for n=1:1:length(filelist)
  %%% Ignore hidden files
  if (filelist(n).name(1) == '.')
    continue;
  end    
  copyfile(fullfile(inputconfigdir,filelist(n).name),fullfile(inputfolder,filelist(n).name));
end   

%%% Copy time step file over
copyfile(fullfile(inputconfigdir,'TIME_STEP'),'.');