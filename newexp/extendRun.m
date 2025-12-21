%%%
%%% extendRun.m
%%%
%%% Extends an MITgcm simulation by restarting from the last checkpoint for a selected period of time.
%%%
%%% expdir - full path to folder containing the experiment
%%% expname - name of the experiment
%%% exttime - new end time for the simulation, in seconds
%%% use_shelfice_diags - set true for SHELFICE diagnostics
%%% use_flux_diags - set true for SHELFICE diagnostics
%%% use_eddy_diags - set true for SHELFICE diagnostics
%%% use_layers_diags - set true for SHELFICE diagnostics
%%%
function extendRun (expdir,expname,newEndTime,use_shelfice_diags,use_flux_diags,use_eddy_diags,use_layers_diags)

  %%% Define data format
  defineGrid;

  %%% List terminator character for parameter files - may be '/' or '&'
  %%% depending on operating system
  listterm = '&';

  %%% Directories containing simulation files  
  inputdir = fullfile(expdir,expname,'input');  
  resultsdir = fullfile(expdir,expname,'results');
  codedir = fullfile(expdir,expname,'code');  
  datafile_path = fullfile(inputdir,'data');
  paramsfile = 'params.m';
  tmpparamsfile = 'tmp_params.m';
  diagparamsfile = 'diag_params.m';
  paramsfile_path = fullfile(inputdir,paramsfile);
  tmpparamsfile_path = fullfile(inputdir,tmpparamsfile);
  diagparamsfile_path = fullfile(inputdir,diagparamsfile);
  
  %%% Find the last pickup file in the directory
  resultsFiles = dir(resultsdir);
  newStartIter = -1;
  for m = 1:length(resultsFiles)
    pname = resultsFiles(m).name;
    if (startsWith(pname,'pickup.'))
      startIter = str2num(pname(8:17));
      if (startIter > newStartIter)
        newStartIter = startIter;
      end
    end
  end
  if (newStartIter == -1)
    error(['No numbered pickup files found in ',resultsdir]);
  end

  %%% Set nIter0 to timestepnumber in the simulation's input data file
  fid = fopen(datafile_path,'r');
  if (fid == -1)
    error(['Could not open ',datafile_path]);
  end
  datastr = '';
  tline = fgetl(fid);
  while (ischar(tline))
    if (~isempty(strfind(tline,'nIter0')))
      datastr = [datastr,' nIter0=',num2str(newStartIter),',\n'];
    else
      if (~isempty(strfind(tline,'endTime')))
        datastr = [datastr,' endTime=',num2str(newEndTime),',\n'];
      else
        datastr = [datastr,tline,'\n'];
      end
    end      
    tline = fgetl(fid);
  end
  fclose(fid);
  fid = fopen(datafile_path,'w');
  fprintf(fid,datastr);
  fclose(fid);
  
  %%% Update the simulation's params.m file to reflect the new simulation
  %%% end time and delete diagnostics parameters
  fid = fopen(paramsfile_path,'r');  
  if (fid == -1)
    error(['Could not open ',paramsfile_path]);
  end
  tfid = fopen(tmpparamsfile_path,'w');
  if (tfid == -1)
    error(['Could not open ',tmpparamsfile_path]);
  end
  tline = fgetl(fid);
  while (ischar(tline))
    if (~isempty(strfind(tline,'endTime')))
      paramsstr = ['endTime=',num2str(newEndTime),';'];
    else
      paramsstr = tline;
    end
    if (~startsWith(tline,'diag_'))
      fprintf(tfid,'%s\n',paramsstr);
    end
    tline = fgetl(fid);
  end
  fclose(fid);





  %%% Reset diagnostic parameters %%%
  diag_freq_avg = t1month;
  diag_freq_inst = t1day;

  %%% Configure diagnostic parameter choices
  [DIAG_PARM,DIAG_MATLAB_PARM,ndiags] = setDiagParams(use_shelfice_diags,use_flux_diags,use_eddy_diags,use_layers_diags,diag_freq_avg,diag_freq_inst);

  %%% Replace the data.diagnostics file
  write_data_diagnostics(inputdir,DIAG_PARM,listterm,realfmt);
  
  %%% Replace the DIAGNOSTICS_SIZE.h file
  createDIAGSIZEh(codedir,ndiags,Nr);
  
  %%% Generate a new diag_params.m file that contains only the diagnostics
  %%% parameters
  write_matlab_params(inputdir,diagparamsfile,[DIAG_MATLAB_PARM],realfmt);

  %%% Append new diagnostics parametes to the end of the temp_params.m
  %%% file
  dfid = fopen(diagparamsfile_path,'r');
  if (dfid == -1)
    error(['Could not open ',diagparamsfile_path]);
  end
  fprintf(tfid,'\n');
  dline = fgetl(dfid);
  while (ischar(dline))   
    fprintf(tfid,'%s\n',dline);    
    dline = fgetl(fid);
  end
  fclose(dfid);
  fclose(tfid);  




  %%% Finally, replace old params.m file
  copyfile(tmpparamsfile_path,paramsfile_path);
  delete(tmpparamsfile_path);
  delete(diagparamsfile_path);

    
end


