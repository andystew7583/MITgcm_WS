%%%
%%% write_data_cal
%%%
%%% Writes the 'data.cal' input file from 
%%% the cell array of CAL_PARM parmlist object
%%% 'realfmt' specifies the output format for real values,
%%% and 'listterm' specifies the fortran NAMELIST termination character.
%%%
function write_data_cal (dirname,CAL_PARM,listterm,realfmt)

  %%% Open the 'data.kpp' file for writing
  fname = 'data.cal';
  fid = fopen(fullfile(dirname,fname),'w');
  if (fid == -1)
    error(['Could not open ',fname,' file']);
  end
  
  %%% Write out file title text
  titletext = [
    '# ===================\r\n' ...
    '# | CAL parameters  |\r\n' ...
    '# ===================\r\n' ...
    '\r\n' ];
  fprintf(fid,titletext);

  
  %%% Parameter section titles
  titles = {...
    'Calendar parameters'};
  
  %%% For each parameter section 
  for nparam=1:1:length(CAL_PARM)
    
    %%% Write section header
    fprintf(fid,['# ',titles{nparam},'\r\n']);
    fprintf(fid,[' &CAL_NML','\r\n']);   
        
    %%% Write each parameter out to the 'data' file
    nextparmlist = CAL_PARM{nparam};
    for n=1:1:nextparmlist.getLength()      
      writeParam(fid,nextparmlist.getParm(n),realfmt);
    end    
    
    %%% Write section footer
    fprintf(fid,[' ',listterm,'\r\n']);
    fprintf(fid,'\r\n');
    
  end
   
  
  %%% Close the file when we're finished
  fclose(fid);

end