%%%
%%% setDiagParams.m
%%%
%%% Configures input parameters for the MITgcm diagnostics package.
%%%
function [DIAG_PARM,DIAG_MATLAB_PARM,ndiags] = setDiagParams(use_shelfice,use_flux_diags,use_eddy_diags,use_layers,diag_freq_avg,diag_freq_inst)


  %%% Get parameter type definitions
  paramTypes;
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  diag_parm01 = parmlist;
  diag_parm02 = parmlist;
  DIAG_PARM = {diag_parm01,diag_parm02};
  diag_matlab_parm01 = parmlist;
  DIAG_MATLAB_PARM = {diag_matlab_parm01}; %%% Matlab parameters need to be different
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%% Stores total number of diagnostic quantities
  ndiags = 0;
  
  
  diag_fields_avg = ...
  { ...
    'UVEL','VVEL','WVEL','THETA','SALT',... %%% 3D fields
    'ETAN', ... 
    'SIheff','SIarea','SIhsnow','SItices','SIhsalt' ... %%ice diagnostics
    'SIuice','SIvice' ...
    'oceFWflx','oceSflux','oceQnet','oceTAUX','oceTAUY','TFLUX','SFLUX',...
    'EXFtaux','EXFtauy','EXFlwnet','EXFswnet','EXFlwdn','EXFswdn','EXFqnet','EXFhs','EXFhl','EXFevap','EXFpreci','EXFatemp' ...
    'SIqnet','SIqsw','SIatmQnt','SItflux','SIaaflux','SIhl','SIqneto','SIqneti','SIempmr','SIatmFW','SIsnPrcp','SIactLHF','SIacSubl'
  };

  if (use_shelfice)
    diag_fields_avg = {diag_fields_avg{:}, ...
    ...
      'SHIfwFlx','SHIhtFlx','SHIForcS','SHIForcT' ...
      % 'SHIUDrag','SHIVDrag'...
    };
  end

  if (use_flux_diags)
    diag_fields_avg = {diag_fields_avg{:}, ...
      'ADVr_TH','ADVr_TH','ADVr_TH', ... %%% Advective heat fluxes
      'ADVr_SLT','ADVr_SLT','ADVr_SLT', ... %%% Advective salt fluxes  
      'DFrI_TH','DFrE_TH', ... %%% Vertical diffusive fluxes of heat
      'DFrI_SLT','DFrE_SLT', ... %%% Vertical diffusive fluxes of salt      
    };
  end
    
  if (use_eddy_diags)
    diag_fields_avg = {diag_fields_avg{:}, ...
      'UVELSLT','VVELSLT','WVELSLT', ... %%% Salt fluxes
      'UVELTH','VVELTH','WVELTH', ... %%% Temperature fluxes  
      'UVELSQ','VVELSQ','WVELSQ', ... %%% For Kinetic Energy
      'UV_VEL_Z','WU_VEL','WV_VEL', ... %%% Momentum fluxes
      'SALTSQ','THETASQ','THSLT' %%% Thermodynamic variances
    };
  end     
  
  if (use_layers)
    diag_fields_avg = {diag_fields_avg{:}, ...
      'LaUH1RHO', 'LaVH1RHO', ... %%% LAYERS fluxes
      'LaHw1RHO', 'LaHs1RHO' ... %%% LAYERS thicknesses
    };
  end  
  
  numdiags_avg = length(diag_fields_avg);  
  diag_phase_avg = 0;    
     
  diag_parm01.addParm('diag_mnc',true,PARM_BOOL);  
  for n=1:numdiags_avg    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields{1,',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames{',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL);  
    
  end
  
%   diag_fields_inst = ...
%   {...      
%     'ETAN', ...
%     'SIheff','SIarea','SIhsnow', ... %%ice diagnostics
% %     'UVEL','VVEL','WVEL'...      
% %     'THETA' ... 
% %     'SALT', ... 
% %     'SIuice','SIvice' ...
% %     'PHIBOT','oceFWflx','oceSflux','oceQnet'...
% %     'SHIfwFlx','SHIhtFlx','SHIUDrag','SHIVDrag' ...
%   }; 

  diag_fields_inst = ...
  { ...
    % 'UVEL','VVEL','WVEL'...      
    % 'THETA' ... 
    % 'SALT','ETAN' ... 
    % 'SIheff','SIarea','SIhsnow','SItices','SIhsalt' ... %%ice diagnostics
    % 'SIuice','SIvice' ...
    % 'oceFWflx','oceSflux','oceQnet','oceTAUX','oceTAUY','TFLUX','SFLUX',...
    % 'EXFtaux','EXFtauy','EXFlwnet','EXFswnet','EXFlwdn','EXFswdn','EXFqnet','EXFhs','EXFhl','EXFevap','EXFpreci','EXFatemp' ...
    % 'SIqnet','SIqsw','SIatmQnt','SItflux','SIaaflux','SIhl','SIqneto','SIqneti','SIempmr','SIatmFW','SIsnPrcp','SIactLHF','SIacSubl', ...    
    % 'SHIfwFlx','SHIhtFlx','SHIUDrag','SHIVDrag' ...    
  };
  
  numdiags_inst = length(diag_fields_inst);    
  diag_phase_inst = 0;
  
  for n=1:numdiags_inst    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);  
    
  end

end

