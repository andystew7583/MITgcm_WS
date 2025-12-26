%%%
%%% setLayersParams.m
%%%
%%% Configures input parameters for the MITgcm LAYERS package.
%%%
function [LAYERS_PARM,Nlayers] = setLayersParams ()

  %%% Get parameter type definitions
  paramTypes;
  
  %%% To store parameter names and values
  layers_parm01 = parmlist;
  LAYERS_PARM = {layers_parm01};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Define parameters for layers package %%%
  
  %%% Specify potential temperature
  layers_name = 'RHO';  
  
  %%% Potential temperature bounds for layers  
  % layers_bounds = [30.5:1:36.5 36.6:0.1:36.8 36.9:0.01:37.4 37.42:0.02:37.6] - 9.38; 
  layers_bounds = [21:1:24 ...
                   24.5:25 ...
                   25.1:.1:27 ...
                   27.05:0.05:27.5 ...
                   27.51:.01:27.75 ...
                   27.755:0.005:27.83 ...
                   27.8325:0.0025:27.9 ...
                   27.91:0.01:28 ...
                   28.02:0.02:28.2 ...
                   28.25:0.05:28.4];
  Nlayers = length(layers_bounds) - 1;

  %%% Reference level for calculation of potential density
  layers_krho = 1;    
  
  %%% If set true, the GM bolus velocity is added to the calculation
  layers_bolus = false;  
   
  %%% Layers
  layers_parm01.addParm('layers_bounds',layers_bounds,PARM_REALS); 
  layers_parm01.addParm('layers_krho',layers_krho,PARM_INT); 
  layers_parm01.addParm('layers_name',layers_name,PARM_STR); 
  layers_parm01.addParm('layers_bolus',layers_bolus,PARM_BOOL); 
  

end

