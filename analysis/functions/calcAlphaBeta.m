%%%
%%% calcAlphaBeta
%%%
%%% Approximates the Boussinesq thermal expansion and haline contractions 
%%% coefficients, and their derivatives with respect to temperature, 
%%% salinity and depth.
%%%
function [alpha,beta] = calcAlphaBeta(SP,PT,PP)

  dPT = 1e-3;
  dSP = 1e-3;
  alpha = - (densjmd95(SP,PT+0.5*dPT,PP) ...
                - densjmd95(SP,PT-0.5*dPT,PP) ) ...
                ./ densjmd95(SP,PT,PP) / dPT;
  beta = (densjmd95(SP+0.5*dSP,PT,PP) ...
                - densjmd95(SP-0.5*dSP,PT,PP) ) ...
                ./ densjmd95(SP,PT,PP) / dSP;
          
end