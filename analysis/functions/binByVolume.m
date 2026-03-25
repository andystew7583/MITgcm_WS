%%%
%%% binByVolume.m
%%%
%%% Computes total volume in a series of arbitrary property-property bins.
%%%
%%% S3D,T3D - 3D matrices (of same size) of properties to bin
%%% W3D - 3D matrix of weights to multiply the volumes. If empty then ones
%%%       will be used.
%%% Smin,Smax,dS,Tmin,Tmax,dT - Define grid of bins in S/T space
%%% RAC,DRF,hFacC - MITgcm grids that define cell volumes
%%%
function TSvol = binByVolume(S3D,T3D,W3D, ...
                  Smin,Smax,dS,Tmin,Tmax,dT, ...
                  RAC,DRF,hFacC)

  %%% Physical grid dimensions
  Nx = size(S3D,1);
  Ny = size(S3D,2);
  Nr = size(S3D,3);
  
  %%% T/S grids
  SS = Smin:dS:Smax;
  TT = Tmin:dT:Tmax;
  NS = length(SS);
  NT = length(TT);
  
  %%% Weights
  if (isempty(W3D))
    W3D = ones(Nx,Ny,Nr);
  end

  %%% Do volume binning
  TSvol = zeros(NS,NT);
  for i = 1:Nx
    i
    for j = 1:Ny
      for k = 1:Nr    

        %%% Extract properties in this cell
        Sval = S3D(i,j,k);
        Tval = T3D(i,j,k);

        %%% Check this is a wet cell and that properties lie within ranges of
        %%% the TS diagram
        if ( isnan(Sval) || isnan(Tval) ...
            || (Sval<SS(1)-0.5*dS) ...
            || (Sval>SS(end)+0.5*dS) ...
            || (Tval<TT(1)-0.5*dT) ...
            || (Tval>TT(end)+0.5*dT) )
          continue;
        end

        %%% Add this volume to the appropriate T/S gridcell
        ns = find((Sval>SS-0.5*dS) & (Sval<SS+0.5*dS),1);
        nt = find((Tval>TT-0.5*dT) & (Tval<TT+0.5*dT),1);
        TSvol(ns,nt) = TSvol(ns,nt) + W3D(i,j,k)*RAC(i,j)*DRF(k)*hFacC(i,j,k);   

      end
    end    
  end  

end

