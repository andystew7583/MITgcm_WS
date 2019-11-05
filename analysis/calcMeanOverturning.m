%%%
%%% calcMeanOverturning.m
%%%
%%% Calculates the Eulerian-mean overturning circulation.
%%%

%%% To store Eulerian-mean streamfunction
psimean = zeros(Ny+1,Nr+1);

%%% Integrate in z to calculate psim
for j=1:Ny  
  
  %%% Calculate topographic depths at v-gridpoint
  hFacS_col = squeeze(hFacS(1,j,:));  
  hb_v = sum(delR.*squeeze(hFacS(1,j,:))');         
  kmax = length(hFacS_col(hFacS_col>0));    
  
  %%% Do the integration
  psimean(j,1) = 0;
  for k=1:Nr
    psimean(j,k+1) = psimean(j,k) + vv_avg(j,k) * delR(k)*hFacS_col(k);        
  end  
  
  %%% If not all vertical grid cells are used at this latitude then set
  %%% those that lie within topography to NaN.
  if (kmax < Nr)
    psimean(j,kmax+2:Nr+1) = NaN;
  end
  
end
psimean(Ny+1,:) = psimean(1,:);

