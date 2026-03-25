

% loadExp;

nPx = 8; %%% no. of processors in x-direction
nPy = 48; %%% no. of processors in y-direction
sNx = Nx / nPx;
sNy = Ny / nPy;

nWet = zeros(nPx,nPy);
for i=1:nPx
  for j=1:nPy
    nWet(i,j) = sum(sum(sum(ceil(hFacC((i-1)*sNx+1:i*sNx,(j-1)*sNy+1:j*sNy,:)))));
  end
end

nDry = length(find(nWet==0))
nPx*nPy-nDry
(sNx+6)*(sNy+6)
