%%%
%%% smooth2D.m
%%%
%%% Smooths 2D data using a Gaussian kernel function.
%%%
%%% XX,YY       meshgrid of physical locations
%%% FF          2D data to smooth, must be the same size as XX and YY
%%% wx,wy       Physical widths of the smoothing filter
%%%
function FF_smooth = smooth2D(XX,YY,FF,wx,wy)

  Nx = size(XX,1);
  Ny = size(YY,2);
  
  FF_smooth = FF;
  
  for j=1:Ny     
    res = ksr(XX(:,j),FF_smooth(:,j),wx,Nx);
    FF_smooth(:,j) = res.f;
  end
  
  for i=1:Nx
    res = ksr(YY(i,:),FF_smooth(i,:),wy,Ny);
    FF_smooth(i,:) = res.f;
  end
  
end

