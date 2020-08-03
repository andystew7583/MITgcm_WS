%%%v
%%% smooth2D.m
%%%
%%% Smooths 2D data using a Gaussian kernel function.
%%%
%%% XC,YC       meshgrid of physical locations
%%% DXG,DYG     Dimensions of grid cells
%%% FF          2D data to smooth, must be the same size as XX and YY
%%% sigma       Physical width of smoothing filter
%%%
function FF_smooth = smooth2D(XC,YC,DXG,DYG,FF,sigma)

  %%% NaN grid spacings where FF is NaN so that the weighting function is
  %%% calculated properly
  DXG(isnan(FF)) = NaN;
  DYG(isnan(FF)) = NaN;
  
  %%% Grid sizess
  Nx = size(XC,1);
  Ny = size(XC,2);
  
  %%% Smooth the streamfunction using a Gaussian kernel
  FF_smooth = 0*FF;
  for i=1:Nx    
    i
    for j=1:Ny
      
      if (isnan(FF(i,j)))
        continue;
      end
      
      XX_ker = XC - XC(i,j); %%% Coordinates for the kernel function
      YY_ker = YC - YC(i,j);
      
      kerf = (1./(pi*sigma^2)).*exp(-(XX_ker.^2+YY_ker.^2)./(sigma^2)); %%% Defines the filter
      
      tmp = kerf.*DXG.*DYG;
      FF_smooth(i,j) = 1 ./ nansum(tmp(:)); %%% Denominator, ensures proper weighting of the average
      
      tmp = FF.*tmp;
      FF_smooth(i,j) = FF_smooth(i,j) .* nansum(tmp(:)); %%% Numerator, convolution of FF with the exponential filter    
      
    end
  end

end

