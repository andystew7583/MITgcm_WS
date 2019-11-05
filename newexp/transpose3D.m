function outmat = transpose3D( inmat )
%TRANSPOSE3D Takes the 2D transpose of a 3D matrix

  outmat = zeros(size(inmat,2),size(inmat,1),size(inmat,3));
  for k=1:size(inmat,3)
    outmat(:,:,k) = inmat(:,:,k)';
  end

end

