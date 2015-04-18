function num=numLocalRels(dim)
  n=2^dim;
  num = 0;
  
  for i=0:n-1
    for j=0:i
      if bitand(i,j) == 0
        [i j];
        num = num + 1;
      end
    end
  end
  num