%
% Johannes Willkomm 2008
function res = dumpascii(v, outf)
 
  nx = size(v, 1);
  ny = size(v, 2);
  nz = size(v, 3);
  
  for i=1:nx
    for j=1:ny
      fprintf(outf, '%g ', v(i,j,k));
    end
    fprintf(outf, '\n');
  end
