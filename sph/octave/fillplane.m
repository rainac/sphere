% -*- octave -*-
% Johannes Willkomm 2008
function res = fillplane(delta, pnt_ul, pnt_or)

  dims = find(pnt_ul != pnt_or);

  lx = pnt_or(dims(1)) - pnt_ul(dims(1));
  ly = pnt_or(dims(2)) - pnt_ul(dims(2));
  
  rx = 0:delta:lx;
  ry = 0:delta:ly;
  
  nx = length(rx);
  ny = length(ry);
  
  res = zeros(nx * ny, 3);
  
  for i=1:nx
    for j=1:ny
      res((i-1)*ny + j, :) = pnt_ul;
      res((i-1)*ny + j, dims(1)) = pnt_ul(dims(1)) + rx(i);
      res((i-1)*ny + j, dims(2)) = pnt_ul(dims(2)) + ry(j);
    end
  end
