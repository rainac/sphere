% function res = genplane(delta, sizes [, axes])
%  generate 3D plane of size sizes centered at (0, 0, 0)
%    delta    -  density of points
%    sizes(2) -  length and width of plane
%    axes(2)  -  axes indizes (between 1 and 3), default: [1 2]
%
% Johannes Willkomm 2008
%
function res = genplane(delta, sizes, axes)
  if (nargin < 3) 
    axes = [1 2];
  end

  if (length(axes) ~= 2)
    error('length(axes) must be 2!')
  end
  
  if (axes(1) == axes(2))
    error('entries of axes must be different!')
  end

  lx = sizes(1);
  ly = sizes(2);
  
  nx = round(lx/delta);
  ny = round(ly/delta);

  dx = lx/nx;
  dy = ly/ny;

  nx=nx+1;
  ny=ny+1;

  res = zeros(nx*ny, 3);
  
  xs = (0:dx:lx) - lx/2;
  ys = (0:dy:ly) - ly/2;

  res(:, axes(1)) = repmat(xs', ny, 1);
  res(:, axes(2)) = reshape(repmat(ys, nx, 1), nx*ny, 1);

% -*- octave -*-
