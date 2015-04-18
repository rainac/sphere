% -*- octave -*-
function r = plotkernel(fname, kname, ndim, H, Z)
  ktitle = sprintf("%dD %s kernel, H = %g", ndim, kname, H);
  if (ndim == 3) 
    ktitle = [ktitle sprintf(", z = %g", Z)];
  endif
  s = load(fname);
  N = size(s, 1);
  n = sqrt(N);
  s = reshape(s, n, n, size(s, 2));
  mesh(s(1,:,1), s(1,:,1), s(:,:,4));
  title(sprintf('%s', ktitle));
  xlabel('X-Coord');
  ylabel('Y-Coord');
  z = zlabel('W(x,y,z)');
  set(z, 'rotation', 90);
  print([fname '.f.eps'], '-deps2c');
  print([fname '.f.png'], '-dpng');
  figure
  mesh(s(1,:,1), s(1,:,1), s(:,:,5));
  title(sprintf('%s, gradient (x-Coord)', ktitle));
  xlabel('x-Coord');
  ylabel('y-Coord');
  z = zlabel('d W(x,y,z) / d x');
  set(z, 'rotation', 90);
  print([fname '.gx.eps'], '-deps2c');
  print([fname '.gx.png'], '-dpng');
  figure
  mesh(s(1,:,1), s(1,:,1), s(:,:,6));
  title(sprintf('%s, gradient (y-Coord)', ktitle));
  xlabel('x-Coord');
  ylabel('y-Coord');
  z = zlabel('d W(x,y,z) / d y');
  set(z, 'rotation', 90);
  print([fname '.gy.eps'], '-deps2c');
  print([fname '.gy.png'], '-dpng');
  figure
  mesh(s(1,:,1), s(1,:,1), s(:,:,7));
  title(sprintf('%s, gradient (z-Coord)', ktitle));
  xlabel('x-Coord');
  ylabel('y-Coord');
  z = zlabel('d W(x,y,z) / d z');
  set(z, 'rotation', 90);
  print([fname '.gz.eps'], '-deps2c');
  print([fname '.gz.png'], '-dpng');
