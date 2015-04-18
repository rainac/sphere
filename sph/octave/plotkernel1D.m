% -*- octave -*-
function r = plotkernel1D(fname, kname, H)
  s = load(fname);
  ktitle = sprintf("1D %s kernel, H = %g", kname, H);
  n = size(s, 1);
  plot(s(:,1), s(:,4));
  title(ktitle);
  xlabel('x-Coord');
  ylabel('W(x)');
  print([fname '.f.eps'], '-deps2c');
  print([fname '.f.png'], '-dpng');
  figure
  plot(s(:,1), s(:,5));
  title([ktitle ', gradient']);
  xlabel('x-Coord');
  ylabel('d W(x) / d x');
  print([fname '.gx.eps'], '-deps2c');
  print([fname '.gx.png'], '-dpng');
  figure
  hold on
  plot(s(:,1), s(:,4), 'b;kernel;');
  plot(s(:,1), s(:,5), 'r;gradient;');
  title([ktitle ' and gradient']);
  xlabel('x-Coord');
  ylabel('W(x), d W(x) / d x');
  print([fname '.all.eps'], '-deps2c');
  print([fname '.all.png'], '-dpng');
