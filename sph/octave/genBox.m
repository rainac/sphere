% -*- octave -*-
% Johannes Willkomm 2008
% generate 3-part channel: botton and two sides
function res = genBox(delta, Lx, Ly, Lz)
 
  s1 = genplane(delta, [Lx Ly], [1 2]);
  s2 = genplane(delta, [Lz Ly], [3 2]);
  s3 = s1;
  s4 = s2;

  s1 = translate(s1, [0 0 -Lz/2]);
  s3 = translate(s3, [0 0 Lz/2]);

  s2 = translate(s2, [-Lx/2 0 0]);
  s4 = translate(s4, [Lx/2  0 0]);

  b = genplane(delta, [Lx, Lz], [1 3]); % boden
%  b = rotate(b, [1 0 0], pi/2);
  b = translate(b, [0 -Ly/2 0]);

  res = s1;
  res = [s1' s2' s3' s4' b']';
