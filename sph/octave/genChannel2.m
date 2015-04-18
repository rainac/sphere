% -*- octave -*-
% Johannes Willkomm 2008
% generate 3-part channel: botton and two sides
function res = genChannel2(delta, length, depth, angle)
 
  s1 = genplane(delta, [length depth]);
  s2 = s1;

  s1 = rotate(s1, [1 0 0], 2*pi -angle/2);
  s2 = rotate(s2, [1 0 0], angle/2);

  s2 = translate(s2, [0 0 depth*sin(angle/2)]);

  res = s1;
  res = [s1' s2']';

