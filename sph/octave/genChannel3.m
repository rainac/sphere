% -*- octave -*-
% Johannes Willkomm 2008
% generate 3-part channel: botton and two sides
function res = genChannel3(delta, length, width, depth)
 
  boden = genplane(delta, [length width]);
  s1 = genplane(delta, [length depth]);

  s2 = translate(s1, [0 0 width]);

  boden = rotate(boden, [1 0 0], radians(90));
  s2 = translate(s1, [0 0 width]);

  res = [boden' s1' s2']';

