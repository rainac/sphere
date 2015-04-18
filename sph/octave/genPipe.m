% -*- octave -*-
% Johannes Willkomm 2008
% generate 3-part channel: botton and two sides
function res = genPipe(delta, length, depth)
 
  s1 = genplane(delta, [length depth]);
  s2 = genplane(delta, [length depth], [1 3]);
  s3 = s1;
  s4 = s2;

  s1 = translate(s1, [0  0  -depth/2]);
  s3 = translate(s3, [0  0  depth/2 ]);

  s2 = translate(s2, [0  -depth/2 0]);
  s4 = translate(s4, [0  depth/2  0]);

  res = [s1' s2' s3' s4']';

  