% -*- octave -*-
% Johannes Willkomm 2008
% generate 3-part channel: botton and two sides
function res = genPipe90(delta, depth)
 
  s1 = genplane(delta, [depth depth]);
  s2 = genplane(delta, [depth depth], [2 3]);
  s3 = s1;
  s4 = genplane(delta, [depth depth], [1 3]);

  s1 = translate(s1, [0  0  -depth/2]);
  s3 = translate(s3, [0  0  depth/2 ]);

  s2 = translate(s2, [depth/2 0  0]);
  s4 = translate(s4, [0  -depth/2 0]);

  res = [s1' s2' s3' s4']';
