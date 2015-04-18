% -*- octave -*-
% Johannes Willkomm 2008
% generate 3-part channel: botton and two sides
function res = genBoxLochUnten(delta, width, length, depth, holeSize, holePos)
 
  if (nargin < 5)
    holepos=[0 0]
  end

  s1 = genplane(delta, [length depth]);
  s2 = genplane(delta, [width depth], [3 2]);
  s3 = s1;
  s4 = s2;

  s1 = translate(s1, [0 0 -width/2]);
  s3 = translate(s3, [0 0 width/2]);

  s2 = translate(s2, [-length/2 0 0]);
  s4 = translate(s4, [length/2  0 0]);

  b = genplaneWithHole(delta, [length width], [1 3], holeSize, holePos);
  b = translate(b, [0 -depth/2 0]);

  res = s1;
  res = [s1' s2' s3' s4' b']';

