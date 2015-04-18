% -*- octave -*-
% Johannes Willkomm 2008
% generate 3-part channel: botton and two sides
function res = genChannel2Capped(delta, length, depth, angle)

  ch = genChannel2(delta, length, depth, angle);

  cap = genplane(delta, [2*depth*sin(angle/2) 2*depth*sin(angle/2)]);
  cap = rotate(cap, [0 1 0], radians(90));
  cap = translate(cap, [-length/2 0 depth*sin(angle/2)/2], radians(90));

  res = [ch' cap']';
