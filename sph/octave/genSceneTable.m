% -*- octave -*-
% Johannes Willkomm 2008

function r = genSceneTable(delta, size, offset, boxL, boxW)


  b1 = genplane(delta, [boxW, boxL], [1 3]);
  b2 = translate(b1, [0 -delta 0]);

  r = [ b1' b2' ]';

  r = scale(r, size);
  r = translate(r, offset);
