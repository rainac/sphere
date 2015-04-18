% -*- octave -*-
% Johannes Willkomm 2008

function r = genSceneBox(delta, size, offset, boxLx, boxLy, boxLz)

  b = genBox(delta, boxLx, boxLy, boxLz);

  r = [ b']';

  r = scale(r, size);
  r = translate(r, offset);
