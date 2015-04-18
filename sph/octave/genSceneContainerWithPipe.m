% -*- octave -*-
% Johannes Willkomm 2008

function r = genSceneContainerWithPipe(delta, size, offset)

  b = genBoxLochRechts(delta, 0.2, 0.2, 0.9, 0.05, 0.05);
%  b = translate(b, [0.1 0.55 0.5]);

  p1 = genPipe(delta, 0.3, 0.05);
  p2 = genPipe90(delta, 0.05);
  p3 = genPipe(delta, 0.1, 0.05);

  p3 = rotate(p3, [0 0 1], pi/2);

  p1 = translate(p1, [0.25 -0.425 0]);

  p2 = translate(p2, [0.425 -0.425 0]);

  p3 = translate(p3, [0.425 -0.35 0]);

  r = [b' p1' p2' p3']';

  r = scale(r, size);
  r = translate(r, offset - min(r));
