% -*- octave -*-
% Johannes Willkomm 2008

function r = genSceneContainerWithPipe2(delta, size, offset)

  pipeWidth = 0.07;

  b = genBoxLochUnten(delta, 0.45, 0.45, 0.9, [pipeWidth pipeWidth], [0.16 0]);
  pipe = genUTube(delta, 1, [0 0 0], pipeWidth, 0.35, 0.1, 0.15);

  pipe = translate(pipe, [0.3, -0.45 - 0.1 + pipeWidth/2, 0]);

  r = [b; pipe];
  r = scale(r, size);
  r = translate(r, offset - min(r));
