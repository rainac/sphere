% -*- octave -*-
% Johannes Willkomm 2008

function pipe = genUTube(delta, size, offset, pipeWidth, uWidth, heightLeft, heightRight)

  p1Length = heightLeft - pipeWidth;
  p2Length = uWidth - 2*pipeWidth;
  p3Length = heightRight - pipeWidth;

  p1 = genPipe(delta, p1Length, pipeWidth);
  p2 = genPipe(delta, p2Length, pipeWidth);
  p3 = genPipe(delta, p3Length, pipeWidth);

  c1 = genPipe90(delta, pipeWidth);
  c2 = genPipe90(delta, pipeWidth);

  p1 = rotate(p1, [0 0 1], pi/2);
  p3 = rotate(p3, [0 0 1], pi/2);

  c1 = rotate(c1, [0 0 1], -pi/2);

  p1 = translate(p1, [-p2Length/2 - pipeWidth/2,  p1Length/2 + pipeWidth/2,  0]);
  p3 = translate(p3, [ p2Length/2 + pipeWidth/2,  p3Length/2 + pipeWidth/2,  0]);

  c1 = translate(c1, [-p2Length/2 - pipeWidth/2, 0, 0]);
  c2 = translate(c2, [ p2Length/2 + pipeWidth/2, 0, 0]);

  pipe = [p1; c1; p2; c2; p3];

  pipe = scale(pipe, size);
  pipe = translate(pipe, offset);
