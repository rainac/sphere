% -*- octave -*-
% Johannes Willkomm 2008

function r = genSceneArtesianWellWithBottom(delta, size, offset)

  well = genSceneContainerWithPipe(delta, 1, [0 0 0]);
  
  bottom = genBox(delta, 0.9, 0.9, 0.1);
  bottom = translate(bottom, [0.4 0 0.1]);

  r = [well' bottom']';

  r = scale(r, size);
  r = translate(r, offset - min(r));
