% -*- octave -*-
% Johannes Willkomm 2008

function r = genSceneContainerWith3Slides(delta, size, offset, channelTheta=120, channelNeigung=2)

  channels = genScene3Slides(delta, 1, [0.35 0 0.42], channelTheta, channelNeigung);
  
  b = genBoxLochUnten(delta, 0.35, 0.35, 0.45, 0.03, 0.03);
  b = translate(b, [0.5 0.65 0.5]);

  r = [channels' b']';

  r = scale(r, size);
  r = translate(r, offset);
