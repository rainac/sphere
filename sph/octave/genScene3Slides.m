% -*- octave -*-
% Johannes Willkomm 2008

function r = genScene3Slides(delta, size, offset, channelTheta=120, channelNeigung=2)

  c1 = genChannel2Capped(delta, 0.5, 0.1, radians(channelTheta)); 
  c2 = genChannel2Capped(delta, 0.5, 0.1, radians(channelTheta)); 
  c3 = genChannel2Capped(delta, 0.5, 0.1, radians(channelTheta)); 

  c1 = rotate(c1, [0 0 1], radians(360-channelNeigung));
  c2 = rotate(c2, [0 0 1], radians(360-channelNeigung));
  c2 = rotate(c3, [0 0 1], radians(360-channelNeigung));
  
  c2 = rotate(c2, [0, 1, 0], 3*pi/2);

  c3 = rotate(c3, [0, 1, 0], pi);

  c1 = translate(c1, [0,   0.2,  0   ]);
  c2 = translate(c2, [0.3, 0.1,  0.2 ]);
  c3 = translate(c3, [0.1, 0,    0.5 ]);

  r = [c1', c2', c3']';

  r = scale(r, size);
  r = translate(r, offset - min(r));
