% -*- octave -*-
% Johannes Willkomm 2008

function r = genSceneBoxTriang(delta, size, offset, box, triangleBox, trianglePos)

  b = genBox(delta, box(1), box(2), box(3));

  width = triangleBox(1);
  height = triangleBox(2);
  w_2 = width / 2;
  wallLength = sqrt(height*height + w_2*w_2);
  slope = asin(height/wallLength);

  p1 = genplane(delta, [wallLength, box(3)], [1, 3]);
  p2 = p1;
  
  p1 = rotate(p1, [0, 0, 1], slope);
  p2 = rotate(p2, [0, 0, 1], -slope);

  p1 = translate(p1, [(trianglePos(1)-box(1)/2 + w_2/2)        -box(2)/2+triangleBox(2)/2  trianglePos(2)]);
  p2 = translate(p2, [(trianglePos(1)-box(1)/2 + w_2/2 + w_2)  -box(2)/2+triangleBox(2)/2  trianglePos(2)]);

  r = [ b' p1' p2']';

  r = scale(r, size);
  r = translate(r, offset);
