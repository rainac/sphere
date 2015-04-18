% -*- octave -*-
% Johannes Willkomm 2008

function r = genSceneStairs(delta, size, offset, numStairs, angle, ...
                            stairLength, stairHeight, stairWidth, wallWidth, ...
                            sceneSize)

  r = [];

  stairOffset = zeros(1, 3);

  for i=1:numStairs
    b = genStair(delta, stairLength, stairHeight, stairWidth, wallWidth, mod(i,2) == 0);
    b = translate(b, stairOffset);
    r = [ r; b];
    stairOffset = stairOffset + [stairLength -stairHeight 0];
  end
  
  r = translate(r, [stairLength 0 stairWidth]*0.5);
  r = rotate(r, [0 0 1], d2r(angle));

  bleft = genplane(delta, [sceneSize(1) sceneSize(2)], [1 2]);
  bright = translate(bleft, [0 0 sceneSize(3)]);
  bbottom = genplane(delta, [sceneSize(1) sceneSize(3)], [1 3]);
  bbottom = translate(bbottom, [0 -sceneSize(2)*0.5 sceneSize(3)*0.5]);

  box = [];
  box = [bleft; bright; bbottom];
  box = translate(box, [sceneSize(1) sceneSize(2) 0]*0.5);

  r = [r; box];

  r = scale(r, size);
  r = translate(r, offset);

end

function r = d2r(a)
  r = (a / 360) * 2 * pi;
end
