% -*- octave -*-
% Johannes Willkomm 2008

function r = genStair(delta, stairLength, stairHeight, stairWidth, ...
                      wallWidth, odd)

  base = genplane(delta, [stairLength stairWidth], [1 3]);
  down = genplane(delta, [stairHeight stairWidth], [2 3]);
  down = translate(down, [stairLength*0.5 -stairHeight*0.5 0]);
  
  up = [];
  
  if odd
    up = genplane(delta, [stairHeight stairWidth], [2 3]);
    up = translate(up, [stairLength*0.5 stairHeight*0.5 0]);
  end
  
  r = [base; down; up];
  
