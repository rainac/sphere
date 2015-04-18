% -*- octave -*-
% Johannes Willkomm 2008
function res = genb12(fname, delta)
 
  outf = fopen(fname, 'w');

  genBox(outf, [0.35, 0.35, 0.35], delta);

  fclose(outf);
