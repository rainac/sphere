% -*- octave -*-
% Johannes Willkomm 2008
function res = genb10(fname, delta)
 
  outf = fopen(fname, 'w');

  r = fillplane(delta, [0.1, 0.1, 0.1], [0.9, 0.1, 0.9], 1, 3);
  dumpascii(r, outf);
  
  r = fillplane(delta, [0.1, 0.1, 0.1], [0.1, 0.3, 0.9], 2, 3);
  dumpascii(r, outf);

  r = fillplane(delta, [0.9, 0.1, 0.1], [0.9, 0.3, 0.9], 2, 3);
  dumpascii(r, outf);

  r = fillplane(delta, [0.1, 0.1, 0.1], [0.9, 0.3, 0.1], 1, 2);
  dumpascii(r, outf);

  r = fillplane(delta, [0.3, 0.1, 0.9], [0.9, 0.3, 0.9], 1, 2);
  dumpascii(r, outf);

  fclose(outf);
