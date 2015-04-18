% -*- octave -*-
% Johannes Willkomm 2008
function res = genb11(fname, delta)
 
  outf = fopen(fname, 'w');

#   r = fillplane(delta, [0, 0.95, 0], [1, 0.95, 1], 1, 3);
#   dumpascii(r, outf);

  offset = [0.4, 0.4, 0];
  genTreppe(outf, offset, delta);
  
  offset = [0.6, 0.2, 0];
  genTreppe(outf, offset, delta);

#   offset = [0.4, 0.4, 0];
#   genTreppe(outf, offset, delta);

#   offset = [0.6, 0.2, 0];
#   genTreppe(outf, offset, delta);

  fclose(outf);
