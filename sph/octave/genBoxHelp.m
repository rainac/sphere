% -*- octave -*-
% Johannes Willkomm 2008
function res = genBoxHelp(outf, offset, delta)
 
  len = 0.4;
  
  r = fillplane(delta, [offset(1), offset(2), offset(3)], [offset(1) + len, offset(2), offset(3) + len], 1, 3);
  dumpascii(r, outf);

  r = fillplane(delta, [offset(1), offset(2), offset(3)], [offset(1) + len, offset(2) + len, offset(3)], 1, 2);
  dumpascii(r, outf);

  r = fillplane(delta, [offset(1), offset(2), offset(3) + len], [offset(1) + len, offset(2) + len, offset(3) + len], 1, 2);
  dumpascii(r, outf);

  r = fillplane(delta, [offset(1), offset(2), offset(3)], [offset(1), offset(2) + len, offset(3) + len], 2, 3);
  dumpascii(r, outf);

  r = fillplane(delta, [offset(1) + len, offset(2), offset(3)], [offset(1) + len, offset(2) + len, offset(3) + len], 2, 3);
  dumpascii(r, outf);

