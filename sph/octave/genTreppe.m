% -*- octave -*-
% Johannes Willkomm 2008
function res = genTreppe(outf, offset, delta)
 
  r = fillplane(delta, [offset(1), offset(2), 0.1], [offset(1) + 0.18, offset(2), 0.9], 1, 3);
  dumpascii(r, outf);

  r = fillplane(delta, [offset(1) + 0.18, offset(2) + 0.02, 0.1], [offset(1) ...
                      + 0.2, offset(2) + 0.02, 0.9], 1, 3);
  dumpascii(r, outf);

  r = fillplane(delta, [offset(1) + 0.18, offset(2), 0.1], ...
                [offset(1) + 0.18, offset(2) + 0.02, 0.9], 2, 3);
  dumpascii(r, outf);

  r = fillplane(delta, [offset(1) + 0.2, offset(2), 0.1], ...
                [offset(1) + 0.2, offset(2) + 0.02, 0.9], 2, 3);
  dumpascii(r, outf);


  % linke, rechte seite
  r = fillplane(delta, [offset(1), offset(2), 0.1], [offset(1) + ...
                      0.2, offset(2) + 0.1, 0.1], 1, 2);
  dumpascii(r, outf);
  r = fillplane(delta, [offset(1), offset(2), 0.9], [offset(1) + ...
                      0.2, offset(2) + 0.1, 0.9], 1, 2);
  dumpascii(r, outf);

  r = fillplane(delta, [offset(1), offset(2), 0.1], [offset(1), offset(2) + 0.2, 0.9], 2, 3);
  dumpascii(r, outf);
