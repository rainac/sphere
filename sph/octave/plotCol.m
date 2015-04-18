% -*- octave -*-
% Johannes Willkomm 2008
% 
function res = plotCol(pfun, fname, dcol, pcol, title)
  data = load(fname);
  res = gcf;
  ind=data(:, 1);
  pfun(ind, data(:,dcol), [pcol]);
