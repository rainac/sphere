% -*- octave -*-
% Johannes Willkomm 2008
% 
function res = getCol2(version, runName, proc, comp, m, salgo, dcol)
  ftmp = [sprintf('../results/performance/%d', version) ...
          '/%s/%s/m%d/%s/' ...
          sprintf('%s-efficiency.txt', runName)];
  fname=sprintf(ftmp, proc, comp, m, salgo);
  data = load(fname);
  res = data(:,dcol);
  