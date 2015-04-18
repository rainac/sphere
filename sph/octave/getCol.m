% -*- octave -*-
% Johannes Willkomm 2008
% 
function res = getCol(version, runName, proc, comp, m, salgo, ialgo, dcol)
  ftmp = [sprintf('results/performance/%d', version) ...
          '/%s/%s/m%d/%s/%s/' ...
          sprintf('%s-efficiency.txt', runName)];
  fname=sprintf(ftmp, proc, comp, m, salgo, ialgo);
  data = load(fname);
  res = data(:,dcol);
  