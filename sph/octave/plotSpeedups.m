% -*- octave -*-
% Johannes Willkomm 2008
% 
function res = plotSpeedups(fname)
  data = load(fname);
  res = gcf;
  hold on
  ind=data(:, 1).*data(:, 2);
  loglog(ind, ind, 'k-;ideal Speedup;');
  loglog(ind, data(:,4), 'r*-;Speedup overall;');
%  loglog(ind, data(:,12), 'g+-;Speedup update-index;');
%  loglog(ind, data(:,15), 'bx-;Speedup summation;');
%  loglog(ind, data(:,18), 'ro-;Speedup eval-state;');
%  loglog(ind, data(:,21), 'g^-;Speedup update-state;');
