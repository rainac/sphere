% -*- octave -*-
% Johannes Willkomm 2008
% 
function res = plotTimes(fname)
  data = load(fname);
  res = gcf
  hold on
  ind=data(:, 1).*data(:, 2);
  loglog(ind, data(:,5), 'r*-;Time;');
%  loglog(ind, data(:,13), 'g+-;Time update-index;');
%  loglog(ind, data(:,16), 'bx-;Time summation;');
%  loglog(ind, data(:,19), 'ro-;Time eval-state;');
%  loglog(ind, data(:,22), 'g^-;Time update-state;');
