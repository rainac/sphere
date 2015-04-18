% -*- octave -*-
% Johannes Willkomm 2008
% 
function res = plotEfficiencies(fname)
  data = load(fname);
  res = gcf
  hold on
  ind=data(:, 1).*data(:, 2);
  loglog(ind, data(:,3), 'r*-;Efficiency;');
%  plot(ind, data(:,11), 'g+-;Efficiency update-index;');
%  plot(ind, data(:,14), 'bx-;Efficiency summation;');
%  plot(ind, data(:,17), 'ro-;Efficiency eval-state;');
%  plot(ind, data(:,20), 'g^-;Efficiency update-state;');
