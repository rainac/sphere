% -*- octave -*-
% Johannes Willkomm 2008
function p = schwerpunkt(data)

  sum = zeros(size(data(1, :)));

  for i=1:size(data,1)
    sum = sum + data(i, :);
  endfor

  p = sum / size(data,1)
