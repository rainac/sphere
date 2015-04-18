% -*- octave -*-
% Johannes Willkomm 2008
function data = translate(data, v)
  for i=1:size(data,1)
    data(i, :) = data(i, :) +  v;
  endfor
