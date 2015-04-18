% -*- octave -*-
% Johannes Willkomm 2008
function num = numPrefixes(c1, c2)
  num = 0;
  for i=1:length(c1)
    found = max(strncmp(c1{i}, c2, length(c1{i})));
    num = num + found;
  end
end
