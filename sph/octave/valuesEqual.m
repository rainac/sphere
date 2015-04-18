% -*- octave -*-
% Johannes Willkomm 2008
function eq = valuesEqual(v1, v2)
  if (v1 == v2)
    eq = true(1);
  else
    if (isnan(v1) && isnan(v2))
      eq = true(1);
    else
      eq = false(1);
    end
  end
