% -*- octave -*-
% Johannes Willkomm 2008
function res = genplaneWithHole(delta, sizes, axes, hole, holepos)

  if (nargin < 5)
    holepos=[0 0]
  end

  b = genplane(delta, sizes, axes);

  ix = find(b(:,axes(1)) < holepos(1)-hole(1)/2 | b(:,axes(1)) > holepos(1)+hole(1)/2);
  iy = find(b(:,axes(2)) < holepos(2)-hole(2)/2 | b(:,axes(2)) > holepos(2)+hole(2)/2);

  res = b(unique([ix' iy']), :);
