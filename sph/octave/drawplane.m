% -*- octave -*-
% Johannes Willkomm 2008
function r = drawplane(delta, p1, p2, n)

  dv = p2 - p1;

  n = n ./ norm(n);

  r = genplane(delta, sizes);

  rotAxis = cross(n, [0, 0, 1])
  rotAngle = acos(norm(n*rotAxis') / norm(n)*norm(rotAxis))

  mitte = sizes ./ 2;
  mitte = [ mitte 0]
    
  if (norm(rotAxis))
    rotAxis = rotAxis ./ norm(rotAxis)

    qrot = quaternion(rotAxis, rotAngle);
    
    for i=1:size(r,1)
      v = r(i, :) - mitte;
      vsz = norm(v);
      if (vsz)
	q = quaternion(v ./ vsz, rotAngle);
	q = qmult(qmult(qinv(qrot), q), qrot);
	[v, theta] = quaternion(q);
	r(i, :) = v .* vsz + mitte;
      endif
    endfor
  endif

  for i=1:size(r,1)
    r(i, :) = r(i, :) + p1 - mitte;
  endfor

