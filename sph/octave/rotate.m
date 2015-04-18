% -*- octave -*-
% Johannes Willkomm 2008
function data = rotate(data, rotAxis, rotAngle, mitte=[0 0 0])

%  mitte = (max(data) - min(data)) ./ 2;
%  mitte = schwerpunkt(data);

  if (norm(rotAxis))
    rotAxis = rotAxis ./ norm(rotAxis);

    qrot = quaternion(rotAxis, rotAngle);
    
    for i=1:size(data,1)
      v = data(i, :) - mitte;
      vsz = norm(v, 2);
      if (vsz)
	v = v ./ vsz;
	q = quaternion(v(1), v(2), v(3), 0);
	q = qtrans(q, qrot);
	[v, theta] = quaternion(q);
	data(i, :) = v .* vsz + mitte;
      endif
    endfor
  else
    printf('warning: rotation axis is 0')
  endif
