function res = dnorm3(v, dv)
   res = zeros(size(v,2), 1);
   for i=1:size(v,2)
      res(i) = (v(:,i)'*dv(:,i)) / norm(v(:,i));
   end
