function res = norm3(v)
   res = zeros(size(v,2), 1);
   for i=1:size(v,2)
      res(i) = norm(v(:,i));
   end
