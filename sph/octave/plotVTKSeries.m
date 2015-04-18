% -*- matlab -*-
% Johannes Willkomm 2008

function stem = plotVTKSeries(namePat, steps, pId, compIds)

  if (steps(2) < steps(1))
    error('upper bound lower then lower bound')
  end
  
  len = steps(2) - steps(1) + 1;
  
  alldata = zeros(len, length(compIds));
  
  for i=steps(1):steps(2)
    fname = sprintf(namePat, i)
    fdata = load(fname);
    for j=length(compIds)
      alldata(i-steps(1)+1,j) = fdata(pId, compIds(j));
    end
  end
  
  plotStyles = {'-+', '-*', '-x', '-o'};
  length(compIds)
  compIds
  figure;
  hold on
  for i=1:length(compIds)
    ctitle = sprintf(';comp: %d;', compIds(i));
    plot(steps(1):steps(2), alldata(:,i), [plotStyles{i} ctitle]);
  end
