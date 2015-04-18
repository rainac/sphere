function m=plotadres(fnames, column)
  markers={'+', '*', 'o', 'x', '^'};
  colors={'r', 'g', 'b'};
  for i=1:length(fnames)
    pcolor=colors{mod(i,length(colors))+1};
    pmark=markers{mod(i,length(markers))+1};
    fname=fnames{i}
    data=load(fname);
    plotFunctionAndDerivative(data(:,1), data(:,column), ...
    data(:,column + 2));
  end
