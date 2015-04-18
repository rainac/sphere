function m = plotFunctionAndDerivative(xdata, ydata, deriv, derivLength=0.01)
  plot(xdata, ydata);
  for i=1:length(deriv)
    dvec = [1 deriv(i)];
    dvec = dvec .* (derivLength / norm(dvec));
    line([xdata(i), xdata(i) + dvec(1)], [ydata(i), ydata(i) + dvec(2)], "color", "red")
  end
  
  