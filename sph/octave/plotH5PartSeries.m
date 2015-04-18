function r = plotH5PartSeries(fnames, prefix, extractField, extractFunction=@absmean)

  colors = {'r', 'g', 'b', 'k', 'y', 'm'};

  hold on;

  for i=1:length(fnames)
    fname = [prefix '/' fnames{i}];
    t = fnames{i};
    
    d = load(fname);
    series = extractH5PartSeries(d, extractField, extractFunction);
    
    semilogy(1:length(series), series, [colors{mod(i,length(colors))+1} ...
                        ';' t ';']);
  end

  hold off;
