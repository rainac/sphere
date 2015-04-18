% -*- octave -*-
% Johannes Willkomm 2008
% 

function res=plotPerfResult(pfun, version, runName, m, comp, algo, ialgo, ...
                            proc, coli, plotSpec, titleSpec)
  
  ftmp = [sprintf('results/performance/%d', version) ...
          '/%s/%s/m%d/%s/%s/' ...
          sprintf('%s-efficiency.txt', runName)];
  fname=sprintf(ftmp, proc, comp, m, algo, ialgo);
  plotCol(pfun, fname, coli, plotSpec, titleSpec);
