% -*- octave -*-
% Johannes Willkomm 2008
% 

function res=plotPerfResults(pfun, versions, runName, ms, compiler, algos, ...
                             algoNames, procs, oss, coli, colors, markers)
  
ms = [6,7];
i=1;
hold on
for version=versions
  for pi=1:length(procs)
    proc=procs{pi};
    for ci=1:length(compiler)
      comp=compiler{ci};
      if strcmp(comp, "icpc") && strcmp(proc, "ultraT2")
        continue
      end
      if (strcmp(comp, "studio") && strcmp(proc, "amd64")) ...
            || (strcmp(comp, "studio") && strcmp(proc, "intel64"))
        continue
      end
      for m=ms
        for ai=1:length(algos)
          col=colors{mod(i,length(colors)) + 1};
          mark=markers{mod(i,length(markers)) + 1};
          plotSpec=[col mark '-'];
          plotPerfResult(pfun, version, runName, m, comp, algos{ai}, ...
                         algoNames{ai}, proc, coli, plotSpec);
          i=i+1;
        end
      end
    end
  end
end
