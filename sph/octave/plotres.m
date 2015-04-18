                                % -*- octave -*-
% Johannes Willkomm 2008
% 

versions = [2199];
runName = "run";

%plotSize = '-S1728,972';
plotSize = '-S1296,729';
%plotSize = '-S864,486'; 

%plotResIndexTimesBar(versions(1), runName, 16, "Index update time Counting Sort", "csort")
%print('-dpng', plotSize, 'img/index-time-csort.png');

%plotResIndexTimesBar(versions(1), runName, 16, "Index update time GNU std::sort", "stdsort_par_index")
%print('-dpng', plotSize, 'img/index-time-stdsort.png');

plotResSpeedupBar(versions(1), runName, 3, "Total speedup")
print('-dpng', plotSize, sprintf('img/total-speedup-%d.png', versions));
print('-depsc2', sprintf("img/total-speedup-%d.eps", versions));

plotResSpeedupBar(versions(1), runName, 18, "Summation speedup")
print('-dpng', plotSize, sprintf('img/summation-speedup-%d.png', versions));
print('-depsc2', sprintf("img/summation-speedup-%d.eps", versions));

plotResSpeedupBar(versions(1), runName, 15, "Index update speedup")
print('-dpng', plotSize, sprintf('img/index-summation-speedup-%d.png', versions));

%return

plotResRatiosBar(versions(1), runName, 4, "Algorithm speedup, total time");
print('-dpng', plotSize, sprintf('img/algo-speedup-total-%d.png', versions));
print('-depsc2', sprintf("img/algo-speedup-total-%d.eps", versions));

plotResRatiosBar(versions(1), runName, 19, "Algorithm speedup, summation time");
print('-dpng', plotSize, sprintf('img/algo-speedup-summation-%d.png', versions));
print('-depsc2', sprintf("img/algo-speedup-summation-%d.eps", versions));

plotResRatiosBar(versions(1), runName, 22, "Algorithm speedup, eval time");
print('-dpng', plotSize, sprintf('img/algo-speedup-eval-%d.png', versions));
print('-depsc2', sprintf("img/algo-speedup-eval-%d.eps", versions));

%return

plotResEfficiencyBar(versions(1), runName, 2, "Total run time efficiency")
print('-dpng', plotSize, sprintf('img/total-efficiency-%d.png', versions));
print('-depsc2', sprintf("img/total-efficiency-%d.eps", versions));

plotResEfficiencyBar(versions(1), runName, 17, "Summation time efficiency")
print('-dpng', plotSize, sprintf('img/summation-efficiency-%d.png', versions));
print('-depsc2', sprintf("img/summation-efficiency-%d.eps", versions));

plotResTimesBar(versions(1), runName, 4, "Total run time")
print('-dpng', plotSize, sprintf('img/total-time-%d.png', versions));
print('-depsc2', sprintf("img/total-time-%d.eps", versions));

plotResTimesBar(versions(1), runName, 19, "Summation time")
print('-dpng', plotSize, sprintf('img/summation-time-%d.png', versions));
print('-depsc2', sprintf("img/summation-time-%d.eps", versions));

return

ms = [6,7];
compiler = {"gcc", "icpc", "studio"};
%algos = {"symmetric", "sym_naive", "naive"};
%algoNames = {"symmetric", "sym-naive", "naive"};
algos = {"sym-naive", "naive"};
%algos = {"sym-mb-tb"};
algoNames = {"symmetric", "naive"};
procs = {"ultraT2", "amd64", "intel64"};
oss = {"sunos", "linux", "linux"};

colors = {'r', 'g', 'b', 'k'};
markers = {'+', 'o', '*', 'x', 'd'};

cols=1:27;
colNames = {"threads", ...
            "efficiency-time", "speedup-time", "time", ... 
            "efficiency-time-steps", "speedup-time-steps", "time-steps", ...
            "efficiency-steps", "speedup-steps", "steps", ...
            "memory", "rels", "dtests", ...
            "efficiency-update-index", "speedup-update-index", "update-index", ...
            "efficiency-summation", "speedup-summation", "summation", ...
            "efficiency-eval-state", "speedup-eval-state", "eval-state", ...
            "efficiency-update-state", "speedup-update-state", "update-state", ...
            "clear-scratch", "reduce-scratch"};

whichCols=1:length(cols);
whichCols=2:4;
%whichCols=14:16;
whichCols=17:19;
%whichCols=26:27;
%whichCols=20:22;
%whichCols=2:4;

for coli=whichCols
  figure;
  hold on;
  [tok,rem] = strtok(colNames{coli}, "-");
  plotfun=@loglog;
  if (strcmp(tok, "efficiency") || strcmp(tok, "speedup"))
    plotfun=@plot;
  end
  plotPerfResults(plotfun, versions, runName, ms, compiler, algos, ...
                  algoNames, procs, oss, cols(coli), colors, markers);
  title([runName ' ' colNames{coli}]);
end

# figure;
# hold on;
# %algos = {"sym-naive"};
# plotPerfResults(@semilogx, version, runName, ms, compiler, algos, ...
#                 algoNames, procs, oss, 21, colors, markers);
# plotPerfResults(@semilogx, version, runName, ms, compiler, algos, ...
#                 algoNames, procs, oss, 22, colors, markers);
# title("Time U/R");

# figure;
# hold on;
# plotPerfResult(@semilogx, version, runName, 6, "icpc", "sym-naive", ...
#                 "symmetric", "amd64", 2, "r+-");
# plotPerfResult(@semilogx, version, runName, 7, "icpc", "naive", ...
#                 "naive", "amd64", 2, "g*-");
# plotPerfResult(@semilogx, version, runName, 6, "gcc", "sym-naive", ...
#                 "symmetric", "ultraT2", 2, "bx-");
# plotPerfResult(@semilogx, version, runName, 7, "gcc", "naive", ...
#                 "naive", "ultraT2", 2, "ko-");
# title("Efficiency");

# figure;
# hold on;
# plotPerfResult(@loglog, version, runName, 6, "icpc", "sym-naive", ...
#                 "symmetric", "amd64", 4, "r+-");
# plotPerfResult(@loglog, version, runName, 7, "icpc", "naive", ...
#                 "naive", "amd64", 4, "g*-");
# plotPerfResult(@loglog, version, runName, 6, "gcc", "sym-naive", ...
#                 "symmetric", "ultraT2", 4, "bx-");
# plotPerfResult(@loglog, version, runName, 7, "gcc", "naive", ...
#                 "naive", "ultraT2", 4, "ko-");
# title("Time");
