% -*- octave -*-
% Johannes Willkomm 2008
% 

versions = [2177];
runName = "run";

ms = [6,7];
compiler = {"icpc"};
%algos = {"symmetric", "sym_naive", "naive"};
%algoNames = {"symmetric", "sym-naive", "naive"};
algos = {"sym-naive", "naive"};
%algos = {"sym-mb-tb"};
algoNames = {"symmetric", "naive"};
procs = {"amd64", "ultraT2"};
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

tripels = { ...
           {"icpc", "intel64"} ; ...
           {"icpc", "amd64"} ; ...
           {"gcc", "ultraT2"} ; ...
};

timeCol = 19;

alldata = zeros(9, 4);
legends = cell(6, 1);
for i=1:size(tripels,1)
  trip=tripels{i,:};
  threadsNaive = getCol(versions, runName, trip{2}, trip{1}, ...
                        6, "naive", 1);
  alldata(1:length(threadsNaive), 1) = threadsNaive;
  timesNaive6 = getCol(versions, runName, trip{2}, trip{1}, ...
                       6, "naive", timeCol);
  timesNaive7 = getCol(versions, runName, trip{2}, trip{1}, ...
                       7, "naive", timeCol);
  timesSym6 = getCol(versions, runName, trip{2}, trip{1}, ...
                     6, "sym-naive", timeCol);
  timesSym7 = getCol(versions, runName, trip{2}, trip{1}, ...
                  7, "sym-naive", timeCol);
%  semilogx(threadsNaive, timesNaive6 ./ timesSym6, ...
%           sprintf("rx-;%s, %s, naive(m=6) over sym(m=6);", trip{2}, trip{1}));
%  semilogx(threadsNaive, timesNaive7 ./ timesSym7, "g+-;naive, m=7 over sym, m=7;");
%  semilogx(threadsNaive, timesNaive7 ./ timesSym6, ...
%           sprintf("bo-;%s, %s, naive(m=7) over sym(m=6);", trip{2}, trip{1}));
  alldata(1:length(timesNaive7), 2*(i-1)+2) = timesNaive6 ./ timesSym6;
  legends{2*(i-1)+1} = ...
      sprintf("%s, %s, naive(m=6) over sym(m=6)", trip{2}, trip{1});
  alldata(1:length(timesNaive7), 2*(i-1)+3) = timesNaive7 ./ timesSym6;
  legends{2*(i-1)+2} = ...
      sprintf("%s, %s, naive(m=7) over sym(m=6)", trip{2}, trip{1});
end
f1=figure;
hold on;
title(["algorithm speedup, " colNames{timeCol}]);
h1=bar(1:9, alldata(:,2:7), ...
       sprintf("bo-;%s, %s, naive(m=7) over sym(m=6);", trip{2}, trip{1}), 0.8, "grouped");
legends
legend(legends);
set(gca, "ylim", [0, 5])
set(gca, "key", "on")
set(gca, "xtick", 1:9);
labels = cell(9, 1);
for i=1:length(alldata(:,1))
  labels{i} = sprintf("%d", alldata(i,1));
end
set(gca, "xticklabel", labels)
set(h1(1), "facecolor", "red")
set(h1(1), "edgecolor", [0.3, 0.3, 0.3])
set(h1(1), "linewidth", 2)
set(h1(2), "facecolor", [0.75, 0, 0])
set(h1(2), "edgecolor", [0, 0, 0])
set(h1(2), "linewidth", 2)
set(h1(3), "facecolor", "green")
set(h1(4), "facecolor", [0, 0.75, 0])
set(h1(5), "facecolor", "blue")
set(h1(6), "facecolor", [0, 0, 0.75])
ylabel('Ratio of naive/symmetric');
xlabel('Threads');

f2=figure;
hold on;
title(["times, " colNames{timeCol}]);
for i=1:size(tripels,1)
  trip=tripels{i,:};
  threadsNaive = getCol(versions, runName, trip{2}, trip{1}, ...
                        6, "naive", 1);
  timesNaive6 = getCol(versions, runName, trip{2}, trip{1}, ...
                       6, "naive", timeCol);
  timesNaive7 = getCol(versions, runName, trip{2}, trip{1}, ...
                       7, "naive", timeCol);
  timesSym6 = getCol(versions, runName, trip{2}, trip{1}, ...
                     6, "sym-naive", timeCol);
  timesSym7 = getCol(versions, runName, trip{2}, trip{1}, ...
                  7, "sym-naive", timeCol);
  loglog(threadsNaive, timesNaive6, ...
         sprintf("r*-;%s, %s, naive, m=6;", trip{2}, trip{1}));
  loglog(threadsNaive, timesNaive7, ...
         sprintf("gx-;%s, %s, naive, m=7;", trip{2}, trip{1}));
  loglog(threadsNaive, timesSym6, ...
         sprintf("b+-;%s, %s, sym, m=6;", trip{2}, trip{1}));
  loglog(threadsNaive, timesSym7, ...
         sprintf("ko-;%s, %s, sym, m=7;", trip{2}, trip{1}));
end
set(gca, "ylim", [1, 1000])
ylabel('Time / s');
xlabel('Threads');

timeCol=timeCol-1;
f3=figure;
hold on;
title([colNames{timeCol}]);
for i=1:size(tripels,1)
  trip=tripels{i,:};
  threadsNaive = getCol(versions, runName, trip{2}, trip{1}, ...
                        6, "naive", 1);
  timesNaive6 = getCol(versions, runName, trip{2}, trip{1}, ...
                       6, "naive", timeCol);
  timesNaive7 = getCol(versions, runName, trip{2}, trip{1}, ...
                       7, "naive", timeCol);
  timesSym6 = getCol(versions, runName, trip{2}, trip{1}, ...
                     6, "sym-naive", timeCol);
  timesSym7 = getCol(versions, runName, trip{2}, trip{1}, ...
                  7, "sym-naive", timeCol);
  semilogx(threadsNaive, timesNaive6, ...
           sprintf("r*-;%s, %s, naive, m=6;", trip{2}, trip{1}));
  semilogx(threadsNaive, timesNaive7, ...
           sprintf("gx-;%s, %s, naive, m=7;", trip{2}, trip{1}));
  semilogx(threadsNaive, timesSym6, ...
           sprintf("b+-;%s, %s, sym, m=6;", trip{2}, trip{1}));
  semilogx(threadsNaive, timesSym7, ...
           sprintf("ko-;%s, %s, sym, m=7;", trip{2}, trip{1}));
end
ylabel('Speedup');
xlabel('Threads');
%set(gca, "ylim", [1, 1000])
