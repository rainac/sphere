% -*- octave -*-
% Johannes Willkomm 2008
% 

function res=plotResIndexTimesBar(versions, runName, timeCol, titelText, ialgo)
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
             {"gcc", "intel64"} ; ...
             {"gcc", "amd64"} ; ...
             {"gcc", "ultraT2"} ; ...
             };

  threadsUltra = getCol(versions, runName, "ultraT2", "gcc", ...
                        6, "naive", ialgo, 1);

  numTests = size(threadsUltra,1);

  alldata = zeros(numTests, 7);
  legends = cell(6, 1);
  for i=1:size(tripels,1)
    trip=tripels{i,:};
    threadsNaive = getCol(versions, runName, trip{2}, trip{1}, ...
                          6, "naive", ialgo, 1);
    alldata(1:length(threadsNaive), 1) = threadsNaive;
%    timesNaive6 = getCol(versions, runName, trip{2}, trip{1}, ...
%                         6, "naive", "csort", timeCol);
    timesNaive7 = getCol(versions, runName, trip{2}, trip{1}, ...
                         7, "naive", ialgo, timeCol);
    timesSym6 = getCol(versions, runName, trip{2}, trip{1}, ...
                       6, "symmetric", ialgo, timeCol);
%    timesSym7 = getCol(versions, runName, trip{2}, trip{1}, ...
%                       7, "sym-naive", "csort", timeCol);
                                %  semilogx(threadsNaive, timesNaive6 ./ timesSym6, ...
%    alldata(1:length(timesNaive7), 2*(i-1)+2) = timesNaive6;
    alldata(1:length(timesNaive7), 2*(i-1)+2) = timesNaive7;
    alldata(1:length(timesNaive7), 2*(i-1)+3) = timesSym6;
%    alldata(1:length(timesNaive7), 2*(i-1)+5) = timesSym7;
%    legends{2*(i-1)+1} = sprintf("%s, %s, naive(m=6)", trip{2}, trip{1});
    legends{2*(i-1)+1} = sprintf("%s, %s, naive(m=7)", trip{2}, trip{1});
    legends{2*(i-1)+2} = sprintf("%s, %s, symmetric(m=6)", trip{2}, trip{1});
%    legends{2*(i-1)+4} = sprintf("%s, %s, sym(m=7)", trip{2}, trip{1});
  end
  res=figure;
  hold on;
%  title(["algorithm speedup, " colNames{timeCol}]);
  title(titelText);
  h1=bar(1:numTests, alldata(:,2:7), 
         sprintf("bo-;%s, %s, naive(m=7) over sym(m=6);", trip{2}, trip{1}), 0.8, "grouped");
%  legends
  legend(legends);
  set(gca, "ylim", [0.1, 100])
  set(gca, "yscale", "log")
  set(gca, "key", "on")
  set(gca, "xtick", 1:numTests);
  labels = cell(numTests, 1);
  for i=1:numTests
    labels{i} = sprintf("%d", alldata(i,1));
  end
  set(gca, "xticklabel", labels)
  color=[1, 0, 0];
  for i=0:2
    set(h1(2*i+1), "facecolor", color);
%    set(h1(2*i+2), "facecolor", 0.85*color);
    set(h1(2*i+2), "facecolor", 0.7*color);
%    set(h1(2*i+4), "facecolor", 0.55*color);
    color=shift(color, 1);
  end
  ylabel('Time / s');
  xlabel('Threads');

