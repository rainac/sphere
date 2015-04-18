% -*- octave -*-
% Johannes Willkomm 2008
% 

versions = [2375];
runName = 'run';

textFontsize = 8;
legendFontsize = 8;

%csortName = 'counting_sort';
%stdsortName = 'std_sort_par_index';
%symAlgoName = 'sym-naive';

csortName = 'csort';
stdsortName = 'stdsort_par_index';
symAlgoName = 'symmetric';

paperPosition = [0 0 12 6.75];

figure;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', paperPosition);

hold on;
set(gca, 'fontsize', textFontsize);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
plotPerfResult(@loglog, versions, runName, ...
               6, 'gcc', symAlgoName, csortName, 'ultraT2', 16, 'rx-', 'M=64, Symmetric, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               7, 'gcc', symAlgoName, csortName, 'ultraT2', 16, 'g+-', 'M=128, Symmetric, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               6, 'gcc', symAlgoName, stdsortName, 'ultraT2', 16, 'bx-', 'M=64, Symmetric, GNU std::sort')
plotPerfResult(@loglog, versions, runName, ...
               7, 'gcc', symAlgoName, stdsortName, 'ultraT2', 16, 'k+-', 'M=128, Symmetric, GNU std::sort')
title('Index update times on SUN Niagara 2');
xlabel('Threads');
ylabel('Time (s)');
hl = legend('M=64, Symmetric, Counting Sort', ...
            'M=128, Symmetric, Counting Sort', ...
            'M=64, Symmetric, GNU std::sort', ...
            'M=128, Symmetric, GNU std::sort');
set(hl, 'fontsize', legendFontsize);
print('-depsc2', '-r1200', 'img/index-times-ultraT2-3d.eps');
print('-dpng', '-r600', 'img/index-times-ultraT2-3d.png');

figure;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', paperPosition);

hold on;
set(gca, 'fontsize', textFontsize);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
plotPerfResult(@loglog, versions, runName, ...
               6, 'gcc', symAlgoName, csortName, 'ultraT2', 14, 'rx-', 'M=64, Symmetric, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               7, 'gcc', symAlgoName, csortName, 'ultraT2', 14, 'g+-', 'M=128, Symmetric, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               6, 'gcc', symAlgoName, stdsortName, 'ultraT2', 14, 'bx-', 'M=64, Symmetric, GNU std::sort')
plotPerfResult(@loglog, versions, runName, ...
               7, 'gcc', symAlgoName, stdsortName, 'ultraT2', 14, 'k+-', 'M=128, Symmetric, GNU std::sort')
title('Index update efficiency on SUN Niagara 2');
xlabel('Threads');
ylabel('Efficiency');
hl = legend('M=64, Symmetric, Counting Sort', ...
            'M=128, Symmetric, Counting Sort', ...
            'M=64, Symmetric, GNU std::sort', ...
            'M=128, Symmetric, GNU std::sort', 'location', 'southwest');
set(hl, 'fontsize', legendFontsize);
print('-depsc2', '-r1200', 'img/index-efficiency-ultraT2-3d.eps');
print('-dpng', '-r600', 'img/index-efficiency-ultraT2-3d.png');

figure;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', paperPosition);

hold on;
set(gca, 'fontsize', textFontsize);
set(gca, 'xscale', 'log');
plotPerfResult(@semilogx, versions, runName, ...
               6, 'gcc', symAlgoName, csortName, 'ultraT2', 2, 'rx-', 'M=64, Symmetric, Counting Sort')
plotPerfResult(@semilogx, versions, runName, ...
               7, 'gcc', symAlgoName, csortName, 'ultraT2', 2, 'g+-', 'M=128, Symmetric, Counting Sort')
plotPerfResult(@semilogx, versions, runName, ...
               6, 'gcc', 'naive', csortName, 'ultraT2', 2, 'b*-', 'M=64, Naive, Counting Sort')
plotPerfResult(@semilogx, versions, runName, ...
               7, 'gcc', 'naive', csortName, 'ultraT2', 2, 'ko-', 'M=128, Naive, Counting Sort')
title('SPH efficiency on SUN Niagara 2');
xlabel('Threads');
ylabel('Efficiency');
hl = legend('M=64, Symmetric, Counting Sort', ...
            'M=128, Symmetric, Counting Sort', ...
            'M=64, Naive, Counting Sort', ...
            'M=128, Naive, Counting Sort', 'location', 'southwest');
set(hl, 'fontsize', legendFontsize);
print('-depsc2', '-r1200', 'img/efficiency-ultraT2-3d.eps');
print('-dpng', '-r600', 'img/efficiency-ultraT2-3d.png');


figure;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', paperPosition);

hold on;
set(gca, 'fontsize', textFontsize);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
plotPerfResult(@loglog, versions, runName, ...
               6, 'gcc', symAlgoName, csortName, 'ultraT2', 4, 'rx-', 'M=64, Symmetric, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               7, 'gcc', symAlgoName, csortName, 'ultraT2', 4, 'g+-', 'M=128, Symmetric, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               6, 'gcc', 'naive', csortName, 'ultraT2', 4, 'b*-', 'M=64, Naive, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               7, 'gcc', 'naive', csortName, 'ultraT2', 4, 'ko-', 'M=128, Naive, Counting Sort')
title('SPH times on SUN Niagara 2');
xlabel('Threads');
ylabel('Time (s)');
hl = legend('M=64, Symmetric, Counting Sort', ...
            'M=128, Symmetric, Counting Sort', ...
            'M=64, Naive, Counting Sort', ...
            'M=128, Naive, Counting Sort');
set(hl, 'fontsize', legendFontsize);
print('-depsc2', '-r1200', 'img/times-ultraT2-3d.eps');
print('-dpng', '-r600', 'img/times-ultraT2-3d.png');

close all
return

figure;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'paperposition', paperPosition);

hold on;
set(gca, 'fontsize', 16);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
plotPerfResult(@loglog, versions, runName, ...
               6, 'gcc', symAlgoName, csortName, 'ultraT2', 3, 'rx-', 'M=6, Symmetric, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               7, 'gcc', symAlgoName, csortName, 'ultraT2', 3, 'g+-', 'M=128, Symmetric, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               6, 'gcc', 'naive', csortName, 'ultraT2', 3, 'b*-', 'M=6, Naive, Counting Sort')
plotPerfResult(@loglog, versions, runName, ...
               7, 'gcc', 'naive', csortName, 'ultraT2', 3, 'ko-', 'M=128, Naive, Counting Sort')
title('SPH speedup on SUN Niagara 2');
xlabel('Threads');
ylabel('Speedup');
print('-depsc2', '-r1200', 'img/index-times-ultraT2-3d.eps');
print('-dpng', '-r600', 'img/speedup-ultraT2-3d.png');


ms = [6,7];
compiler = {'gcc', 'icpc', 'studio'};
%algos = {symAlgoName, 'sym_naive', 'naive'};
%algoNames = {symAlgoName, 'sym-naive', 'naive'};
algos = {'sym-naive', 'naive'};
%algos = {'sym-mb-tb'};
algoNames = {'symmetric', 'naive'};
procs = {'ultraT2', 'amd64', 'intel64'};
oss = {'sunos', 'linux', 'linux'};

colors = {'r', 'g', 'b', 'k'};
markers = {'+', 'o', '*', 'x', 'd'};

cols=1:27;
colNames = {'threads', ...
            'efficiency-time', 'speedup-time', 'time', ... 
            'efficiency-time-steps', 'speedup-time-steps', 'time-steps', ...
            'efficiency-steps', 'speedup-steps', 'steps', ...
            'memory', 'rels', 'dtests', ...
            'efficiency-update-index', 'speedup-update-index', 'update-index', ...
            'efficiency-summation', 'speedup-summation', 'summation', ...
            'efficiency-eval-state', 'speedup-eval-state', 'eval-state', ...
            'efficiency-update-state', 'speedup-update-state', 'update-state', ...
            'clear-scratch', 'reduce-scratch'};

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
  [tok,rem] = strtok(colNames{coli}, '-');
  plotfun=@loglog;
  if (strcmp(tok, 'efficiency') || strcmp(tok, 'speedup'))
    plotfun=@plot;
  end
  plotPerfResults(plotfun, versions, runName, ms, compiler, algos, ...
                  algoNames, procs, oss, cols(coli), colors, markers);
  title([runName ' ' colNames{coli}]);
end

% figure;
% hold on;
% %algos = {'sym-naive'};
% plotPerfResults(@semilogx, version, runName, ms, compiler, algos, ...
%                 algoNames, procs, oss, 21, colors, markers);
% plotPerfResults(@semilogx, version, runName, ms, compiler, algos, ...
%                 algoNames, procs, oss, 22, colors, markers);
% title('Time U/R');

% figure;
% hold on;
% plotPerfResult(@semilogx, version, runName, 6, 'icpc', 'sym-naive', ...
%                 'symmetric', 'amd64', 2, 'r+-');
% plotPerfResult(@semilogx, version, runName, 7, 'icpc', 'naive', ...
%                 'naive', 'amd64', 2, 'g*-');
% plotPerfResult(@semilogx, version, runName, 6, 'gcc', 'sym-naive', ...
%                 'symmetric', 'ultraT2', 2, 'bx-');
% plotPerfResult(@semilogx, version, runName, 7, 'gcc', 'naive', ...
%                 'naive', 'ultraT2', 2, 'ko-');
% title('Efficiency');

% figure;
% hold on;
% plotPerfResult(@loglog, version, runName, 6, 'icpc', 'sym-naive', ...
%                 'symmetric', 'amd64', 4, 'r+-');
% plotPerfResult(@loglog, version, runName, 7, 'icpc', 'naive', ...
%                 'naive', 'amd64', 4, 'g*-');
% plotPerfResult(@loglog, version, runName, 6, 'gcc', 'sym-naive', ...
%                 'symmetric', 'ultraT2', 4, 'bx-');
% plotPerfResult(@loglog, version, runName, 7, 'gcc', 'naive', ...
%                 'naive', 'ultraT2', 4, 'ko-');
% title('Time');
