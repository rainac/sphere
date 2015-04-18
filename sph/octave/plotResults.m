% -*- octave -*-
% Johannes Willkomm 2008
% 
function res = plotResults(fnames, titles, pTitle='test', pPref='', ...
                           colors={'r', 'g', 'b', 'k', 'y'})
  figure
  hold on
  title([pTitle ' Times']);
  xlabel('P');
  ylabel('Time/s');
  for i=1:length(fnames)
    plotCol(fnames{i}, 5, colors{i}, ['Time ' titles{i}]);
  end
  print('-dpng', '-S640,480', [pPref 'times.png']);
  
  figure
  hold on
  title([pTitle ' Speedups']);
  xlabel('P');
  ylabel('Speedup');
  for i=1:length(fnames)
    plotCol(fnames{i}, 4, colors{i}, ['Speedup ' titles{i}]);
  end
  print('-dpng', '-S640,480', [pPref 'speedup.png']);

  figure
  hold on
  title([pTitle ' Efficiencies']);
  xlabel('P');
  ylabel('Efficiency');
  for i=1:length(fnames)
    plotCol(fnames{i}, 3, colors{i}, ['Efficiency ' titles{i}]);
  end
  print('-dpng', '-S640,480', [pPref 'efficiency.png']);
