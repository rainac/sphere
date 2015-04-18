fnames = {'res-2e-5-0.01.txt', 'res-2e-5-0.02.txt', 'res-2e-5-0.03.txt', ...
          'res-4e-5-0.01.txt', 'res-4e-5-0.02.txt', 'res-4e-5-0.03.txt', ...
          'res-8e-5-0.02.txt', 'res-8e-5-0.03.txt'};

fnames = {'res-2e-5-0.01.txt'};

fnames = {'res-0_10.txt'};
fnames = {'res-1.txt'};
fnames = {'res-2.txt'};

fnames = {'res-pres-mon.txt'};
fnames = {'res-vel-mon.txt'};
fnames = {'res-dens-mon.txt'};
fnames = {'res-v100-vel-fer.txt'};

figure; hold on;
plotadres(fnames, 4);

figure; hold on;
plotadres(fnames, 8);
