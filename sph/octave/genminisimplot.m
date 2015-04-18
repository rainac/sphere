%
%

titleSuffix = sprintf('Y=%g, Visc=%s, %s', ...
  Y, visc, particle);

testname = sprintf('Y%g_DT%g_VS%s_%s', Y, dt, visc, particle)

indata = [ 'src/minisim_' testname '.h5' ];

load(indata);

basedir = ['img/minisim/' testname ];
mkdir(basedir)

T=(1:size(p1.density, 2))*1e-5;

addpath('octave/minisimplot')

msp1
msp2
msp3
msp4
msp5
msp6

close all
