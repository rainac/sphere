
close all;

d1 = load "results/performance/1556/sun/gcc/hr/m6/symmetric/run-3D-hr-raw.txt";
d2 = load "results/performance/1801/ultraT2/gcc/hr/m6/symmetric/run-3D-hr-raw.txt";
d3 = load "results/performance/1962/ultraT2/gcc/hr/m6/symmetric/run-3D-hr-raw.txt";

loglog(d1(:, 1), d1(:,2), 'g-x;Version 1556(Monaghan);')
hold on
loglog(d2(:, 1), d2(:,2), 'r-+;Version 1801(Ferrari);')
loglog(d3(:, 1), d3(:,2), 'b-*;Version 1962(Ferrari);')
xlabel('Threads');
ylabel('time/s');

print('-S640,480', '-dpng', 'img/vcomp_1556_1801_1962.png');
