
prefix = '/home/jw642797/work/sph-results/batch/2286-linux/visctest';

% fnames = { ...
% '2e-5-0.01-v0-vel-fer/treppe1_0.01_61.0.h5', ...
% '2e-5-0.01-v1-vel-fer/treppe1_0.01_61.0.h5', ...
% '2e-5-0.01-v10-vel-fer/treppe1_0.01_61.0.h5', ...
% '2e-5-0.01-v100-vel-fer/treppe1_0.01_61.0.h5', ...
%  };

fnames = { ...
'treppe1_0.01_70_v0.001.h5', ...
'treppe1_0.01_70_v0.01.h5', ...
'treppe1_0.01_70_v0.1.h5', ...
'treppe1_0.01_70_v1.h5', ...
'treppe1_0.01_70_v10.h5', ...
'treppe1_0.01_70_v100.h5', ...
 };

extractField = 'ad_dichte_d0';
extractFunction = @absmean;

figure;
plotH5PartSeries(fnames, prefix, extractField, extractFunction);
title('mean of |ad\_dichte\_d0| over time steps');
xlabel('time step');
ylabel('mean of |ad\_dichte\_d0|');
print('-dpng', '-S1024,576', 'img/derivOverTime_wehr2_dichte.png');

figure;
extractField = 'ad_coord_d0_0';
plotH5PartSeries(fnames, prefix, extractField, extractFunction);
title('mean of |ad\_coord\_d0\_0| over time steps');
xlabel('time step');
ylabel('mean of |ad\_coord\_d0\_0|');
print('-dpng', '-S1024,576', 'img/derivOverTime_wehr2_xcoord.png');

figure;
extractField = 'ad_coord_d0_1';
plotH5PartSeries(fnames, prefix, extractField, extractFunction);
title('mean of |ad\_coord\_d0\_1| over time steps');
xlabel('time step');
ylabel('mean of |ad\_coord\_d0\_1|');
print('-dpng', '-S1024,576', 'img/derivOverTime_wehr2_ycoord.png');

