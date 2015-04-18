
fnames = { ...
'box1-v0-d5e-5.h5', ...
'box1-v0.0001-d5e-5.h5', ...
'box1-v0.01-d5e-5.h5', ...
'box1-v0.1-d5e-5.h5', ...
'box1-v1-d5e-5.h5', ...
'box1-v10-d5e-5.h5' };

extractField = 'ad_dichte_d0';
extractFunction = @absmean;

figure;
plotH5PartSeries(fnames, extractField, extractFunction);
title('mean of |ad\_dichte\_d0| over time steps');
xlabel('time step');
ylabel('mean of |ad\_dichte\_d0|');
print('-dpng', '-S1024,576', 'img/derivOverTime_box_dichte.png');

figure;
extractField = 'ad_coord_d0_0';
plotH5PartSeries(fnames, extractField, extractFunction);
title('mean of |ad\_coord\_d0\_0| over time steps');
xlabel('time step');
ylabel('mean of |ad\_coord\_d0\_0|');
print('-dpng', '-S1024,576', 'img/derivOverTime_box_xcoord.png');

figure;
extractField = 'ad_coord_d0_1';
plotH5PartSeries(fnames, extractField, extractFunction);
title('mean of |ad\_coord\_d0\_1| over time steps');
xlabel('time step');
ylabel('mean of |ad\_coord\_d0\_1|');
print('-dpng', '-S1024,576', 'img/derivOverTime_box_ycoord.png');

