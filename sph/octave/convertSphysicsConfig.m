% -*- matlab -*-
% Johannes Willkomm 2008
% Convert H5Part to individual VTK files
% Usage: h52vtk('res.h5', 'outdir')

function res = convertSphysicsConfig(name, outname, usead, minPos, maxPos)

  indatName = [name '/INDAT'];
  indatData = load(indatName);

  leFile = fopen(outname, 'w');
  
  fprintf(leFile, '# -*- sh -*-\n');
  fprintf(leFile, '# from SPHYSICS Case-Directory "%s"\n', name);

  fprintf(leFile, 'export SPH_PARAM_QUADER_P2="%.16g %.16g"\n', ...
          maxPos(1) + 0.01, maxPos(2) + 0.01);

  parB = indatData(12);
  fprintf(leFile, 'export SPH_PARAM_B="%.16g"\n', parB);
  
  parGamma = indatData(13);
  fprintf(leFile, 'export SPH_PARAM_GAMMA="%d"\n', parGamma);

  parRho0 = indatData(16);
  fprintf(leFile, 'export SPH_PARAM_RHO0="%.16g"\n', parRho0);

  parEpsilon = indatData(15);
  fprintf(leFile, 'export SPH_PARAM_EPSILON="%.16g"\n', parEpsilon);


  dx = indatData(22);
  dy = indatData(23);
  dz = indatData(24);
%  parH = sqrt(dx*dx + dy*dy + dz*dz);

  parH = indatData(25);
  fprintf(leFile, 'export SPH_PARAM_H="%.16g"\n', parH);
  fprintf(leFile, 'export SPH_PARAM_LJ_R0="%.16g"\n', parH/3);

  parNumSensors = indatData(29);
  fprintf(leFile, 'export SPH_PARAM_SENSORS="%d"\n', parNumSensors);

  sensorDataName = [name '/sensors.txt'];
  sensorData = load(sensorDataName);

  for i=0:parNumSensors-1
    fprintf(leFile, 'export SPH_PARAM_SENSOR_S%d_P1="%.16g %.16g %.16g"\n', ...
            i, sensorData(2*i + 1, 1) - minPos(1), sensorData(2*i + 2, 1) - minPos(2), 0);
    fprintf(leFile, 'export SPH_PARAM_SENSOR_S%d_P2="%.16g %.16g %.16g"\n', ...
            i, sensorData(2*i + 1, 2) - minPos(1), sensorData(2*i + 2, 2) - minPos(2), 0);
  end
  
  parDT = indatData(31);
  fprintf(leFile, 'export SPH_PARAM_DT="%.16g"\n', parDT);

  parTmax = indatData(32);
  fprintf(leFile, 'export SPH_PARAM_TMAX="%.16g"\n', parTmax);

  if (usead)
    fprintf(leFile, ['export SPH_PARAM_LOAD_FIELDS=' ...
                     '"masse dichte waerme coord vel class color ad_coord_d0 ad_coord_d1"\n']);
    fprintf(leFile, ['export SPH_PARAM_SAVE_FIELDS=' ...
                     '"masse dichte waerme druck coord vel class color' ...
                     ' id sumOfNeighbourIds neighbours' ...
                     ' ad_coord ad_vel ad_dichte ad_druck' ...
                     '"\n']);
  end
  
  fclose(leFile);

