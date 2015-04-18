% -*- matlab -*-
% Johannes Willkomm 2008
% Save a struct of vectors to a VTK file (Version 3.0 binary format)
% function res = saveH5StepBinary(stepdata, indizes, outname, binary)
%    stepdata: struct with vectors of dimension 1xN
%    indizes:  an integer vector of indices which vector components
%              to save
%    outname:  name of VTK outfile, without .vtk suffix
%    binary:   use VTK binary or ASCII format

function res = saveVTK(stepdata, indizes, outname, binary, bBox)

  coordTemplate  = 'coord_%d';

  vectorFieldNames = {'vel', 'accel', 'ad_coord_d0', 'ad_coord_d1', 'ad_coord_d2', 'ad_vel_d0', 'ad_vel_d1', 'ad_vel_d2'};

  outname = [outname '.vtk'];
  
  [outfile, fopenMsg] = fopen(outname, 'w');
  if outfile == -1 
    error(['failed to open file ''' outname ''': ' fopenMsg]);
  end

  numPoints = length(indizes);
  stepfieldnams = fieldnames(stepdata);
  numFields = length(stepfieldnams);
  numVectorFieldsPresent = numPrefixes(vectorFieldNames, stepfieldnams);

  fprintf(outfile, '# vtk DataFile Version 3.0\n');
  if (binary)
    fprintf(outfile, 'BINARY\n');
  else
    fprintf(outfile, 'ASCII\n');
  end
  fprintf(outfile, 'DATASET UNSTRUCTURED_GRID\n');
  
  % save the point coordinates
  fprintf(outfile, 'POINTS %d double\n', numPoints);
  coordField1 = sprintf(coordTemplate, 0);
  coordField2 = sprintf(coordTemplate, 1);
  coordField3 = sprintf(coordTemplate, 2);
  coords1 = stepdata.(coordField1);
  coords1 = coords1(indizes);
  numDims = 1;
  if (isfield(stepdata, coordField2))
    coords2 = stepdata.(coordField2);
    coords2 = coords2(indizes);
    numDims = 2;
  else
    coords2 = zeros(size(coords1));
  end
  if (isfield(stepdata, coordField3))
    coords3 = stepdata.(coordField3);
    coords3 = coords3(indizes);
    numDims = 3;
  else
    coords3 = zeros(size(coords1));
  end
  
  if (size(size(coords1)) ~= 2)
    error('error: coords1: expecting a vector')
  end
  if (size(coords1,1) ~= 1 && size(coords1,2) ~= 1)
    error('error: coords1: expecting a vector')
  end
  if (size(coords1, 2) == 1)
    fprintf(1, ['warning: transposing vector coords1\n'])
    coords1 = coords1';
  end
  if (size(coords2, 2) == 1)
    fprintf(1, ['warning: transposing vector coords2\n'])
    coords2 = coords2';
  end
  if (size(coords3, 2) == 1)
    fprintf(1, ['warning: transposing vector coords3\n'])
    coords3 = coords3';
  end

  allCoords = reshape([coords1; coords2; coords3], numPoints*3, 1);
  for i=1:min(100, numPoints)
    if (~ isequal(coords1(i), allCoords((i-1)*3 + 1)))
      error('error: value in coords1(%d) should be the sames as allCoords(%d): %g %g', ...
            i, (i-1)*3 + 1, coords1(i), allCoords((i-1)*3 + 1))
    end
    if (~ isequal(coords2(i), allCoords((i-1)*3 + 2)))
      error('error: value in coords2(%d) sould be the sames as allCoords(%d): %g %g', ...
            i, (i-1)*3 + 2, coords1(i), allCoords((i-1)*3 + 1))
    end
    if (~ isequal(coords3(i), allCoords((i-1)*3 + 3)))
      error('error: value in coords2(%d) sould be the sames as allCoords(%d): %g %g', ...
            i, (i-1)*3 + 3, coords1(i), allCoords((i-1)*3 + 1))
    end
  end

  writeVector(outfile, allCoords, binary, 'double', 'ieee-be');
  fprintf(outfile, '\n');

  % specify the cells: here 1 coord <-> 1 cell
  fprintf(outfile, 'CELLS %d %d\n', numPoints, numPoints*2);
  cells1 = ones(numPoints, 1);
  cells2 = 0:numPoints-1;
  cellData = zeros(numPoints*2, 1, 'int32');
  cellData(:) = reshape([cells1'; cells2], numPoints*2, 1);
  writeVector(outfile, cellData, binary, 'int32', 'ieee-be');
  fprintf(outfile, '\n');

  % specify the cell types: all are type 1
  fprintf(outfile, 'CELL_TYPES %d\n', numPoints);
  celltypes = ones(numPoints, 1, 'int32');
  writeVector(outfile, celltypes, binary, 'int32', 'ieee-be');
  fprintf(outfile, '\n');

  % save remaining point data fields
  fprintf(outfile, 'CELL_DATA %d\n', numPoints);
  fprintf(outfile, 'POINT_DATA %d\n', numPoints);
  numFields = length(stepfieldnams) - numVectorFieldsPresent*(numDims-1) - numDims;
  fprintf(outfile, 'FIELD FieldData %d\n', numFields);
  fprintf(1, '%s: %d additional fields: ', outname, numFields);
  nField = 0;
  for j=1:length(stepfieldnams)
    datafield = stepfieldnams{j};
    if (strncmp(datafield, coordTemplate, 6))
      continue
    end
    datafieldBegin = datafield(1:max(strfind(datafield, '_'))-1);
    if (max(strcmp(datafieldBegin, vectorFieldNames))) 
      if (strcmp(datafield, [datafieldBegin '_0']) == 0)
        continue
      end
      if (nField > 0)  
        fprintf(1, ', ', datafieldBegin);
      end
      fprintf(1, '%s[]', datafieldBegin);
      fprintf(outfile, '%s 3 %d double\n', datafieldBegin, numPoints);
      %	fprintf(outfile, 'LOOKUP_TABLE default\n');
      vX = stepdata.(datafield);
      vX = vX(indizes);
      vecField2 = [datafieldBegin '_1'];
      if (isfield(stepdata, vecField2))
        vY = stepdata.(vecField2);
        vY = vY(indizes);
      else
        vY = zeros(size(vX));
      end
      vecField3 = [datafieldBegin '_2'];
      if (isfield(stepdata, vecField3))
        vZ = stepdata.([datafieldBegin '_2']);
        vZ = vZ(indizes);
      else
        vZ = zeros(size(vX));
      end
      if (size(vX, 2) == 1)
        fprintf(1, 'warning: transposing vector "%s"\n', datafield)
        vX = vX';
      end
      if (size(vY, 2) == 1)
        fprintf(1, 'warning: transposing vector "%s"\n', vecField2)
        vY = vY';
      end
      if (size(vZ, 2) == 1)
        fprintf(1, 'warning: transposing vector "%s"\n', vecField3)
        vZ = vZ';
      end
      vecData = [vX; vY; vZ];
      vecData = reshape(vecData, numPoints*3, 1);
      writeVector(outfile, vecData, binary, 'double', 'ieee-be');
      fprintf(outfile, '\n');
    else
      if (nField > 0)  
        fprintf(1, ', ', datafieldBegin);
      end
      fprintf(1, '%s', datafield);
      fprintf(outfile, '%s 1 %d double\n', datafield, numPoints);
      v = stepdata.(datafield);
      v = v(indizes);
      writeVector(outfile, v, binary, 'double', 'ieee-be');
      fprintf(outfile, '\n');
    end
    nField = nField + 1;
  end
  fprintf(1, '\r');

  if (nField ~= numFields)
    error('error: wrote %d fields, but thought would write %d', ...
          nField, numFields)
  end

  fclose(outfile);
