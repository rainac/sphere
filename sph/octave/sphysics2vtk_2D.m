% -*- matlab -*-
% Johannes Willkomm 2008
% Convert H5Part to individual VTK files
% Usage: h52vtk('res.h5', 'outdir')

function res = sphysics2vtk_2D(name, outname, binary, usead)

  ipartName = [name '/IPART'];
  indatName = [name '/INDAT'];
  
  ipartData = load(ipartName);
%  whos ipartData
  
  indatData = load(indatName);
%  whos indatName

  nParticles      = indatData(26);
  nBoundary       = indatData(27);
  nBoundaryFixed  = indatData(28);

  if size(ipartData, 1) ~= nParticles
    error(['np read from "%s" is not the same as dimension of matrix ' ...
           ' in "%s": %d %d'], indatName, ipartName, nParticles, size(ipartData, 1))
  end
  
  shiftPos = min(ipartData(:,1:2));
  if (shiftPos(1) < 0)
    ipartData(:,1) = ipartData(:,1) - shiftPos(1);
  end
  if (shiftPos(2) < 0)
    ipartData(:,2) = ipartData(:,2) - shiftPos(2);
  end

  minPos = min(ipartData(:,1:2));
  maxPos = max(ipartData(:,1:2));

  names = {'coord_0', 'coord_1', 'vel_0', 'vel_1', 'dichte', 'druck', 'masse', 'vorty'};
 
  stepdata = struct();
  for i=1:length(names)
    stepdata.(names{i}) = ipartData(:,i)';
  end
  
  classVector = zeros(1, nParticles);
  classVector(1, 1:nBoundary) = 1;
  classVector(1, nBoundary+1:nParticles) = 0;
  stepdata.('class') = classVector;
  
  stepdata.('color') = zeros(1, nParticles);
  stepdata.('waerme') = zeros(1, nParticles);

  if (usead)
    ipartDiffName = [name '/IPART_D'];
%    fprintf(1, 'found derivatives in ''%s''\n', ipartDiffName);
    if exist(ipartDiffName, 'file')
      ipartDiffData = load(ipartDiffName);
      if size(ipartDiffData, 1) ~= nParticles
        error(['np read from "%s" is not the same as dimension of matrix ' ...
               ' in "%s": %d %d'], indatName, ipartDiffName, nParticles, size(ipartDiffData, 1))
      end
      
      namesDiff = {'ad_coord_d0_0', 'ad_coord_d0_1', 'ad_coord_d1_0', 'ad_coord_d1_1'};
      for i=1:length(namesDiff)
        stepdata.(namesDiff{i}) = ipartDiffData(:,i)';
      end
    end

    ipartDiffName = [name '/IPART_L'];
    if exist(ipartDiffName, 'file')
%      fprintf(1, 'found derivatives in ''%s''\n', ipartDiffName);
      ipartDiffData = load(ipartDiffName);
      if size(ipartDiffData, 1) ~= nParticles
        error(['np read from "%s" is not the same as dimension of matrix ' ...
               ' in "%s": %d %d'], indatName, ipartDiffName, nParticles, size(ipartDiffData, 1))
      end
      
      namesDiff = {'ad_coord_d0_0', 'ad_coord_d0_1'};
      for i=1:length(namesDiff)
        stepdata.(namesDiff{i}) = ipartDiffData(:,i)';
      end
      
      ipartDiffName = [name '/IPART_S'];
      ipartDiffData = load(ipartDiffName);
      if size(ipartDiffData, 1) ~= nParticles
        error(['np read from "%s" is not the same as dimension of matrix ' ...
               ' in "%s": %d %d'], indatName, ipartDiffName, nParticles, size(ipartDiffData, 1))
      end
      
      namesDiff = {'ad_coord_d1_0', 'ad_coord_d1_1'};
      for i=1:length(namesDiff)
        stepdata.(namesDiff{i}) = ipartDiffData(:,i)';
      end
    end
  end

  indizes = 1:nParticles;
  saveVTK(stepdata, indizes, outname, binary);

  convertSphysicsConfig(name, [outname '.env'], usead, shiftPos, maxPos);
  