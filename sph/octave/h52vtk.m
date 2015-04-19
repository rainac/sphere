% -*- matlab -*-
% Johannes Willkomm 2008
% Convert H5Part to individual VTK files
% Usage: h52vtk('res.h5', 'outdir')

function res = h52vtk(name, outdir, binary, matlab, boundingBox, saveFnkt)
  stepTemplate   = 'Step_%d';
  coordTemplate  = 'coord_%d';
  flagsFieldName = 'flags';
  classFieldName = 'class';
  boundaryBit    = 8;
  if (nargin < 3)
    binary         = 1;
  end
  if (nargin < 4)
    matlab         = 0;
  end
  if (nargin < 5)
    boundingBox    = [0 0 0 1 1 1]; % only needed for saveIfrit
  end
  if (nargin < 6)
    saveFnkt       = @saveVTK;
  end

  % fallunterscheidung: wie HDF5 laden?
  if (matlab == 0)
    % 1) octave: einfach load() verwenden... :)
    hdfStruct = load(name);
    if (size(fieldnames(hdfStruct)) == 0)
      res = false(1)
      return
    end
  else
    % 2) MATLAB: hier wird's umstaendlich... :(
    %   hdf5load habe ich selbst geschrieben
    hdfStruct = hdf5load(name);
  end

%  if (length(outdir) ~= 0)
%    outdir = [outdir '/'];
%  end
  
  if ~exist(outdir, 'dir')
    fprintf(1, 'creating output directory "%s"\n', outdir);
    mkdir(outdir);
  end
  
  if ~exist(outdir, 'dir')
    error('output directory "%s" does not exist\n', outdir);
  end

  if ~exist(name, 'file')
    error('input file "%s" does not exist\n', name);
  end

  [infileDir infileStem infileSuffix] = fileparts(name);
%  infileStem = filename_stem(filename_notdir(name));
  
  outdir = [outdir '/' infileStem];

  if ~exist(outdir, 'dir')
    fprintf(1, 'creating output directory "%s"\n', outdir);
    r = mkdir(outdir);
    if ~exist(outdir, 'dir')
      error('failed to create output directory "%s"\n', outdir);
    end
  end

  outnameTemplate = [outdir '/res_%05d'];
  
  coordNames = {'coord_0', 'coord_1', 'coord_2'};
  
  fieldnams = fieldnames(hdfStruct);

  numSteps = length(fieldnams);

  for step=0:numSteps-1
    outname = sprintf(outnameTemplate, step);
    field = sprintf(stepTemplate, step);
    stepdata = hdfStruct.(field);

    if (isfield(stepdata, flagsFieldName)) 
      % compute indizes of non-boundary particles
      flagsField = stepdata.(flagsFieldName);
      numAllPoints = length(flagsField);
      % das cast() hier ist nur wegen MATLAB :(
      indizes = find(not(bitand(cast(flagsField, 'double'), boundaryBit)));
      boundaryIndizes = find(bitand(cast(flagsField, 'double'), boundaryBit));
    else
%      fprintf(1, ['warning: the field "flags" was not saved apparently, using ' ...
%       'field "class" to separate boundary from non-boundary particles\n'])
      if (isfield(stepdata, 'coord_0')) 
        classField = stepdata.(classFieldName);
        indizes = find((classField < hex2dec('40') | classField == 0) ...
                       & classField ~= 1);
        boundaryIndizes = find(classField >= hex2dec('40') | classField == 1);
%        indizes = 1:length(stepdata.coord_0);
%        boundaryIndizes = [];
      else
        % no flagsField and no positions present: no particles left
        break;
      end
    end

    if (step == 0) 
      saveFnkt(stepdata, boundaryIndizes, [outdir '/boundary'], binary, boundingBox);
    end

    saveFnkt(stepdata, indizes, outname, binary, boundingBox);

  end
