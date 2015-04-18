% -*- matlab -*-
% Johannes Willkomm 2008
% Convert H5Part to individual VTK files
% Usage: h52vtk('res.h5', 'res%04d.vtk')

function stem = filename_stem(name)
  infileStemEndPos = find(name == '.', 1, 'last');
  stem = name(1:infileStemEndPos - 1);
 
