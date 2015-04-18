% -*- matlab -*-
% Johannes Willkomm 2008

function stem = filename_dir(name)
  fileDirEndPos = find(name == '/', 1, 'last');
  stem = substring(name, 0, fileDirEndPos - 2);
