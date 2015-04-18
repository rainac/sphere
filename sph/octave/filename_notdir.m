% -*- matlab -*-
% Johannes Willkomm 2008

function stem = filename_notdir(name)
  fileDirEndPos = find(name == '/', 1, 'last');
  if (length(fileDirEndPos) == 0)
    stem = name;
  else
    stem = name(fileDirEndPos+1:length(name)-1);
  end
