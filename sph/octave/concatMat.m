% -*- matlab -*-
% Johannes Willkomm 2008

function alldata = concatMat(flist)

  len = length(flist)
  
  data = load(flist{1});
  
  alldata = data.d;
  whos alldata
  
  for i=2:len
    data = load(flist{i});
    alldata = [alldata; data.d];
    whos alldata
  end
