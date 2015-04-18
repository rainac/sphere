% -*- matlab -*-
% Johannes Willkomm 2008
% load H5Part with matlab :(

function res = hdf5load(filename)
  hinfo = hdf5info(filename);
  
  res = struct();
  
  ngroups = length(hinfo.GroupHierarchy.Groups);
  
  for i=1:ngroups
    grp = hinfo.GroupHierarchy.Groups(i);

    ndatasets = length(grp.Datasets);
    grpres = struct();
    
    for j=1:ndatasets
      dset = grp.Datasets(j);
      dsetres = hdf5read(filename, dset.Name);
      name = substring(dset.Name, 1);
      [A, name] = strtok(name, '/');
      name = substring(name, 1);
      name = strrep(name, '(', '_');
      name = strrep(name, ')', '_');
      name = strrep(name, '/', '_');
      grpres.(name) = dsetres'; % transpose, since octave appears
                                % to do the same
    end
  
    [stepStr, number] = strtok(grp.Name, '#');
    groupIndexStr = substring(number, 1);
    res.(sprintf('Step_%s', groupIndexStr)) = grpres;
end

