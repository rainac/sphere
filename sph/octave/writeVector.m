% Johannes Willkomm 2008
function wr = writeVector(outfile, vector, binary, dataType, machineFormat)
  if (binary) 
    wr = fwrite(outfile, vector, dataType, 0, machineFormat);
  else
    numItems = 6; % per line
    for k=1:length(vector)
      fprintf(outfile, '%.16g ', vector(k));
      if (mod(k, numItems) == 0)
        fprintf(outfile, '\n');
      end
    end
    wr = length(vector);
  end
  
