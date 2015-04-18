% -*- matlab -*-
% Johannes Willkomm 2008

function res = extractH5PartSeries(data, fieldName, selectFun=@max)
  res = zeros(length(fieldnames(data)), 1);
  for i=1:length(fieldnames(data))
    stepName = sprintf('Step_%d', i-1);
    stepData = data.(stepName);
    fieldData = stepData.(fieldName);
    res(i, 1) = selectFun(fieldData);
  end
