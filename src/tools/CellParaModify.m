function CellParameters = CellParaModify(DefaultParameters,ParaModify)
% ParaModify(1).ParaName = 'DivisionAxis'
% ParaModify(1).Mutation = {'EMS',ones(1,3);'P2',ones(1,3)};

% ParaModify(2).ParaName = 'Volume'
% ParaModify(2).Mutation = {'EMS',1000};

CellParameters = DefaultParameters;
for n = 1:length(ParaModify)
  if iscell(eval(sprintf('DefaultParameters.%s',ParaModify(n).ParaName)))
    str = 'CellParameters.%s{%d}=%s;';
  else
    str = 'CellParameters.%s(%d)=%s;';
  end
	for m = 1:size(ParaModify(n).Mutation,1)
		CellIdx = find(strcmp(CellParameters.Name,ParaModify(n).Mutation{m,1}));
		eval(sprintf(str,ParaModify(n).ParaName,CellIdx,'ParaModify(n).Mutation{m,2}'));
	end
end