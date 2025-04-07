 function SimParameters = SimParaModify(DefaultParameters,CellIdx,CellParameters,ParaModify)
% ParaModify(1).ParaName = 'gamma';
% ParaModify(1).ParaValue = 0.1;
% ParaModify(1).Motif = {'EMS',0.25;'P2',0.25};
% 
% ParaModify(2).ParaName = 'sigma';
% ParaModify(2).ParaValue = [0.9,0.2];
% ParaModify(2).Motif = {'P2',0.2;{'ABp','ABa'},0.9};
SimParameters = DefaultParameters;
CellNum = length(CellIdx);
if CellNum > 1
	SimParameters.gamma = DefaultParameters.gamma*ones(CellNum,1);
	SimParameters.sigma = DefaultParameters.sigma*(ones(CellNum)-eye(CellNum));
end

if nargin > 2
	for n = 1:length(ParaModify)
		switch ParaModify(n).ParaName
      case 'M'
        SimParameters.M = VecPara(ParaModify(n).ParaValue,CellIdx,CellParameters,ParaModify(n).Motif);
			case 'gamma'
				SimParameters.gamma = VecPara(ParaModify(n).ParaValue,CellIdx,CellParameters,ParaModify(n).Motif);
			case 'c'
				SimParameters.c = ParaModify(n).ParaValue;
			case 'sigma'
        SimParameters.sigma = ArrayPara(ParaModify(n).ParaValue,CellIdx,CellParameters,ParaModify(n).Motif,ParaModify(n).Noise);
        if any(SimParameters.sigma<0,'all')
          numstr = num2str(SimParameters.sigma,'%.1f');
          sigmaStr = sprintf('%s\n',numstr(1,:));
          for m = 2:size(numstr,1)-1
            sigmaStr = cat(2,sigmaStr,sprintf('%s\n',numstr(m,:)));
          end
          sigmaStr = cat(2,sigmaStr,sprintf('%s',numstr(m+1,:)));
          warning('sigma=\n[%s]<0, make it positive\n',sigmaStr)
          SimParameters.sigma = abs(SimParameters.sigma);
        end
		end
	end
end
end

%% subfunctions
function Para = VecPara(ParaValue,CellIdx,CellParameters,Motif)
CellNum = length(CellIdx);
Para = ParaValue(1)*ones(CellNum,1);
if ~isempty(Motif)
	for n = 1:size(Motif,1)
		i = CellIdx == find(strcmp(CellParameters.Name,Motif{n,1}));
		Para(i) = Motif{n,2};
	end
end
end

function para = ArrayPara(ParaValue,CellIdx,CellParameters,Motif,Noise)
if size(ParaValue,1)>1
  % directly assign
  para = ParaValue;
else
  % ParaValue(1) for non-sister attractions; ParaValue(2) for sister attractions
  CellNum = length(CellIdx);
  para = ParaValue(1)*(ones(CellNum)-eye(CellNum));
  % binary: sister weak; non-sister strong
  if length(ParaValue) > 1
    for n = 1:CellNum
      SisterIdx = CellParameters.DaughtersIdx{CellParameters.MotherIdx{CellIdx(n)}};
      para(n,CellIdx == SisterIdx(SisterIdx~=CellIdx(n))) = ParaValue(2);
    end
  end
  % motif
  if ~isempty(Motif)
    for n = 1:size(Motif,1)
      if iscell(Motif{n,1})
        i = CellIdx == find(strcmp(CellParameters.Name,Motif{n,1}{1}));
        j = CellIdx == find(strcmp(CellParameters.Name,Motif{n,1}{2}));
        para(i,j) = Motif{n,2};
        para(j,i) = Motif{n,2};
      else
        i = CellIdx == find(strcmp(CellParameters.Name,Motif{n,1}));
        para(i,:) = Motif{n,2};
        para(:,i) = Motif{n,2};
        para(i,i) = 0;
      end
    end
  end
  % noise
  for n = 1:CellNum
    for m = n+1:CellNum
      r = 2*Noise*(rand-0.5);
      para(n,m) = para(n,m)+r;
      para(m,n) = para(m,n)+r;
    end
  end
end
end