function [time, idx, FileList, SimName] = GetTime(varargin)
% Output: [time, idx, FileList, SimName]
temp = ismember(varargin,'IC');% bug caused by " "
EnableIC = any(temp);
DirNameIdx = [find(~temp,1),find(all(temp))];
switch length(varargin)
  case 0
    DirName = '.';
  case 1
    DirName = varargin{1};
    if EnableIC
      DirName = '.';
      warning('get IC time in current dir')
    end
  case 2
    DirName = varargin{DirNameIdx};
end
if isfolder(DirName)
% 	FileNameChar = ls(DirName);
  DirInfo = dir(DirName);
else
	time = [];
  idx = [];
  FileList = [];
  warning('dir %s do not exist',DirName)
	return
end
FileList = string({DirInfo.name}');
[timeVar,simNameVar] = cellfun(@StrFun,FileList,'UniformOutput',0);
time = cell2mat(timeVar);
% time = cell2mat(cellfun(@StrFun,FileList,'UniformOutput',0));
IsCellPara = cellfun(@(x)contains(x,'CellParameter'),FileList);
IsIC = cellfun(@(x)contains(x,'IC'),FileList);
IsMat = cellfun(@(x)contains(x,'.mat'),FileList);
if EnableIC
  IsPhi = IsMat&~IsCellPara&IsIC&~isnan(time);
else
  IsPhi = IsMat&~IsCellPara&~(isnan(time)|IsIC);
end
SimName = simNameVar(IsPhi);
[time,I] = sort(time(IsPhi));
idx = find(IsPhi);
idx = idx(I);
SimName = SimName(I);
USimName = unique(SimName);
if length(USimName)==1
  SimName = USimName{1};
else
  warning('multiple SimName in %s',DirName)
end
end

function [num,str] = StrFun(in)
idx1 = find(in=='_',1,'last');
idx2 = find(in=='.');
temp = in(idx1+1:idx2-1);
temp(temp == 'd') = '.';
str = in(1:idx1-1);
num = str2double(temp);
end
