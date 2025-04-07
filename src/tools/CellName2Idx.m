function [CellIdx,existIdx] = CellName2Idx(CellParameters,CellName)
% [CellIdx,existIdx] = CellName2Idx(CellParameters,CellName)
% output:
%   existIdx = 1 or 0 represents the existence or absence of a cell. 
%     CellName(existIdx) is a list of Existing Cells
%   CellIdx is the index for Existing Cells

if iscell(CellName)
  CellIdx = cellfun(@(x)find(ismember(CellParameters.Name,x)),CellName(:),'UniformOutput',0);
  existIdx = ~cellfun(@isempty,CellIdx);
  CellIdx = cell2mat(CellIdx(existIdx));
else
  CellIdx = find(ismember(CellParameters.Name,CellName));
  existIdx = ~isempty(CellIdx);
end