function [AdjecMatrix, AdjecMatrix_n] = IsContact(PhiStack,threshold)
CellNum = length(PhiStack);
% AdjecMatrix = false(CellNum);
AdjecMatrix_n = zeros(CellNum);
SE = strel('sphere',1);
CellMask = cell(CellNum,1);
CellRange = zeros(CellNum,6);
N = size(PhiStack{1});
for n = 1:CellNum
  CellMask{n} = imdilate(PhiStack{n}>threshold,SE);
  ind = find(CellMask{n});
  [y,x,z] = ind2sub(N,ind);
  sub = [y,x,z];
  CellRange(n,:) = [min(sub),max(sub)];
end
for n = 1:CellNum
  for m = n+1:CellNum
    lower = max(CellRange([n,m],1:3));
    upper = min(CellRange([n,m],4:6));
    if all(upper>=lower)
      InterfaceRange = {lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)};
      AdjecMatrix_n(m,n) = sum(CellMask{n}(InterfaceRange{:}).*CellMask{m}(InterfaceRange{:}),'all');
      AdjecMatrix_n(n,m) = AdjecMatrix_n(m,n);
%       AdjecMatrix(m,n) = any(CellMask{n}(InterfaceRange{:}).*CellMask{m}(InterfaceRange{:}),'all');
%       AdjecMatrix(n,m) = AdjecMatrix(m,n);
    end
  end
end
AdjecMatrix = AdjecMatrix_n>0;