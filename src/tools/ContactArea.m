function [ContArea,MemArea] = ContactArea(PhiStack,dr,th,CellIdx,CellList,CellParameters)
if nargin<3
  th = 0.34;
end
CellNum = length(PhiStack);
N = size(PhiStack{1});
L = N([2,1,3])*dr;
x = GridInit(L,N([2,1,3]));
fv = cellfun(@(phi)isosurface(x{1},x{2},x{3},phi,th),PhiStack,'UniformOutput',0);
ContArea = zeros(CellNum);
MemArea = zeros(CellNum,1);
for n = 1:CellNum
  verts = fv{n}.vertices;
  faces = fv{n}.faces;
  facesVec1 = verts(faces(:,2),:)-verts(faces(:,1),:);
  facesVec2 = verts(faces(:,3),:)-verts(faces(:,1),:);
  s = vecnorm( cross(facesVec1,facesVec2)/2 ,2,2);
  MemArea(n) = sum(s);
  vint = round(verts/dr+1);
  vind = sub2ind(N,vint(:,2),vint(:,1),vint(:,3));
  for m = [1:n-1,n+1:CellNum]
    vContIdx = find(PhiStack{m}(vind)>=0.2);
    fContIdx = any(ismember(faces,vContIdx),2);
    ContArea(n,m) = sum(s(fContIdx));
  end
end
ContArea = (ContArea+ContArea')/2;

% sorting according to CellList
if nargin>3
  if nargin<6
    CellParameters = GetCellPara;
  end
  CellIdxS = cellfun(@(x)find(ismember(CellParameters.Name,x)),CellList);
  [a,I] = intersect(CellIdx,CellIdxS);
  if length(a)~=CellNum
    warning('CellList mismatch\n')
  else
    MemArea = MemArea(I);
    ContArea = ContArea(I,I);
  end
end