function CellParameters = GenCellPara(r,CelliNum,ColorPat)
% CellParameters = GenCellPara(r,CelliNum,ColorPat)
% cell radius
% ColorPat = [0.5,1,0.5;1,0,0];
CellNum = sum(CelliNum);
TypeList = ones(CellNum,1);
a = 1;
for n = 1:length(CelliNum)
  TypeList(a:a+CelliNum(n)-1) = n;
  a = a+CelliNum(n);
end
CellParameters.Color = ColorPat(TypeList,:);
CellParameters.CellType = TypeList;
CellParameters.Volume = 4*pi/3*r^3*ones(CellNum,1);