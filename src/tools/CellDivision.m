function [DaughterCellPhi,DaughterCellIdx] = CellDivision(CellParameters,MotherCellIdx,phi,xx,yy,zz)
DivisionAxis = CellParameters.DivisionAxis{MotherCellIdx};
DaughterCellIdx = CellParameters.DaughtersIdx{MotherCellIdx};
ObjRatio = CellParameters.Volume(DaughterCellIdx(1))/CellParameters.Volume(DaughterCellIdx(2));
x_c = sum(phi(:).*xx(:))/sum(phi(:));
y_c = sum(phi(:).*yy(:))/sum(phi(:));
z_c = sum(phi(:).*zz(:))/sum(phi(:));
b = fzero(@(x)OptObj(ObjRatio,phi,DivisionAxis,x,xx,yy,zz,x_c,y_c,z_c),0);
[phi1,phi2] = Separation(phi,DivisionAxis,b,xx,yy,zz,x_c,y_c,z_c);
DaughterCellPhi = {phi1;phi2};
end

function [phi1,phi2] = Separation(phi,DivisionAxis,b,xx,yy,zz,x_c,y_c,z_c)
% div_surf > 0 -> phi1
% div_surf < 0 -> phi2
% DivisionAxis points from phi2 to phi1
div_surf = (xx-x_c)*DivisionAxis(1)+(yy-y_c)*DivisionAxis(2)+(zz-z_c)*DivisionAxis(3)+b;
phi1 = phi.*(tanh(div_surf)+1)/2;  
phi2 = phi.*(tanh(-div_surf)+1)/2; 
end

function ratio = OptObj(ObjRatio,phi,DivisionAxis,b,xx,yy,zz,x_c,y_c,z_c)
[phi1,phi2] = Separation(phi,DivisionAxis,b,xx,yy,zz,x_c,y_c,z_c);
ratio = sum(phi1(:))/sum(phi2(:))-ObjRatio;
end
