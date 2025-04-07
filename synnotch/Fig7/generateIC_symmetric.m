clear
% celltypeNum = [200,200];
celltypeNum = [120,280];
CellNum = sum(celltypeNum);
% CellNum = 500;
% synNotch
L = 75*ones(1,3);
dr = 0.5;
r = 5;
deltaL = 8;
% old framework
% L = [60,40,30];
% dr = 0.5;
% r = 4;
% deltaL = 2;
rc = PositionInit(CellNum,L-deltaL,r)+deltaL/2;
phi_save = cell(CellNum,1);
[xx,yy,zz] = meshgrid(0:dr:L(1)-dr,0:dr:L(2)-dr,0:dr:L(3)-dr);
for n = 1:CellNum
  d = r-sqrt((xx-rc(n,1)).^2+(yy-rc(n,2)).^2+(zz-rc(n,3)).^2);
  phi_save{n} = (tanh(d)+1)/2;
end
CellIdx = (1:CellNum)';
% CellParameters = GenCellPara(r,celltypeNum,[0.5,1,0.5;1,0,0]);
CellParameters = GenCellPara(r,celltypeNum,[0.5,1,0.5;0,0,1]);
StageIdx=1;
%% cell position check
CellType = CellParameters.CellType;
rc_mean = mean(rc);
d2center = vecnorm(rc-rc_mean,2,2);
dtype = [mean(d2center(CellType==1)),mean(d2center(CellType==2))];
%% vis
viewVec = [-50, 30];
DomainRange = [zeros(1,3),L];
[x,y,z] = meshgrid(DomainRange([1,4]),DomainRange([2,5]),DomainRange([3,6]));
v = [x(:),y(:),z(:)];

% plotted = zeros(8,1);
% for n = 1:8
%   idx = find(sum(v(n,:)==v,2)==2);
%   for m = find(plotted(idx)<3)'
%     plot3(v([n,idx(m)],1),v([n,idx(m)],2),v([n,idx(m)],3),'k');
%   end
%   plotted(n) = 3;
%   plotted(idx) = plotted(idx)+1;
% end

savename = sprintf('N%dr%sdx%s_1_0.mat',CellNum,dot2d(r),dot2d(dr));
cellParamName = sprintf('CellParameters_%s_1.mat',sprintf('%d',celltypeNum));
save(savename,'phi_save','CellIdx','StageIdx','-v7.3')
save(cellParamName,'CellParameters')
