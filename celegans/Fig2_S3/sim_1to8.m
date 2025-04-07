load('../CellParameters.mat','CellParameters')
DefaultSimParameters.tau = 2.62;
DefaultSimParameters.M = 8;
DefaultSimParameters.gamma = 0.25;
DefaultSimParameters.g = 1.6;
DefaultSimParameters.g_shell = 16;
DefaultSimParameters.sigma = 0;
DefaultSimParameters.c = 18;
%% eggshell
eggshellConfig = {'axisLength',[27.7846,18.3778,14.7022],'compress',10.1188,'epsilon',0.1};%'loadEggshell',path
%% setting
dx = 0.25;
dt = 0.1;
S1 = 12;
S2 = 0;

sigma1 = 0.5;
sigma2 = 0;
%% StageEnd
% 1-cell
StageEnd(1).DividingCell = {};
StageEnd(1).SimParaModify = [];
StageEnd(1).StageEndDef = {'SS',1e-4};
% 2-cell
StageEnd(2).DividingCell = {'P0'};
StageEnd(2).SimParaModify = [];
StageEnd(2).StageEndDef = {'SS',1e-4};
% 3-cell
StageEnd(3).DividingCell = {'AB'};
StageEnd(3).SimParaModify = struct('ParaName','sigma','ParaValue',[sigma1,sigma2],'Motif',{{}},'Noise',0);
StageEnd(3).StageEndDef = {'SS',1e-4};
% 4-cell with attraction
StageEnd(4).DividingCell = {'P1'};
StageEnd(4).SimParaModify = struct('ParaName','sigma','ParaValue',sigma1,'Motif',{{{'EMS','P2'},sigma2}},'Noise',0);
StageEnd(4).StageEndDef = {'SS',1e-4};
% 6-cell
StageEnd(5).DividingCell = {'ABa';'ABp'};
StageEnd(5).SimParaModify = struct('ParaName','sigma','ParaValue',[sigma1,sigma2],'Motif',{{}},'Noise',0);
StageEnd(5).StageEndDef = {'SS',1e-4};
% 7-cell
StageEnd(6).DividingCell = {'EMS'};
StageEnd(6).SimParaModify = struct('ParaName','sigma','ParaValue',[sigma1,sigma2],'Motif',{{}},'Noise',0);
StageEnd(6).StageEndDef = {'SS', 1e-4;
                           'QSS', 1};
% 8-cell
StageEnd(7).DividingCell = {'P2'};
StageEnd(7).SimParaModify = struct('ParaName','sigma','ParaValue',[sigma1,sigma2],'Motif',{{{'ABpl','E'},sigma2}},'Noise',0);
StageEnd(7).StageEndDef = {'SS', 1e-4;
                           'QSS', 1};
% end
StageEnd(end+1).DividingCell = {};
StageEnd(end).SimParaModify = [];
StageEnd(end).StageEndDef = {'SS',1e-5};
N = ceil([60,40,30]/dx);
L = N*dx;
tspan = [0,50000];
SimNum = 1;
for SimIdx = 1:SimNum
  SimName = sprintf('1to8dx%sdt%sadh_1',dot2d(dx),dot2d(dt));
  SimStageEnd = StageEnd;
  mainBDF2VC(SimName,L,N,dt,DefaultSimParameters,CellParameters,SimStageEnd,tspan,...
    'EggshellConfig',eggshellConfig,...
    'StabilizationPara',[S1,S2],...
    'GPU',...
    'SaveTrajectory',...
    'SavePhi',{'step',1000});
end
