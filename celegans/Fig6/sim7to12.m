addpath('/bicmr/home/kuangxy/bin')
load('CellParameters.mat','CellParameters')
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
% 4-cell
StageEnd(4).DividingCell = {'P1'};
StageEnd(4).SimParaModify = struct('ParaName','sigma','ParaValue',sigma1,'Motif',{{{'EMS','P2'},sigma2}},'Noise',0);
StageEnd(4).StageEndDef = {'Duration',12140};
% 6-cell
StageEnd(5).DividingCell = {'ABa';'ABp'};
StageEnd(5).SimParaModify = struct('ParaName','sigma','ParaValue',[sigma1,sigma2],'Motif',{{}},'Noise',0);
StageEnd(5).StageEndDef = {'Duration',16420};
% 7-cell
StageEnd(6).DividingCell = {'EMS'};
StageEnd(6).SimParaModify = struct('ParaName','sigma','ParaValue',[sigma1,sigma2],'Motif',{{}},'Noise',0);
StageEnd(6).StageEndDef = {'Duration',360};
% 8-cell
StageEnd(7).DividingCell = {'P2'};
StageEnd(7).SimParaModify = struct('ParaName','sigma','ParaValue',[sigma1,sigma2],'Motif',{{{'ABpl','E'},sigma2}},'Noise',0);
StageEnd(7).StageEndDef = {'Duration',1340};
% 12-cell
StageEnd(8).DividingCell = {'ABal';'ABar';'ABpl';'ABpr'};
StageEnd(8).SimParaModify = struct('ParaName','sigma','ParaValue',[sigma1,sigma2],'Motif',{{{'MS','E'},sigma1}},'Noise',0);
StageEnd(8).StageEndDef = {'Duration',10000};
% end
StageEnd(end+1).DividingCell = {};
StageEnd(end).SimParaModify = [];
StageEnd(end).StageEndDef = {'Duration',10000};

mut.ParaName = 'DivisionAxis';
% mut.Mutation = {'EMS',[0,-1,0]}
mut.Mutation = {'ABar',CellParameters.DivisionAxis{24}};
CellParameters = CellParaModify(CellParameters,mut);

N = ceil([60,40,30]/dx);
L = N*dx;
% tspan = [23360,50000];
tspan = [24700,50000];
SimNum = 1;%length(gammaABpl);
for SimIdx = 1:SimNum
  % SimName = 'div1';
  SimName = 'div2';
  SimStageEnd = StageEnd;
  mainBDF2VC(SimName,L,N,dt,DefaultSimParameters,CellParameters,SimStageEnd,tspan,...
    'EggshellConfig',eggshellConfig,...
    'StabilizationPara',[S1,S2],...
    'GPU',...
    'SaveTrajectory',...
    'LoadDir','1to8dx0d25dt0d1adh_1',...
		'EnforceStage',7,...
    'SavePhi',{'step',500});
	cd(SimName)
  AutoImg('dx',dx,'ImgType',3,'CellPara',CellParameters,'save','close','fig')
  cd ..
end
