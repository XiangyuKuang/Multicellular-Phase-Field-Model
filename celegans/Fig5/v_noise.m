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
StageEnd(3).SimParaModify = [];
StageEnd(3).StageEndDef = {'SS',1e-4};
% 4-cell
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
% 12-cell
StageEnd(8).DividingCell = {'ABal';'ABar';'ABpl';'ABpr'};
StageEnd(8).SimParaModify = struct('ParaName','sigma','ParaValue',[sigma1,sigma2],'Motif',{{{'MS','E'},sigma1}},'Noise',0);
StageEnd(8).StageEndDef = {'Duration',10000};
% end
StageEnd(end+1).DividingCell = {};
StageEnd(end).SimParaModify = [];
StageEnd(end).StageEndDef = {'Duration',10000};
N = ceil([60,40,30]/dx);
L = N*dx;
tspan = [34700,50000];
SimNum = 100;
noise = 1;
DefaultSimParameters.vnoise = noise;
for SimIdx = 76:SimNum
	rng(SimIdx)
  SimName = sprintf('vn%s_%d',dot2d(noise),SimIdx);
  SimStageEnd = StageEnd;
  mainBDF2VC_vnoise(SimName,L,N,dt,DefaultSimParameters,CellParameters,SimStageEnd,tspan,...
    'EggshellConfig',eggshellConfig,...
    'StabilizationPara',[S1,S2],...
    'GPU',...
    'SaveTrajectory',...
    'LoadDir','../Fig2_S3/8adhto12dx0d25dt0d1c18sigma0d5-0_3',...
		'EnforceStage',8,...
    'SavePhi',{'step',1000});
end
