clear
DefaultSimParameters.tau = 2.62;
DefaultSimParameters.M = 8;
DefaultSimParameters.gamma = 0.25;
DefaultSimParameters.g = 1.6;
DefaultSimParameters.c = 2;
%% eggshell
AxisLength = [];

%% setting
dx = 0.5;
dt = 1.25;

%% StageEnd
StageEnd(1).DividingCell = {};
StageEnd(1).SimParaModify = [];
StageEnd(1).StageEndDef = {'SS',1e-4};
StageEnd(2).DividingCell = {};
StageEnd(2).SimParaModify = [];
StageEnd(2).StageEndDef = {'Duration',100000};

%% simulation
tspan = [0,100000];
N = ceil(75*ones(1,3)/dx);
L = N*dx;
S1 = 12;
S2 = 0;
noiseValue = 0.25;
ICDirName = 'N400r5dx0d5_1';
CellParaName = '120280';
DefaultSimParameters.noise = noiseValue;
DefaultSimParameters.sigma = [0.9,0.4; 0.4,0.4;];
SimName = sprintf('Symdt%sCT1_%s_2',dot2d(dt),CellParaName);
load(sprintf('%s/CellParameters_%s_1.mat',ICDirName,CellParaName),'CellParameters')

SimConfig = {'StabilizationPara',[S1,S2],'SavePhi',1000,'GPU','EnableEggshell','EnforceStage',1,'EnableCrop'};
if tspan(1)==0
  SimConfig = cat(2,SimConfig,{'LoadDir',ICDirName});
end
SimReport = mainBDF2VC_synNotch2(SimName,L,N,AxisLength,dt,DefaultSimParameters,CellParameters,StageEnd,tspan,SimConfig{:});