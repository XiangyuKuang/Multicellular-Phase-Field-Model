function SimReport = mainBDF2VC_synNotch2(SimName,L,N,EggshellPara,h,SimParameters,CellParameters,StageEnd,TimeSpan,varargin)
% SimReport = {N,h,StageEndTime,SimEndInfo,escapeT,trajectory,velocity,CellIdx};
%   StageEndTime format: 2x(N-1), N = StageEndNum
%             |StageEnd(2)|,...,|StageEnd(N)|
%            T|           |     |           |
%     TimeCost|           |     |           |

%% input config
dr = L./N;
dv = prod(dr);
[xx,yy,zz] = meshgrid(0:dr(1):L(1)-dr(1),0:dr(2):L(2)-dr(2),0:dr(3):L(3)-dr(3));% periodic, phi(0)=phi(L)
GPU = false;
SavePhi = false;
LoadDir = SimName;
EnableEggshell = false;
EnableCrop = false;
s = [0,0];
EnforceStageIdx = [];
for StrIdx = find(cellfun(@ischar,varargin))
  switch varargin{StrIdx}
    case 'SavePhi'
      SavePhi = true;
      SaveStep = varargin{StrIdx+1};
    case 'Uncompressed'
      EggshellPara(4) = -1;
    case 'GPU'
      GPU = true;
    case 'LoadDir'
      LoadDir = varargin{StrIdx+1};
    case 'StabilizationPara'
      s = varargin{StrIdx+1};
    case 'EnforceStage'
      EnforceStageIdx = varargin{StrIdx+1};
    case 'EnableEggshell'
      EnableEggshell = true;
    case 'EnableCrop'
      EnableCrop = true;
  end
end
if exist(SimName,'dir') ~= 7% make dir named by SimName
  mkdir(SimName)
end
% if GPU && SavePhi% init parpool for asynchronously saveing(ParaSave)
%   p = gcp;
%   p.IdleTimeout = inf;
%   f = parfeval(p,@rand,1);
% end
%% init log
[fileID, msg] = fopen(sprintf('%s/%s_log',SimName,SimName),'a');
if fileID < 0
  error('Failed to open file because: %s', msg);
end
fprintf(fileID,'+%s+\n',repmat('-',1,28));
%% load or generate initial state
% load initial state
[PhiStack,CellIdx,StageBreak,StageEndIdx] = LoadInitState(TimeSpan(1),LoadDir,StageEnd,fileID);
if ~isempty(EnforceStageIdx)
  StageEndIdx = EnforceStageIdx:length(StageEnd)-1;
  if isempty(StageEndIdx)
    error('Enforce StageIdx larger than length(StageEnd)-1')
  end
end
save(sprintf('%s/%s_StageEnd.mat',SimName,SimName),'StageEnd')
CellNum = length(CellIdx);
% load or create analsis var
trajectory = InitAnalysisVar(SimName,fileID,'trajectory',1);
velocity = InitAnalysisVar(SimName,fileID,'velocity',1);
% generate eggshell
if EnableEggshell
  PhiShell = zeros(size(xx));
  PhiShell(1,:,:) = 1;
  PhiShell(:,1,:) = 1;
  PhiShell(:,:,1) = 1;
else
  PhiShell = 0;
end
%% init var
% wavenum
[kx,ky,kz] = meshgrid(fftshift(-floor(N(1)/2):ceil(N(1)/2)-1),fftshift(-floor(N(2)/2):ceil(N(2)/2)-1),fftshift(-floor(N(3)/2):ceil(N(3)/2)-1));
kx = 2*pi/L(1)*1i*kx;
ky = 2*pi/L(2)*1i*ky;
kz = 2*pi/L(3)*1i*kz;
% transfer var to gpu if GPU is true
if GPU
  PhiStack = cellfun(@gpuArray,PhiStack,'UniformOutput',0);
  PhiShell = gpuArray(PhiShell);
  kx = gpuArray(kx);
  ky = gpuArray(ky);
  kz = gpuArray(kz);
  xx = gpuArray(xx);
  yy = gpuArray(yy);
  zz = gpuArray(zz);
end
ksq = kx.^2+ky.^2+kz.^2;
t = TimeSpan(1);% time in simulation
TimeCost = 0;% runing time
SimBreak = false;% terminate simulation if true, updated by IsEnd
AnalysisStep = 10;% do analysis every AnalStep
StageEndTime = zeros(2,StageEndIdx(end));% record t when StageEnd

%% simulation
fprintf(fileID,'simulation %s in progress...\n%0.1f Init\n',SimName,t);
for StageIdx = StageEndIdx
  if SimBreak% terminate simulation if true
    continue
  end
  % init velocity&trajectory
  r = cell(CellNum+1,1);
  v = cell(CellNum+1,1);
  for n = 1:CellNum
    r{n} = gather([sum(PhiStack{n}.*xx,'all'),sum(PhiStack{n}.*yy,'all'),sum(PhiStack{n}.*zz,'all')]/sum(PhiStack{n},'all'));
  end
  r{end} = t;
  % update parameters
  tau = SimParameters.tau;
  M = SimParameters.M;
  gamma = SimParameters.gamma*ones(CellNum,1);
  c = SimParameters.c;
  g = SimParameters.g;
  sigma = SimParameters.sigma;
  noiserdt = SimParameters.noise*sqrt(h);
  dvrV = dv./CellParameters.Volume(CellIdx);
  CellType = CellParameters.CellType(CellIdx);
  Type = unique(CellType)';
  TypeNum = length(Type);
  % init sim var
  [gammaValue,~,gammaIdx] = unique(gamma);
  EggShellRep = g*10*PhiShell;
  StageBreak = false;
  NextStageEnd = StageEnd(StageIdx+1).StageEndDef;
  [DefSwitch, DefIdx]= ismember({'Duration';'SS';'QSS'},NextStageEnd(:,1));
  t0 = t;% time at initiation of stage. t0=TimeSpan(1) for the 1st stage
  countQSS = 0;% count QSS for StageEndDef
  % save ICS
%   if SavePhi
%     StageIdxSave = [StageIdx;StageIdx+1];
%     if GPU
%       fetchNext(f);
%       f = parfeval(p,@ParaSave,0,'%s/%s-IC_%s.mat',{SimName,SimName,dot2d(t)},cellfun(@gather,PhiStack,'UniformOutput',0),CellIdx,SimParameters,StageIdxSave,EnableCrop);
%     else
%       ParaSave('%s/%s-IC_%s.mat',{SimName,SimName,dot2d(t)},PhiStack,CellIdx,SimParameters,StageIdxSave,EnableCrop);
%     end
%   end
  % computation
  PhiXStack = cell(CellNum,1);
  PhiYStack = cell(CellNum,1);
  PhiZStack = cell(CellNum,1);
  FOld = cell(CellNum,1);
  PhiStackOld = PhiStack;
  PhiXSum = cell(TypeNum,1);
  PhiYSum = cell(TypeNum,1);
  PhiZSum = cell(TypeNum,1);
  A = arrayfun(@(x)(tau/h+s(1)+2*x*c-x*ksq).^-1,gammaValue,'UniformOutput',0);
  B = tau/h+s(1);
  C = -0.5*tau/h-s(1);
  D = 2*gamma*c;
  % compute phi^1_i with first order scheme
  tic
  fPhi = fftn(PhiStack{1});
  PhiXStack{1} = real(ifftn(fPhi.*kx));
  PhiYStack{1} = real(ifftn(fPhi.*ky));
  PhiZStack{1} = real(ifftn(fPhi.*kz));
  Rep = PhiStack{1}.^2;
  for n = 2:CellNum
    fPhi = fftn(PhiStack{n});
    PhiXStack{n} = real(ifftn(fPhi.*kx));
    PhiYStack{n} = real(ifftn(fPhi.*ky));
    PhiZStack{n} = real(ifftn(fPhi.*kz));
    Rep = Rep + PhiStack{n}.^2;
  end
  Rep = EggShellRep + g*Rep;
  MRV = M*(1-cellfun(@(x)sum(x,'all'),PhiStack).*dvrV);
  for m = Type
    CellTypeIdx = find(CellType == m)';
    PhiXSum{m} = PhiXStack{CellTypeIdx(1)};
    PhiYSum{m} = PhiYStack{CellTypeIdx(1)};
    PhiZSum{m} = PhiZStack{CellTypeIdx(1)};
    for n = CellTypeIdx(2:end)
      PhiXSum{m} = PhiXSum{m} + PhiXStack{n};
      PhiYSum{m} = PhiYSum{m} + PhiYStack{n};
      PhiZSum{m} = PhiZSum{m} + PhiZStack{n};
    end
  end
  for n = 1:CellNum
    ut = noiserdt*randn(3,1);% noise
    AtrX = sigma(CellType(n),CellType(n))*(PhiXSum{CellType(n)}-PhiXStack{n})+ut(1);
    AtrY = sigma(CellType(n),CellType(n))*(PhiYSum{CellType(n)}-PhiYStack{n})+ut(2);
    AtrZ = sigma(CellType(n),CellType(n))*(PhiZSum{CellType(n)}-PhiZStack{n})+ut(3);
    for m = Type(Type~=CellType(n))
      AtrX = AtrX + sigma(CellType(n),m)*PhiXSum{m};
      AtrY = AtrY + sigma(CellType(n),m)*PhiYSum{m};
      AtrZ = AtrZ + sigma(CellType(n),m)*PhiZSum{m};
    end
    PhiSq = PhiStack{n}.^2;
    FOld{n} = MRV(n)*sqrt(PhiXStack{n}.^2+PhiYStack{n}.^2+PhiZStack{n}.^2)...%volume constrain
      -(Rep-g*PhiSq).*PhiStack{n}...%repulsion
      -AtrX.*PhiXStack{n}-AtrY.*PhiYStack{n}-AtrZ.*PhiZStack{n}...%attraction
      -D(n)*PhiSq.*(2*PhiStack{n}-3);%double well
    fPhi = A{gammaIdx(n)}.*fftn(B.*PhiStack{n}+FOld{n});
    PhiStack{n} = real(ifftn(fPhi));
    PhiXStack{n} = real(ifftn(fPhi.*kx));
    PhiYStack{n} = real(ifftn(fPhi.*ky));
    PhiZStack{n} = real(ifftn(fPhi.*kz));
    MRV(n) = M*(1-real(fPhi(1))*dvrV(n));
  end
  TimeCost = TimeCost + toc;
  tAnalysis = t0;
  t = t+h;
  IntervalLength = round((tAnalysis+AnalysisStep-t)/h);% IntervalLength will be updated by IsEnd
  % compute phi^n_i with BDF2 scheme for n>1
  A = arrayfun(@(x)(1.5*tau/h+2*x*c+s(1)-x*ksq).^-1,gammaValue,'UniformOutput',0);
  B = 2*tau/h+2*s(1);
  while ~(StageBreak||SimBreak)
    tic
    for iter = 1:IntervalLength
      Rep = PhiStack{1}.^2;
      for n = 2:CellNum
        Rep = Rep + PhiStack{n}.^2;
      end
      Rep = EggShellRep + g*Rep;
      for m = Type
        CellTypeIdx = find(CellType == m)';
        PhiXSum{m} = PhiXStack{CellTypeIdx(1)};
        PhiYSum{m} = PhiYStack{CellTypeIdx(1)};
        PhiZSum{m} = PhiZStack{CellTypeIdx(1)};
        for n = CellTypeIdx(2:end)
          PhiXSum{m} = PhiXSum{m} + PhiXStack{n};
          PhiYSum{m} = PhiYSum{m} + PhiYStack{n};
          PhiZSum{m} = PhiZSum{m} + PhiZStack{n};
        end
      end
      for n = 1:CellNum
        ut = noiserdt*randn(3,1);% noise
        AtrX = sigma(CellType(n),CellType(n))*(PhiXSum{CellType(n)}-PhiXStack{n})+ut(1);
        AtrY = sigma(CellType(n),CellType(n))*(PhiYSum{CellType(n)}-PhiYStack{n})+ut(2);
        AtrZ = sigma(CellType(n),CellType(n))*(PhiZSum{CellType(n)}-PhiZStack{n})+ut(3);
        for m = Type(Type~=CellType(n))
          AtrX = AtrX + sigma(CellType(n),m)*PhiXSum{m};
          AtrY = AtrY + sigma(CellType(n),m)*PhiYSum{m};
          AtrZ = AtrZ + sigma(CellType(n),m)*PhiZSum{m};
        end
        PhiSq = PhiStack{n}.^2;
        F =  MRV(n)*sqrt(PhiXStack{n}.^2+PhiYStack{n}.^2+PhiZStack{n}.^2)...%volume constraint
          -(Rep-g*PhiSq).*PhiStack{n}...%repulsion
          -AtrX.*PhiXStack{n}-AtrY.*PhiYStack{n}-AtrZ.*PhiZStack{n}...%attraction
          -D(n)*PhiSq.*(2*PhiStack{n}-3);%double well
        fPhi = A{gammaIdx(n)}.*fftn(B.*PhiStack{n}+C.*PhiStackOld{n}+2*F-FOld{n});
        PhiStackOld{n} = PhiStack{n};
        PhiStack{n} = real(ifftn(fPhi));
        FOld{n} = F;
        PhiXStack{n} = real(ifftn(fPhi.*kx));
        PhiYStack{n} = real(ifftn(fPhi.*ky));
        PhiZStack{n} = real(ifftn(fPhi.*kz));
        MRV(n) = M*(1-real(fPhi(1))*dvrV(n));
      end
    end
    TimeCost = TimeCost + toc;
    t = t+IntervalLength*h;
    tAnalysis = tAnalysis+AnalysisStep;
    % analysis
    for n = 1:CellNum
      rc = gather([sum(PhiStack{n}.*xx,'all'),sum(PhiStack{n}.*yy,'all'),sum(PhiStack{n}.*zz,'all')]/(1-MRV(n)/M)*dvrV(n));
      r{n} = cat(1,r{n},rc);
      v{n} = cat(1,v{n},(r{n}(end,:)-r{n}(end-1,:))/h/IntervalLength);
    end
    r{end} = cat(1,r{end},t);
    v{end} = cat(1,v{end},t);
    % StageEndCheck
    [StageBreak, SimBreak, IntervalLength, StageEndInfo, SimEndInfo, countQSS] = IsEnd(NextStageEnd, DefSwitch, DefIdx, PhiStack, v(1:end-1), countQSS, t, t0, tAnalysis, TimeSpan(2), AnalysisStep, h, fileID);
    if DefSwitch(3)
      if strcmp('QSS',StageEndInfo)% return to last analysis step(t-IntervalLength*h) if QSS
        PhiStack = PhiOld;
        t = tOld;
      end
      PhiOld = PhiStack;
      tOld = t;
    end
    % save
    if SavePhi && (mod(length(v{end}),SaveStep) == 0||StageBreak||SimBreak)
      if StageBreak
        StageIdxSave = StageIdx+1;
      else
        StageIdxSave = [StageIdx;StageIdx+1];
      end
%       if GPU
%         % asynchronously saveing
%         fetchNext(f);
%         f = parfeval(p,@ParaSave,0,'%s/%s_%s.mat',{SimName,SimName,dot2d(t)},cellfun(@gather,PhiStack,'UniformOutput',0),CellIdx,SimParameters,StageIdxSave,EnableCrop);
%       else
      ParaSave('%s/%s_%s.mat',{SimName,SimName,dot2d(t)},cellfun(@gather,PhiStack,'UniformOutput',0),CellIdx,SimParameters,StageIdxSave,EnableCrop);
%       end
    end
  end
  % update analysis result
  StageEndTime(1,StageIdx) = t;
  StageEndTime(2,StageIdx) = TimeCost;
  trajectory = cat(2,trajectory,{r;CellIdx});
  velocity = cat(2,velocity,{v;CellIdx});
end

%% close
SimReport = {N,h,StageEndTime,SimEndInfo,TimeCost,trajectory,velocity,CellIdx};
if ~StageBreak
  ParaSave('%s/%s_%s.mat',{SimName,SimName,dot2d(t)},cellfun(@gather,PhiStack,'UniformOutput',0),CellIdx,SimParameters,StageIdx,EnableCrop)
end
save(sprintf('%s/%s_trajectory.mat',SimName,SimName),'trajectory');
save(sprintf('%s/%s_velocity.mat',SimName,SimName),'velocity');
fprintf(fileID,'simulation elapsed time is %0.3fs\n',TimeCost);
fclose(fileID);
% if SavePhi && GPU
%   fetchNext(f);
% end
end

function [PhiStack,CellIdx,StageBreak,StageEndIdx] = LoadInitState(t0,DirName,StageEnd,fileID)
[tSave, MatIdx, FileList] = GetTime(DirName);
LoadMatIdx = MatIdx(abs(tSave-t0)<1e-4);% floating-point problem
% load phi
if isempty(LoadMatIdx)
  error('.mat at %0.1f do not exist',t0)
end
FileName = sprintf('%s/%s',DirName,FileList{LoadMatIdx});
load(FileName,'phi_save','CellIdx','StageIdx')
fprintf(fileID,'load %s\n',FileName);
StageBreak = length(StageIdx) == 1;
PhiStack = phi_save;
StageEndIdx = StageIdx(1):length(StageEnd)-1;
end

function AnalysisResult = InitAnalysisVar(SimName,fileID,AnalysisName,LoadSwitch)
if exist(sprintf('%s/%s_%s.mat',SimName,SimName,AnalysisName),'file')
  if LoadSwitch
    load(sprintf('%s/%s_%s.mat',SimName,SimName,AnalysisName),AnalysisName);
    AnalysisResult = eval(AnalysisName);
    fprintf(fileID,'load %s_%s.mat\n',SimName,AnalysisName);
  else
    error('%s file already exist !',AnalysisName)
  end
else
  AnalysisResult = {};
  fprintf(fileID,'create %s_%s.mat\n',SimName,AnalysisName);
end
end

function ParaSave(StrF,StrStack,phi_save,CellIdx,SimParameters,StageIdx,EnableCrop)
unsaved = 1;
n = 0;
MatName = sprintf(StrF,StrStack{:});
while unsaved
  try
    if EnableCrop
      [phi_cropped, box] = crop_phi(phi_save,1e-2);
      save(MatName,'phi_cropped','CellIdx','SimParameters','box','StageIdx','-v7.3')
    else
      save(MatName,'phi_save','CellIdx','SimParameters','StageIdx','-v7.3')
    end
    unsaved = 0;
  catch
    unsaved = 1;
    n = n+1;
  end
  if n>10
    warning('%s unsaved\n',MatName)
    break
  end
end
end

function [StageBreak, SimBreak, IntervalLength, StageEndInfo, SimEndInfo, countQSS] = IsEnd(StageEndDef, DefSwitch, DefIdx, PhiStack, v, countQSS, t, t0, tAnalysis, SimEndTime, AnalStep, h, fileID)
StageBreak = 0;
SimBreak = 0;
StageEndInfo = '';
SimEndInfo = '';
T = SimEndTime-t0;% period of stage
Telasped = t-t0;
% SimEnd
if any(cellfun(@(x)any(isnan(x),'all'),PhiStack))% unstable
  SimEndInfo = 'Unstable';
  SimBreak = 1;
end
if any(cellfun(@(x)all(x<0.5,'all'),PhiStack))%CellLoss
  SimEndInfo = 'CellLoss';
  SimBreak = 1;
end
if t >= SimEndTime%SimTimeout
  SimEndInfo = 'Timeout';
  SimBreak = 1;
end
% StageEnd
vRMS = rms(cell2mat(cellfun(@(x)vecnorm(x(max(end-2,1):end,:)',2),v,'UniformOutput',0)));
if DefSwitch(1)%Duration
  T = min(T,StageEndDef{DefIdx(1),2});
  if Telasped>=StageEndDef{DefIdx(1),2}
    StageEndInfo = 'Duration';
    StageBreak = 1;
  end
end
if DefSwitch(2)%SS
  if vRMS(end) <= StageEndDef{DefIdx(2),2}
    StageEndInfo = 'SS';
    StageBreak = 1;
  end
end
if DefSwitch(3)%QSS
  countQSS = countQSS+any(islocalmin(vRMS));
  if countQSS==StageEndDef{DefIdx(3),2}
    StageEndInfo = 'QSS';
    StageBreak = 1;
    t = t-AnalStep;
  end
end
IntervalLength = min(round((tAnalysis+AnalStep-t)/h),ceil((T-Telasped)/h));
if StageBreak
  fprintf(fileID,'%0.1f reach %s\n',t,StageEndInfo);
end
if SimBreak
  fprintf(fileID,'%0.1f simulation reach %s\n',t,SimEndInfo);
  return
end
end
