function AutoImg(varargin)
%Example: AutoImg('name',[159999,172999,183299],'ImgType',3,'TitleFormat','%d','save')
%% Input
ImgSave = false;
Specified = false;
CloseFig = false;
EnableLoadIC = false;
FormatSpec = [];
ContactComp = [];
CloseLegendAxis = [false,false];
% AxesConfig = {'unit','normalized','position',[0,0,0.4,0.6]};
AxesConfig = {'unit','pixels','position',[0,0,2048/3,640]};
th = 0.34;
for StrIdx = find(cellfun(@ischar,varargin))
  switch varargin{StrIdx}
    case 'dx'
      dx = varargin{StrIdx+1};
    case 'ImgType'
      ImgType = varargin{StrIdx+1};
    case 'save'
      ImgSave = true;
    case 't'
      select_t = varargin{StrIdx+1};
      Specified = true;
    case 'TitleFormat'
      FormatSpec = varargin{StrIdx+1};
      if ~contains(FormatSpec,'%')
        error('TitleFormat must be a FormatSpec')
      end
    case 'threshold'
      th = varargin{StrIdx+1};
    case 'close'
      CloseLegendAxis = contains({'legend','axis'},varargin{StrIdx+1});
      CloseFig = contains({'fig'},varargin{StrIdx+1});
    case 'ContactComp'
      ContactComp = varargin{StrIdx+1};
    case 'CellPara'
      CellParameters = varargin{StrIdx+1};
    case 'AxesConfig'
      AxesConfig = varargin{StrIdx+1};
    case 'IC'
      EnableLoadIC = true;
  end
end
if ~exist('CellParameters','var')
  CellParameters = GetCellPara;
end
if EnableLoadIC
  [matT,matIdx,FileList] = GetTime('IC');
else
  [matT,matIdx,FileList] = GetTime;
end
if Specified
  idx = [];
  for n = 1:length(select_t)
    m = find(select_t(n) == matT);
    if isempty(m)
      warning('no t=%0.1f in current dir %s',select_t(n),pwd)
    end
    idx = cat(1,idx,m);
  end
  MatLoadIdx = matIdx(idx);
  t = matT(idx);
else
  MatLoadIdx = matIdx';
  t = matT;
end
if isempty(t)
  error('no phase-field data in current dir %s',pwd)
end
if ImgSave
  AxesConfig = cat(2,AxesConfig,{'toolbar','none'});
end
%% image
FigStack = cell(length(ImgType),2);
TitleText = [];
for n = 1:length(MatLoadIdx)
  MatName = FileList(MatLoadIdx(n));
  load(MatName,'phi_save','CellIdx','SimParameters');
  if ~isempty(FormatSpec)
    TitleText = sprintf(FormatSpec, t(n));
  end
  for m = 1:length(ImgType)
    FigStack{m,1} = figure;
    switch ImgType(m)
      case 1
        set(gcf,AxesConfig{:});
        if CloseLegendAxis(1)
          MultiphaseDisp(phi_save,CellIdx,dx,CellParameters,'title',TitleText,'threshold',th,'LegendOff','AxisFontSize',25.5,'revydir');
        else
          MultiphaseDisp(phi_save,CellIdx,dx,CellParameters,'title',TitleText,'threshold',th,'revydir');
        end
        if CloseLegendAxis(2)
          axis off
        end
      case 2
        set(gcf,'unit','pixels','position',[0,0,2816/3,960]);
        MultiphaseDisp(phi_save,CellIdx,dx,CellParameters,'title',TitleText,'threshold',th,'sigma',SimParameters.sigma,'ContactComp',ContactComp,'revydir');
      case 3
%         set(gcf,'unit','normalized','position',[0,0,0.55,0.9]);
        set(gcf,'unit','pixels','position',[0,0,2816/3,960]);
        MultiphaseDisp(phi_save,CellIdx,dx,CellParameters,'title',TitleText,'threshold',th,'sigma',SimParameters.sigma,'ContactComp',ContactComp,'TripleView','revydir');
      case 4
        set(gcf,AxesConfig{:});
        if CloseLegendAxis(1)
          Celegans3LayerDisp(phi_save,CellIdx,dx,CellParameters,'title',TitleText,'threshold',[0.2,0.34,0.5],'LegendOff','AxisFontSize',25.5);
        else
          Celegans3LayerDisp(phi_save,CellIdx,dx,CellParameters,'title',TitleText,'threshold',[0.2,0.34,0.5]);
        end
        if CloseLegendAxis(2)
          axis off
        end
    end
    FigStack{m,2} = ImgType(m);
  end
  if ImgSave
    for m = 1:length(ImgType)
      lower = find(MatName{1} == '.',1,'last')-1;
      if isempty(ContactComp)
        str = sprintf('(%d',FigStack{m,2});
      else
        str = '_comp';
      end
      print(FigStack{m,1},sprintf('%dcell_%s%s.png',length(phi_save),MatName{1}(1:lower),str),'-dpng','-r800')
    end
  end
  if CloseFig
    close all
  end
end