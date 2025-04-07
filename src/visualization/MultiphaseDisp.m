function PatchStack = MultiphaseDisp(varargin)
PhiStack = varargin{1};
CellIdx = varargin{2};
dr = varargin{3};
CellParameters = varargin{4};

%% image setting
TitleSettings = {[],'FontSize',20,'Fontname','arial','fontweight','bold'};%,'fontweight','bold'};
sigma = [];
ContactComp = [];
FaceAlpha = 0.7;
th = 0.34;
[CellIdx,I] = sort(CellIdx);
PhiStack = PhiStack(I);
color = CellParameters.Color(CellIdx,:);
TripleView = 0;
EnableLegend = 1;
space=0;
EnableBox = 0;
EnablePadding0 = 0;
AxisLimit = [];
ViewVec = [0,0,-1];
AxisFontSize = 20;
LegendFontSize = 17.5;
theta = zeros(1,3);
ydirStr = {};
box = [];
for StrIdx = find(cellfun(@ischar,varargin))
	switch varargin{StrIdx}
		case 'title'
			TitleSettings{1} = varargin{StrIdx+1};
		case 'threshold'
			th = varargin{StrIdx+1};
		case 'sigma'
			sigma = varargin{StrIdx+1};
      sigma = sigma(I,I);
		case 'FaceAlpha'
			FaceAlpha = varargin{StrIdx+1};% scalar, change FaceAlpha of all cell
		case 'TripleView'
			TripleView = true;
    case 'ContactComp'
      ContactComp = varargin{StrIdx+1};
    case 'LegendOff'
      EnableLegend=0;
    case 'LegendSpace'
      space = varargin{StrIdx+1};
    case 'LegendFontSize'
      LegendFontSize = varargin{StrIdx+1};
    case 'color'
      % order of color is the same as order of cells in phi_save
      color = varargin{StrIdx+1}(I,:);
    case 'AxisLimit'
      AxisLimit = varargin{StrIdx+1};
    case 'view'
      ViewVec = varargin{StrIdx+1};
    case 'AddBox'
      EnableBox = 1;
    case 'AxisFontSize'
      AxisFontSize = varargin{StrIdx+1};
      TitleSettings{3} = AxisFontSize;
    case 'GridRot'
      theta = varargin{StrIdx+1};%L3 [-pi/8,-pi/8,pi/3]; PS []
    case 'padding0'
      EnablePadding0 = 1;
    case 'revydir'
      ydirStr = {'ydir','reverse'};
    case 'box'
      box = varargin{StrIdx+1};

	end
end

%% config cell name, color
if any(contains(fieldnames(CellParameters),'FaceAlpha'))
  FaceAlphaList = CellParameters.FaceAlpha(CellIdx);
else
  FaceAlphaList = FaceAlpha*ones(length(CellIdx),1);
end
% if isempty(color)
%   color = CellParameters.Color(CellIdx,:);
% end

if EnableLegend
  CellName = cellfun(@(x)sprintf(' %-*s',space,x),CellParameters.Name(CellIdx),'UniformOutput',false);%
  LegendSettings = {CellName,'FontSize',LegendFontSize,'FontName','arial','box','off'};
end

%% iso-surface
if isempty(box)
  N = size(PhiStack{1});
  L = [zeros(1,3),N([2,1,3]).*dr];
  if EnablePadding0
    PhiStack = PaddingPhi(PhiStack);
    N = size(PhiStack{1});
    L = [zeros(1,3),N([2,1,3])*dr];
  end
  L = [L(1:3)-L(4:6)/2,L(4:6)-L(4:6)/2];
  if isempty(AxisLimit)
    AxisLimit = L;
  end
  [X,Y,Z] = meshgrid(L(1):dr:L(4)-dr,L(2):dr:L(5)-dr,L(3):dr:L(6)-dr);
  if any(theta)
    [X,Y,Z] = rotGrid(X,Y,Z,theta);
  end
  fv = cellfun(@(x)isosurface(X,Y,Z,x,th),PhiStack,'UniformOutput',0);
else
  N = max(box(:,[2,4,6]));%[150,150,150];% to do : add input
  L = [-N([2,1,3]).*dr/2,N([2,1,3]).*dr-N([2,1,3]).*dr/2];
  if isempty(AxisLimit)
    AxisLimit = L;
  end
  [X,Y,Z] = meshgrid(L(1):dr:L(4)-dr,L(2):dr:L(5)-dr,L(3):dr:L(6)-dr);
  fv=cellfun(@(x)isosurface(x,th),PhiStack,'UniformOutput',0);
  for  n = 1:length(PhiStack)
    fv{n}.vertices = fv{n}.vertices+box(n,[3,1,5])+[-2,-2,-2];
    fv{n}.vertices = fv{n}.vertices*dr+L(1:3);
  end
end
%% imaging
if TripleView
	y0 = 0.2;
	x0 = 0.1;
	xw = [15,30,30]*0.75/45;
	yw = [20,20,15]*0.75/45;
	ViewDir = [
    -1,0,0;% DV-LR
    0,0,-1;% AP-DV
    0,-1,0];% AP-LR
  CamRollAngle = [90,0,0];
	ObjStack = cell(3,1);
  for ViewIdx = 1:3
    ObjStack{ViewIdx} = axes('Position',[x0,y0,xw(ViewIdx),yw(ViewIdx)]);
    IsoSurfacePlot(fv,color,FaceAlphaList,20,ViewDir(ViewIdx,:),AxisLimit,X,Y,Z,PhiStack,box,ydirStr);
    camroll(CamRollAngle(ViewIdx))
    if ViewIdx<3
      if xw(ViewIdx) ~= xw(ViewIdx+1)
        x0 = x0+0.1+xw(ViewIdx);
      end
      if yw(ViewIdx) ~= yw(ViewIdx+1)
        y0 = y0+0.1+yw(ViewIdx);
      end
    end
  end
  set(ObjStack{1},'YAxisLocation','right')
	if ~isempty(sigma)
		% legend, adjacent matrix and sigma
    ax = {axes('Position',[0.1,y0,0.25,0.25]),axes('Position',[0.02,y0,0.05*0.25,0.25])};
%     ax = {axes('Position',[0.1,y0,0.25,0.25]),axes('Position',[0.01,y0,0.06*0.25,0.25])};
  else
    % only legend
		h = findobj(gca);
		legend(h(end:-1:3),LegendSettings{:},'position',[0,0.76,0.425,0]);
	end
	sgtitle(TitleSettings{:});
else
	width = 0.32;
	gap = 0.05;
	y0 = 0.5-width/2;
	if isempty(sigma)
    title(TitleSettings{:},'Position',[2,20.5,0]);
  else
    axes('Position',[width+gap*3.2,y0,1.5*width,width]);
  end
  PatchStack = IsoSurfacePlot(fv,color,FaceAlphaList,AxisFontSize,ViewVec,AxisLimit,X,Y,Z,PhiStack,box,ydirStr);
  if EnableBox
    AddBox(AxisLimit)
  end
  if EnableLegend
    if isempty(sigma)
      h = findobj(gca);
      legend(h(end:-1:3),LegendSettings{:},'Location','eastoutside','NumColumns',1+(length(CellIdx)>26));
      set(gca,'Position',[0.12 0.11 0.6528 0.8150])
    else
      % adjacent matrix and sigma
      sgtitle(TitleSettings{:});
      ax = {axes('Position',[1.7*gap,y0,width,width]),axes('Position',[0.01,y0,0.06*width,width])};
    end
  end
end
if ~isempty(sigma)&&EnableLegend
  AdjMatSigmaLegend(IsContact(PhiStack,0.35),ContactComp,sigma,CellName,ax,color) % 0.44 previous
end
end

%% subfunctions
function PatchStack = IsoSurfacePlot(fv,color,FaceAlphaList,FontSize,ViewVec,AxisLimit,X,Y,Z,PhiStack,box,ydirStr)
CellNum = length(fv);
PatchStack = cell(CellNum,1);
if isempty(box)
  for n = 1:CellNum
    PatchStack{n} = patch(fv{n},'FaceColor',color(n,:),'EdgeColor','none','FaceAlpha',FaceAlphaList(n));
    isonormals(X,Y,Z,PhiStack{n},PatchStack{n})
  end
else
  for n = 1:CellNum
    PatchStack{n} = patch(fv{n},'FaceColor',color(n,:),'EdgeColor','none','FaceAlpha',FaceAlphaList(n));
    isonormals( ...
      X(box(n,1):box(n,2),box(n,3):box(n,4),box(n,5):box(n,6)), ...
      Y(box(n,1):box(n,2),box(n,3):box(n,4),box(n,5):box(n,6)), ...
      Z(box(n,1):box(n,2),box(n,3):box(n,4),box(n,5):box(n,6)), ...
      PhiStack{n},PatchStack{n})
  end
end
axis('equal',AxisLimit([1,4,2,5,3,6]))
grid off
ylabel({'\itz\rm / D-V Axis (μm)'});
xlabel({'\itx\rm / A-P Axis (μm)'});
zlabel({'\ity\rm / L-R Axis (μm)'});
set(gca,'FontSize',FontSize,'Fontname','arial','tickdir','out','LineWidth',2,ydirStr{:});
% set(gca,'FontSize',FontSize,'Fontname','arial','tickdir','out','LineWidth',2);
view(ViewVec)
camlight(-20,40)
lighting gouraud
material dull
end

function AddBox(L)
hold on
[x,y,z] = meshgrid(L([1,4]),L([2,5]),L([3,6]));
v = [x(:),y(:),z(:)];
plotted = zeros(8,1);
for n = 1:8
  idx = find(sum(v(n,:)==v,2)==2);
  for m = find(plotted(idx)<3)'
    plot3(v([n,idx(m)],1),v([n,idx(m)],2),v([n,idx(m)],3),'k','linewidth',1);
  end
  plotted(n) = 3;
  plotted(idx) = plotted(idx)+1;
end
end

function PhiStacknew = PaddingPhi(PhiStack)
N = size(PhiStack{1});
PhiStacknew = cell(size(PhiStack));
for n = 1:length(PhiStack)
  PhiStacknew{n} = zeros(N+2);
  PhiStacknew{n}((1:N(1))+1,(1:N(2))+1,(1:N(3))+1) = PhiStack{n};
% PhiStacknew{n} = zeros(N+[0,0,1]);
% PhiStacknew{n}(:,:,1:N(3)) = PhiStack{n};
end
end
% function AddBackground(L)
% [x,y,z] = meshgrid(L([1,4]),L([2,5]),L([3,6]));
% v = [x(:),y(:),z(:)];
% f = [1,3,7,5;1,2,4,3;3,4,8,7];
% patch('Faces',f,'Vertices',v,'FaceColor','w','EdgeColor','none','FaceAlpha',0.5)
% end