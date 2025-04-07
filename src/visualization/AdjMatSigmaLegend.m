function AdjMatSigmaLegend(DispContact,CompContact,sigma,CellName,axStack,CellColor)
% 4 inputs get contact map: AdjMatSigmaLegend(DispContact,CompContact,sigma,CellName)
% 6 inputs get contact map and legend: AdjMatSigmaLegend(DispContact,CompContact,sigma,CellName,ax,CellColor)
if nargin<=4
  axStack{1} = gca;
end
CellNum = length(CellName);
% sigma
if sigma == 0
  sigma = ones(CellNum);
  SigmaValue = 1;
%   c = 0.85;
  c = 1;% exp cell plot
else
  SigmaValue = unique(sigma(:));
  if length(SigmaValue)>2
    c = SigmaValue*-0.8+0.85;
%     c = [1,0.8*2.^-(0:length(SigmaValue)-2)]';
  else
    c = [0.85;0.45];%11/15/23 12cell plot
%     c = [1;0.4];
  end
end
SquarePlot_line(sigma,CellName,[SigmaValue,c*ones(1,3)],axStack{1})
% adjacency matrix 
MarkerSetting = {1, 'k.', 200/CellNum;% conserve
                 2, 'ko', 200/CellNum/4};%non-conserve
% 180/..->200/.. 11/15/23
SquarePlotMarker(DispContact,MarkerSetting,axStack{1})
% highlight broken conserved contact
if ~isempty(CompContact)
  if all(size(DispContact) == size(CompContact))
    if length(unique(CompContact))>2
      SimContact = DispContact;
      ExpContact = CompContact;
    else
      ExpContact = DispContact;
      SimContact = CompContact;
    end
    % conserve 2
    % non-conserve 1
    HighlightSquare(((SimContact-(ExpContact==2))>0&ExpContact~=1)|(SimContact-(ExpContact==2))<0,[1,0,0],axStack{1});
  else
    warning('Contact map size is not consistant')
  end
end
% legend
if nargin > 5
  f = 1:4;
  w = 3*6/max([3,CellNum]);
  for n = 1:CellNum
    vert = [[0;1;1;0],[0;0;1;1]-n];
    patch(axStack{2},'Faces',f,'Vertices',vert,'FaceColor',CellColor(n,:),'EdgeColor','none')
%     patch(ax{2},'Faces',f,'Vertices',vert,'FaceColor',CellColor(n,:),'EdgeColor',[1,1,1],'LineWidth',w)
  end
  hold(axStack{2},'on')
  plot([zeros(CellNum+1,1),ones(CellNum+1,1)]',-[0:CellNum;0:CellNum],'LineWidth',w,'color',[1,1,0.99])
  axis(axStack{2},[0,1,-CellNum,0],'off')
%   set(ax{2},'xticklabels',[],'yticklabels',[])
end