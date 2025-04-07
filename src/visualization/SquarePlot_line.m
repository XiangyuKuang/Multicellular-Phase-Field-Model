function SquarePlot_line(data,CellName,SigmaColor,ax)
if nargin<=3
  ax = gca;
end
CellNum = length(CellName);
f = [1:4;5:8];
for n = 1:CellNum
  for m = n:CellNum 
    v = [[n-[0;1;1;0],m-[0;0;1;1]];[m-[0;1;1;0],n-[0;0;1;1]]];
    if n == m
      c = ones(1,3);
    else
      c = SigmaColor(SigmaColor(:,1) == data(n,m),2:4);
    end
    patch(ax,'Faces',f,'Vertices',v,'FaceColor',c,'EdgeColor',ones(1,3),'LineWidth',20/CellNum)
  end
end
hold(ax,'on')
plot(ax,[0,CellNum],[0,CellNum],'k--','LineWidth',1)
dx = -0.05*max(cellfun(@length,CellName));
if CellNum>8
  FontSize = 144/CellNum;
else
  FontSize = 15;
end
axis(ax,'equal',[0,CellNum,0,CellNum])
set(ax,'FontSize',FontSize,'Fontname','arial','ydir','reverse',...
  'xtick',(0.5:1:0.5+CellNum-1)+dx,'xticklabel',CellName,'xticklabelrotation',60,...
  'ytick',0.5:1:0.5+CellNum-1,'yticklabel',CellName)