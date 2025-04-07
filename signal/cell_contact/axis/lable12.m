[r,c] = size(answer1);
x = 1:c;
y = 1:r;
[xx,yy] = meshgrid(x,y);
yy = flipud(yy);
%scatter(xx(:),yy(:),answer(:),"black");
hold on;
scatter(xx(:),yy(:),answer1(:),"black",'filled');
% 坐标轴美化
axis equal 
%hTitle = title('Contact Area');
%set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
set(gca, 'Box', 'on', ...                                              
         'TickDir', 'in', 'TickLength', [0 0], ...         
         'XMinorTick', 'off', 'YMinorTick', 'off', ...          
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...        
         'XTick', 0:1:c+1,...                                    
         'XLim', [0 c+1],...
         'YTick', 0:1:r+1,...
         'YLim', [0 r+1],...
         'XTickLabel', [], ... % 去掉 X 轴标签
         'YTickLabel', []);    % 去掉 Y 轴标签
set(gca,'FontSize',18,'Fontname', 'Arial');
% 添加 X 轴标签，向左平移一格，并稍微向上移动
x_labels = {'','ABala','ABalp','ABara','ABarp','ABpla','ABplp','ABpra','ABprp','MS','E','C','P3',''};
x_colors = {[0 0 0],[76 76 188]/255,[76 76 233]/255,[76 98 255]/255, [76 143 255]/255,[76 188 255]/255,[76 233 255]/255,[98 255 233]/255,[143 255 188]/255,[210 255 121]/255,[255 210 76]/255,[255 121 76]/255,[210 76 76]/255,[0 0 0]}; 
for i = 1:length(x_labels)
    text(i-1, -0.2, x_labels{i}, 'Color', x_colors{i}, 'FontSize', 18, 'HorizontalAlignment', 'center', 'FontName', 'Arial');
end

% 添加 Y 轴标签，向下平移一格，并稍微向右移动
y_labels = {'','P3','C','E','MS','ABprp','ABpra','ABplp','ABpla','ABarp','ABara','ABalp','ABala',''};
y_colors = {[0 0 0], [210 76 76]/255, [255 121 76]/255, [255 210 76]/255, [210 255 121]/255, [143 255 188]/255, [98 255 233]/255, [76 233 255]/255, [76 188 255]/255, [76 143 255]/255, [76 98 255]/255, [76 76 233]/255, [76 76 188]/255, [0 0 0]
}; 
for i = 1:length(y_labels)
    text(-0.1, i-1, y_labels{i}, 'Color', y_colors{i}, 'FontSize', 18, 'HorizontalAlignment', 'right', 'FontName', 'Arial');
end
a = 0:c+1
b = 0:c+1
b=-a+c+1;
plot(a,b,'k--','LineWidth',2.5)