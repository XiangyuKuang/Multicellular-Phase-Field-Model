[r,c] = size(answer1);
x = 1:c;
y = 1:r;
[xx,yy] = meshgrid(x,y);
yy = flipud(yy);
hold on;
scatter(xx(:),yy(:),answer1(:),"black",'filled');

% 坐标轴美化
axis equal  
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
x_labels = {'','ABa','ABp','P1',''};
x_colors = {[0 0 0], [0 0 1], [0 0 1], [1 0 0], [0 0 0]}; 
for i = 1:length(x_labels)
    text(i-1, -0.2, x_labels{i}, 'Color', x_colors{i}, 'FontSize', 18, 'HorizontalAlignment', 'center', 'FontName', 'Arial');
end

% 添加 Y 轴标签，向下平移一格，并稍微向右移动
y_labels = {'','P1','ABp','ABa',''};
y_colors = {[0 0 0], [1 0 0], [0 0 1], [0 0 1], [0 0 0]}; 
for i = 1:length(y_labels)
    text(-0.1, i-1, y_labels{i}, 'Color', y_colors{i}, 'FontSize', 18, 'HorizontalAlignment', 'right', 'FontName', 'Arial');
end

a = 0:c+1;
b = 0:c+1;
b = -a + c + 1;
plot(a,b,'k--','LineWidth',2.5);
