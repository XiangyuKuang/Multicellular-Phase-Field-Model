clear;
load('E:\data\embryo07_cell_stage_time.mat');
data = xlsread(['E:\data\CD\CD07.csv'] );

for m=1:7
    time=embryo07_cell_stage_time(m);
    selected_data = data(:, [1, 7, 8, 9]);
    index = selected_data(:, 1) == time;

    % 将第9列和第10列乘以0.09，第11列乘以0.42
    selected_data(:, 2:4) = selected_data(:, 2:4) .* [0.42, 0.09, 0.09];
    z_data = selected_data(index, 2); % 第9列数据作为z轴坐标
    x_data = selected_data(index, 3); % 第10列数据作为x轴坐标
    y_data = selected_data(index, 4); % 第11列数据作为y轴坐标

    radius = 30;
    scatter3(x_data, y_data, z_data, 'filled', 'SizeData', radius,'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');

    daspect([1, 1, 1]);
    xlim([5 65])
    ylim([0 40])
    view(0,90);

    axis off;
    grid off;
    set(gcf, 'Color', 'k');
    set(gcf, 'InvertHardcopy', 'off'); % 阻止导出时背景颜色反转
    filename2 = sprintf('%d.png', m);
    print(gcf, filename2, '-dpng', '-r600');
end

