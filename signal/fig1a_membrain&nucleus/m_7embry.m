clear;
close all;
load('E:\data\cellname.mat');
load('E:\data\color_cell_type.mat');
load('E:\data\embryo07_cell_stage_time.mat');
name_dictionary = readtable('E:\data\CShaper_1\name_dictionary.csv');

for m=1:7
    time=embryo07_cell_stage_time(m); 
    % 导入NIfTI文件
    nii = load_nii(['E:\data\CShaper_1\Sample07\SegCell\Sample07_',num2str(time,'%03d'),'_segCell.nii.gz']);
    niiData = nii.img;
    % 获取NIfTI数据中的非零值
    non_zero_values = unique(niiData(niiData ~= 0));

    % 创建网格
    [x, y, z] = ndgrid(1:size(niiData, 1), 1:size(niiData, 2), 1:size(niiData, 3));
    color = zeros(length(non_zero_values),3);

    C=string(Cellname{1,m}); %选择细胞名称
    Colors=colors{m}; %选择各细胞对应颜色

    for i = 1:length(non_zero_values)
        idx = find(name_dictionary{:,1} == non_zero_values(i));
        if ~isempty(idx)
            cell_name = name_dictionary{idx, 2};
            % 在 C 中查找该细胞名称的列索引 m
            m = find(strcmp(C, cell_name));
            if ~isempty(m)
                color(i, :) = Colors(m, :);
            end
        end
    end

    % 循环处理每个非零值对应的细胞空间
    figure;
    hold on;
    for a = 1:length(non_zero_values)
        space_value = non_zero_values(a);

        % 根据不同的空间数值创建三角剖分
        fv = isosurface(x, y, z, niiData == space_value, 0.1); % 当前空间的三角剖分

        % 绘制当前空间的三角剖分，设置透明度为0.7，并去除边缘线
        p = patch(fv);
        set(p, 'FaceColor', color(a,:), 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    end

    hold off;
    view(95,85);
    % 添加光照效果
    light('Position',[0,0,1],'Style','infinite');
    material('dull');

    % 设置绘图属性
    axis equal;
    axis tight;

    daspect([1, 1, 1]);
    xlim([-50 200])
    ylim([0 250])

    %不显示坐标轴和网格
    axis off;
    grid off;

    % 设置图形的背景颜色为黑色
    set(gcf, 'Color', 'k');
    set(gcf, 'InvertHardcopy', 'off'); % 阻止导出时背景颜色反转
    filename2 = sprintf('%d.png', length(C));
    print(gcf, filename2, '-dpng', '-r600');
end
