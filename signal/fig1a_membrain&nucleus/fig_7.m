clear;
load('E:\data\embryo07_cell_stage_time.mat');
for a=1:7
time=embryo07_cell_stage_time(a);
% 读取细胞膜文件
filename = ['E:\data\CShaper_1\Sample07\RawMemb\Sample07_', num2str(time,'%03d'), '_rawMemb.nii.gz'];
niiData = niftiread(filename);
% 读取细胞核文件
filename1 = ['E:\data\RawNuc\Sample07\Sample07_',num2str(time,'%03d'),'_rawNuc.nii.gz'];
niiData1= niftiread(filename1);

% 获取矩阵的大小
[m, n, p] = size(niiData);

% 创建一个网格
[x, y, z] = meshgrid(1:n, 1:m, 1:p);

% 筛选第一组大于150的元素
data=80
filtered_values_group1 = niiData(niiData > data);
filtered_x_group1 = x(niiData > data);
filtered_y_group1 = y(niiData > data);
filtered_z_group1 = z(niiData > data);


% 筛选第二组大于15的元素
data1=10
filtered_values_group2 = niiData1(niiData1 > data1);
filtered_x_group2 = x(niiData1 >data1);
filtered_y_group2 = y(niiData1 >data1);
filtered_z_group2 = z(niiData1 >data1);


% 绘制两组三维荧光图（红色和绿色）
figure;

% 假设 radius_values_group1 是一个与数据点对应的半径向量
radius_values_group = ones(size(filtered_x_group1)) * 0.2;
radius_values_group1 = ones(size(filtered_x_group2)) * 0.3;

scatter3(filtered_x_group1(:), filtered_y_group1(:), filtered_z_group1(:), radius_values_group, 'r', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.12);
% 细胞膜

hold on;

% 细胞核
scatter3(filtered_x_group2(:), filtered_y_group2(:), filtered_z_group2(:),radius_values_group1,'g', 'filled','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerFaceAlpha',.071,'MarkerEdgeAlpha',.071);
hold off;

daspect([1, 1, 1]);
xlim([0 300])
ylim([0 200])
% 设置视角为默认视角（正三视图）
view(0,90);

% 不显示坐标轴和网格
axis off;
grid off;

% 设置图形的背景颜色为黑色
set(gcf, 'Color', 'k');
set(gcf, 'InvertHardcopy', 'off'); % 阻止导出时背景颜色反转
filename2 = sprintf('%d.png', a);
print(gcf, filename2, '-dpng', '-r600');
end