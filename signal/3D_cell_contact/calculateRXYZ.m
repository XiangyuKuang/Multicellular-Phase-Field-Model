function [v, x, y, z, t_stage] = calculateRXYZ(C,Embryo)
v = zeros(17, length(C));
x = zeros(17, length(C));
y = zeros(17, length(C));
z = zeros(17, length(C));
row_count = 1;

for i = Embryo
    CD = readtable(['F:\data\CD\CD', num2str(i,'%02d'), '.csv']);
    CD_cellname = string(CD{:, 2});
    %计算时间

    for n = 1:length(C)
        r = find(strcmp(CD_cellname, C(n)));
        t{n} = CD{r, 3};
        tmax(n) = max(t{n});
    end
        t_stage(row_count) = min(tmax);

    % 细胞核位置
    for m = 1:length(C)
        xyz = CD{:, [3, 9, 10, 11]};
        [r, ~] = find((xyz(:, 1) == t_stage(row_count)) & (strcmp(CD_cellname, C(m))));
        xyz(:, 2:4) = xyz(:, 2:4) .* [0.42, 0.09, 0.09];
        z(row_count, m) = xyz(r, 2);
        x(row_count, m) = xyz(r, 3);
        y(row_count, m) = xyz(r, 4);
        % 读取体积
        volume = readtable(['F:\data\CShaper_1\Sample', num2str(i, '%02d'), '\Sample', num2str(i, '%02d'), '_volume.csv'], 'VariableNamingRule', 'preserve');
        col_data = volume.(C(m));
        v(i-3, m) = col_data(t_stage(row_count), 1);
    end

    if ismember(i, [4, 5, 8, 12, 19])
        y(row_count, :) = 256 - y(row_count, :);
        z(row_count, :) = 256 - z(row_count, :);
    end
    if ismember(i, [6, 10, 14])
        x(row_count, :) = 184 - x(row_count, :);
        z(row_count, :) = 256 - z(row_count, :);
    end
    if ismember(i, [13, 16, 17, 18, 20])
        x(row_count, :) = 184 - x(row_count, :);
        y(row_count, :) = 256 - y(row_count, :);
    end
    row_count = row_count + 1;
end

num_cells=length(C);
save(['F:\code\3D cell contact\', num2str(num_cells), '_x.mat'], 'x');
save(['F:\code\3D cell contact\', num2str(num_cells), '_y.mat'], 'y');
save(['F:\code\3D cell contact\', num2str(num_cells), '_z.mat'], 'z');
save(['F:\code\3D cell contact\', num2str(num_cells), '_v.mat'], 'v');
save(['F:\code\3D cell contact\', num2str(num_cells), '_t.mat'], 't_stage');
end