clear;
close all;
load('F:\data\cellname.mat');
load('F:\data\cell_contact.mat');
load('F:\data\color_cell_labeled.mat');
load('F:\data\embryo.mat');

for A=7%选择细胞期
    C=string(Cellname{1,A}); %选择细胞名称
    Embryo=embryo{A}; %选择合适的胚胎编号
    Colors=colors{A}; %选择各细胞对应颜色

    %选择细胞接触图谱，与colorbar对应
    Matrix=Cell_Cell_Contact{A};
    maxM = max(Matrix(:));
    minM = min(Matrix(Matrix > 0));
    maxM_1=ceil(maxM./100)*100;
    minM_1=floor(minM./10)*10;
    numColors = 256;
    cmap = [linspace(0, 0.8, numColors)', linspace(0, 0.8, numColors)', linspace(0, 0.8, numColors)'];

    %[v, x, y, z, t_stage] = calculateRXYZ(C,Embryo);
    num_cells = length(C);
    load(['F:\code\3D cell contact\', num2str(num_cells), '_t.mat']);
    load(['F:\code\3D cell contact\', num2str(num_cells), '_x.mat']);
    load(['F:\code\3D cell contact\', num2str(num_cells), '_y.mat']);
    load(['F:\code\3D cell contact\', num2str(num_cells), '_z.mat']);
    load(['F:\code\3D cell contact\', num2str(num_cells), '_v.mat']);
    %将坐标移动到0附近
    x0 = x - mean(x, 2);
    z0 = y - mean(y, 2);
    y0 = -z + mean(z, 2);
    %取多个胚胎的平均值坐标/体积，体积化为半径*系数
    for n = 1:length(C)
        V(n) = mean(v(v(:, n) ~= 0, n));
        X(n) = mean(x0(x0(:, n) ~= 0, n));
        Y(n) = mean(y0(y0(:, n) ~= 0, n));
        Z(n) = mean(z0(z0(:, n) ~= 0, n));
    end
    radius = V.^(1/3)*0.12;
    figure(1);
    hold on;
    %调整材料和光源
    light('Position', [-2 -1 1], 'Style', 'infinite');
    material dull;
    lighting gouraud;
    
    set(gcf,'unit','centimeters','position',[10,8,17,11]);
    set(gca,'unit','centimeters','position',[2,2.5,12,8]);
    set(gca, 'Box', 'off', ...
        'TickDir', 'in', 'TickLength', [0 0], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', ...
        'XTick', -25:10:25, ...
        'XLim', [-25 25], ...
        'YTick', -8:8:8, ...
        'YLim', [-10 10], ...
        'ZTick', -10:10:10,...
        'ZLim', [-14 14]);
    %记录legend
    h = gobjects(length(C), 1);
    %画细胞核小球
    for i = 1:length(C)
        s(i)=scatter3(X(i), Y(i), Z(i), 10, 'filled', 'MarkerFaceColor', Colors(i, :),'DisplayName', C{i});
        [sx, sy, sz] = sphere(100);
        h(i)= surf(X(i) + radius(i) * sx, Y(i) + radius(i) * sy, Z(i) + radius(i) * sz, ...
            'FaceColor', Colors(i, :), 'FaceAlpha', 1, 'EdgeColor', 'none');
        set(h(i), 'AmbientStrength', 0.4, 'DiffuseStrength', 0.5, 'SpecularStrength', 0.3);
    end
    %画细胞核连线，并与接触矩阵对应，由小到大颜色从黑变灰，由细变粗
    for a = 1:length(C)
        for b = 1:length(C)
            if (b > a) && (Matrix(a, b) ~= 0)
                lineWidth = Matrix(a, b);
                normalizedWidth = (lineWidth - minM_1) / (maxM_1 - minM_1);
                colorIndex = round(normalizedWidth * (numColors - 1)) + 1;
                lineColor = cmap(colorIndex, :);
                plot3([X(a), X(b)], [Y(a), Y(b)], [Z(a), Z(b)], 'Color', lineColor, 'LineWidth', lineWidth*0.02);
            end
        end
    end

    %坐标轴设置&图例
    legend(h, C);
    if length(C)==12
    set(legend,'unit','centimeters','position',[14.7,3,2,0.525*length(C)]);
    end
    if length(C)==8
    set(legend,'unit','centimeters','position',[14.7,4.1,2,0.525*length(C)]);
    end
    if length(C)==7
    set(legend,'unit','centimeters','position',[14.7,4.4,2,0.525*length(C)]);
    end
    if length(C)==6
    set(legend,'unit','centimeters','position',[14.7,4.8,2,0.525*length(C)]);
    end
    if length(C)==4
    set(legend,'unit','centimeters','position',[14.7,5.6,2,0.525*length(C)]);
    end
    if length(C)==2
    set(legend,'unit','centimeters','position',[14.7,6.3,2,0.525*length(C)]);
    end
    set(legend, 'Box', 'off');
   
    xlabel('{\itx} / A-P Axis (μm)', 'Position', [-0.019479790155017,-16.356157502937435,-15.174631575614], 'Rotation', 3);
    ylabel('{\ity} / D-V Axis (μm)','Position', [-28.581305322049147,-4.180705295824851,-13.436120781498744], 'Rotation', -72);    
    zlabel('{\itz} / L-R Axis (μm)','Position', [-28.950341361076525,16.215367272381172,0.049188056902608], 'Rotation', 90);
    set(gca,'FontSize',14,'Fontname', 'Arial');
    
    view(-7, 22);
    grid off;
    axis equal;
    set(gca, 'XLim', [-25 25], 'YLim', [-10 10], 'ZLim', [[-14 14]]);
    filename1 = sprintf('cell_labeled_%d.png',length(C));
    print(gcf, filename1, '-dpng', '-r600');

    figure(2);
    set(gcf, 'Units', 'centimeters', 'Position', [10, 4, 16, 5]);
    set(gca, 'Units', 'centimeters', 'Position', [3, 0, 12, 1]);
    set(gca, 'Box', 'off', ...
        'TickDir', 'in', 'TickLength', [0 0], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', ...
        'XTick', -15:5:15, ...
        'XLim', [-20 20], ...
        'YTick', -8:4:8, ...
        'YLim', [-1 1]);

    XData = [-15; 15; 15];
    YData = [0;1; -1];
    c = [0;1;1];
    patch(XData, YData, c,'FaceColor','interp','EdgeColor', 'none');
    colormap(cmap);
    text(-16,2.8, sprintf('%d', minM_1),'FontSize',14,'Fontname', 'Arial');
    text(13,2.8,sprintf('%d', maxM_1),'FontSize',14,'Fontname', 'Arial');
    text(4,2.8,sprintf('%.f', (maxM_1+minM_1)./3*2),'FontSize',14,'Fontname', 'Arial');
    text(-6,2.8,sprintf('%.f', (maxM_1+minM_1)./3*1),'FontSize',14,'Fontname', 'Arial');
    text(16.9,3.2,'(μm^2)','FontSize',14,'Fontname', 'Arial');

    grid off;
    axis equal;
    axis off;

    filename2 = sprintf('cell_labeled_%d_1.png',length(C));
    print(gcf, filename2, '-dpng', '-r600');
end