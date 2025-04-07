% 定义信号转导函数
function [max_cell_types, cell_types, signal_pathway] = signal_transduction(cell_types, X, max_iterations,i,col, signal_pathway)
max_cell_types = zeros(1, max_iterations);
initial_max_type = max(cell_types);

for iter = 1:max_iterations
    unique_types = unique(cell_types);
    cell_types_next = cell_types;
    signal_transduced = false;
    
    % 随机选择两组进行信号转导
    idx1 = randi(length(unique_types));
    idx2 = randi(length(unique_types));
    
    while idx1 == idx2
        idx2 = randi(length(unique_types));
    end
    
    group1 = find(cell_types == unique_types(idx1));
    group2 = find(cell_types == unique_types(idx2));
    
    for j = group2
        if any(X(group1, j) > 0)
            % 检查j的姐妹细胞是否与group1中的所有细胞无接触
            sisters = group2(group2 ~= j);
            if any(all(X(sisters, group1) == 0, 2))
                % 更新接收信号的细胞命运
                cell_types_next(j) = max(cell_types) + 1;
                signal_transduced = true;  % 设置标志，表示信号已传导
            end
        end
    end
    

    
    if signal_transduced
        cell_types = cell_types_next;
        signal_pathway{i, col} = struct('cell_types', cell_types, 'idx1', idx1, 'idx2', idx2);
    else signal_pathway{i, col} = struct('cell_types', cell_types, 'idx1', 0, 'idx2', 0);
    end
    
    % 记录细胞类型矩阵的最大值
    max_cell_types(iter) = max(cell_types);
    
    % 如果细胞类型数增加，提前退出
    if max_cell_types(iter) > initial_max_type
        break;
    end
end
end