clc;
clear;
load('all.mat');
%%数据提取画三维命运分布图
Maxcelltype = zeros(20000, 1); 
for r=1:20000
Maxcelltype(r,1)=max(signal_pathway{r, 3}.cell_types);
end

% 找到 Maxcelltype 中所有值为 7 的元素的位置
indices = find(Maxcelltype == 5);

%最终命运数为8的命运分布
five_cell_type = signal_pathway(indices,:); 
final_five_cell_type = cell(size(five_cell_type));  
for i = 1:size(five_cell_type, 1)
    for col=1:3
    final_five_cell_type{i,col} = signal_pathway{i,col}.cell_types;  % 提取第6列的 cell_types 字段
    end
end

%将该命运整理为从1开始
for r=1:size(final_seven_cell_type,1)
    for c=1:size(final_seven_cell_type,2)
        Cell_types = final_seven_cell_type{r,c};
        [unique_cell_types, ~, ~] = unique(Cell_types, 'stable');
        % 初始化一个空数组来存放结果
        sorted_cell_types = zeros(size(Cell_types));
        for i = 1:length(unique_cell_types)
            % 找出当前细胞类型的所有位置
            indices = find(Cell_types == unique_cell_types(i));
            % 将当前细胞类型的所有元素按原始顺序添加到结果数组
            sorted_cell_types(1,indices) = i;
        end
        final_seven_cell_type{r,c}=sorted_cell_types;
    end
end

new_array = {};
indices = [];
for i = 1:size(final_seven_cell_type, 1)
    current_value = final_seven_cell_type{i, 6};

    % 检查该值是否与其他行的第 6 列值相同
    is_unique = true; 
    for j = 1:size(final_seven_cell_type, 1)
        % 如果 i != j 并且当前值与其他行的第 6 列值相同
        if i ~= j && all(current_value == final_seven_cell_type{j, 6})
            is_unique = false;  
            break;
        end
    end

    % 如果该行的值与其他行不同，则将该行的全部数据存入新数组
    if is_unique
        new_array = [new_array; final_seven_cell_type(i, :)]; 
        indices = [indices, i];
    end
    % if length(indices)>20
    %     break
    % end
end

symmetric_400_allsignal = cell(20,6);  
symmetric_400_allsignal(1:20, :) = new_array(1:20, :);