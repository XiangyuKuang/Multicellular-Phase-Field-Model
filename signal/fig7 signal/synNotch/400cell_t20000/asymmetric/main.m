clc;
clear;
load('asymmetric_400_contact.mat');

n = 20000;
DATA=cell(n, 5)
a = 10;
signal_pathway=cell(n, 5);

for i = 1:n
    % t=0
    % X2 =contact{1,1};
    % cell_types2_0 =celltype';
    % [max_cell_types_2, cell_types2_0, signal_pathway] = signal_transduction(cell_types2_0, X2, a, i, 1, signal_pathway);
    
    % t=10000
    cell_types3_0 = celltype';
    X3 =contact{3,1} ;
    [max_cell_types_3, cell_types3_0, signal_pathway] = signal_transduction(cell_types3_0, X3, a, i, 1, signal_pathway);
    
    % t=20000
    cell_types4_0 = cell_types3_0;
    X4 =contact{5,1};
    [max_cell_types_4, cell_types4_0, signal_pathway] = signal_transduction(cell_types4_0, X4, a, i, 2, signal_pathway);
    
    
    % t=30000
    cell_types6_0 = cell_types4_0;
    X6 =contact{7,1};
    [max_cell_types_6, cell_types6_0, signal_pathway] = signal_transduction(cell_types6_0, X6, a, i, 3, signal_pathway);
    
    % t=40000
    cell_types7_0 = cell_types6_0;
    X7 =contact{9,1};
    [max_cell_types_7, cell_types7_0, signal_pathway] = signal_transduction(cell_types7_0, X7, a, i, 4, signal_pathway);
    
    % t=50000
    cell_types8_0 = cell_types7_0;
    X8 =contact{11,1}; 
    [max_cell_types_8, cell_types8_0, signal_pathway] = signal_transduction(cell_types8_0, X8, a, i, 5, signal_pathway);

    %DATA{i,1} = max_cell_types_2;
    DATA{i,1} = max_cell_types_3;
    DATA{i,2} = max_cell_types_4;
    DATA{i,3} = max_cell_types_6;
    DATA{i,4} = max_cell_types_7;
    DATA{i,5} = max_cell_types_8;
    if mod(i,100)==0
    disp(i);
    end
end

B=[];
B=cell2mat(DATA);
[R, C] = size(B);
for i = 1:R
    row = B(i, :);
    b=max(row);
    nonZero = row(row ~= 0);
    num = C - length(nonZero);
    newRow = [nonZero, b*ones(1, num)];
    B(i, :) = newRow;
end
asymmetric_400_allsignal=signal_pathway(2,:);