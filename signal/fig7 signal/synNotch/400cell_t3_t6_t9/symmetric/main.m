clc;
clear;
load('symmetric_400_contact.mat');

r = 20000;
DATA=cell(r, 3)
a = 10;
signal_pathway=cell(r, 3);

for i = 1:r
    % t=30000
    cell_types3_0 = celltype';
    X3 =contact{4,1} ;
    [max_cell_types_3, cell_types3_0, signal_pathway] = signal_transduction(cell_types3_0, X3, a, i, 1, signal_pathway);
    
    % t=60000
    cell_types4_0 = cell_types3_0;
    X4 =contact{7,1};
    [max_cell_types_4, cell_types4_0, signal_pathway] = signal_transduction(cell_types4_0, X4, a, i, 2, signal_pathway);
        
    % t=90000
    cell_types6_0 = cell_types4_0;
    X6 =contact{10,1};
    [max_cell_types_6, cell_types6_0, signal_pathway] = signal_transduction(cell_types6_0, X6, a, i, 3, signal_pathway);  
    
    %DATA{i,1} = max_cell_types_2;
    DATA{i,1} = max_cell_types_3;
    DATA{i,2} = max_cell_types_4;
    DATA{i,3} = max_cell_types_6;

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



