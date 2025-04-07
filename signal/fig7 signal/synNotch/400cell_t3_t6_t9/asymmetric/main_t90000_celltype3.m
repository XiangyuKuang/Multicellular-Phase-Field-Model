clc;
clear;
load('asymmetric_400_contact.mat');

n = 20;
DATA=cell(n, 1)
a = 10;
signal_pathway=cell(n, 1);

for i = 1:n      
    % t=90000
    cell_types_0 = celltype';
    X6 =contact{10,1};
    [max_cell_types, cell_types_0, signal_pathway] = signal_transduction(cell_types_0, X6, a, i, 1, signal_pathway);  
    DATA{i,1} = max_cell_types;
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

asymmetric_400_one_signal_t90000=signal_pathway([1,3],:);