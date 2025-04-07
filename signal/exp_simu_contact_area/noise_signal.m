clear;
close all;
a=5;
%load('sigma_noise_contact_area.mat');
load('velocity_noise_contact_area.mat');
cell_type0=readmatrix('F:\code\exp_simu_contact_area\initial_cell_type_12_cell_stage.xlsx');
file_path = sprintf('F:\\code\\exp_simu_contact_area\\group%d.xlsx', a);  
data = readtable(file_path);  
group1 = table2array(data);  

cell_types12_0=cell_type0(a,:);
m=size(group1,1)/2;
result=zeros(m,1);
num = zeros(1, m);

for n = 3
    ligand = group1(2*n - 1, :);
    receptor = group1(2*n, :);
    valid_ligand = ~isnan(ligand);
    valid_receptor = ~isnan(receptor);
    ligand = ligand(valid_ligand);
    receptor = receptor(valid_receptor);
    ligand_types = cell_types12_0(ligand(1));
    receptor_types = cell_types12_0(receptor(1));

    for i = 1:100
        X = contact_area_sort{i, 1};
        X_ligand_receptor = X(ligand, receptor);
        if all(any(X_ligand_receptor> 0, 1))
            all_ligand = find(cell_types12_0 == ligand_types);
            all_receptor = find(cell_types12_0 == receptor_types);
            sisters = all_receptor;
            sisters = setdiff(sisters, receptor(valid_receptor));
            if all(all(X(all_ligand, sisters) == 0, 2))
                num(n) = num(n) + 1;
            end
        end
    end
    result(n,1) = num(n) / 100;
end


