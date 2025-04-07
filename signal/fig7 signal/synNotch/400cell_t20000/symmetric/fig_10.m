clc;
clear;
load('all.mat');

[R, C] = size(B);
meanValues = zeros(1, C);
stdValues = zeros(1, C);
for w = 1:C
    col = B(:, w);
    meanValues(w) = mean(col);
    stdValues(w) = std(col);
end

meanValues = [2, meanValues];
stdValues = [0, stdValues];
x = [0, 1:C]; 

figure;
hold on;
h1=plot(x, meanValues, '.', 'Color', [8, 126, 139]/255); 
h2=plot(x, meanValues, 'LineWidth', 1, 'Color', [8, 126, 139]/255); 

e = errorbar(x, meanValues, stdValues, 'Color', [8, 126, 139]/255); 
e.LineWidth = 0.5;  
e.CapSize = 0; 

set(gca, 'YTick', 1:1:10);
set(gca, 'XTick', 0:10:60);     
lgd = legend([h2], {'All Programs: 20000 simulations'},'Location', 'north', 'Orientation', 'vertical');
set(lgd, 'FontName', 'Arial', 'FontSize', 11);  % Set legend font size

xlabel('Optimization Step', 'FontSize', 14);  
ylabel('Cell Type Number', 'FontSize', 14);   

set(gca, 'FontSize', 18, 'FontName', 'Arial');  
ylim([1 8.5]);
xlim([0 60]);   
set(gca, 'XTick', 0:10:60);     
set(gca, 'YTick', 2:1:8);  

axis square;
hold off;