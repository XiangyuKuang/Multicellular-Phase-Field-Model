clc;
clear;
load('real_Row_20000.mat');
load('signal_pathway_20000.mat');
load('cell_type_20000.mat');

[R, C] = size(B);
meanValues = zeros(1, C);
stdValues = zeros(1, C);
for w = 1:C
    col = B(:, w);
    meanValues(w) = mean(col);
    stdValues(w) = std(col);
end

real_B=B(Row,:);
[Ra, Ca] = size(real_B);
meanValues_real= zeros(1, Ca);
stdValues_real = zeros(1, Ca);
for w = 1:Ca
    col = real_B(:, w);
    meanValues_real(w) = mean(col);
    stdValues_real(w) = std(col);
end

figure;
hold on;
x = 1:C;
h1 = plot(x, meanValues, '.', 'Color', [8, 126, 139]/255); 
h2 = plot(x, meanValues_real, '.', 'Color', [255, 90, 95]/255); 
h3 = plot(x, meanValues, 'LineWidth', 1, 'Color', [8, 126, 139]/255); 
h4 = plot(x, meanValues_real, 'LineWidth', 1, 'Color', [255, 90, 95]/255); 

e = errorbar(x, meanValues, stdValues, 'Color', [8, 126, 139]/255);
e.LineWidth = 0.5;  
e.CapSize = 0; 
e1 = errorbar(x, meanValues_real, stdValues_real, 'Color', [255, 90, 95]/255); 
e1.LineWidth = 0.5;  
e1.CapSize = 0; 
yline(7, '--b');

% Set legend and its font size
lgd = legend([h3 h4], {'All programs: 20000 simulations', ...
    'C.elegans Program: 386 simulations'}, ...
    'Location', 'north', 'Orientation', 'vertical');
set(lgd, 'FontName', 'Arial', 'FontSize', 11);  % Set legend font size

% Set axis label sizes
xlabel('Optimization Step', 'FontSize', 14);  
ylabel('Cell Type Number', 'FontSize', 14);   

% Set axis properties
set(gca, 'FontSize', 18, 'FontName', 'Arial');  % Set axis tick label size
ylim([2 9]);
set(gca, 'XLim', [0 120]);      
set(gca, 'XTick', 0:30:120);
set(gca, 'YLim', [2 10]);      
set(gca, 'YTick', 2:1:10);  

axis square;
hold off;