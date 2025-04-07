clear;
close all;
load('wt_contact_area.mat');
load('F:\video\multiphase\video_code\video_code\CellParameters.mat');
load('F:\data\cell_contact.mat');

experiment_contact = Cell_Cell_Contact;
simulation_contact = contact_area_sort';

a = [1:7];
figure;
hold on;

colors = jet(78);
X = [];
Y = [];

colorIndex = 1;

% For each pair of cells
for i = a
    expe = experiment_contact{1, i};
    simu = simulation_contact{1, i};
    name = CellName{i, 1};
    c = length(name);
    
    % Loop over each pair of cells
    for m = 1:c
        for n = m+1:c
            if m ~= n && expe(m, n) ~= 0 && simu(m, n) ~= 0
                x = expe(m, n);
                y = simu(m, n);
               
                X(end+1) = x;
                Y(end+1) = y;
                
                % Assign a color to each unique (m,n) pair
                color = colors(mod(colorIndex - 1, size(colors, 1)) + 1, :); % Wrap color index if it exceeds the colormap

                % Plot the points as solid circles using scatter (with unique color for each point)
                scatter(x, y, 30, 'filled', 'MarkerFaceColor', color);
                hold on;
                
                % Increment color index
                colorIndex = colorIndex + 1;
            end
        end
    end
end

k = sum(X .* Y) / sum(X.^2);  

yFit = k * X;  
residuals = Y - yFit;
SSresidual = sum(residuals.^2);

meanY = mean(Y);
SStotal = sum((Y - meanY).^2);
R2 = 1 - (SSresidual / SStotal);

xF = linspace(min(X), max(X), 100);  
yF = k * xF;  
plot(xF, yF, '--k', 'LineWidth', 2);  % '--k' 表示黑色虚线
text(245.268456375839,262.0134228187919,sprintf('R^2'), 'FontSize', 16, 'Color', 'k','FontName', 'Arial','FontAngle', 'italic');
text(276.2248322147653,260.4446308724833,sprintf(' = %.4f', R2), 'FontSize', 16, 'Color', 'k','FontName', 'Arial');

%text(280, 260, sprintf('$R^2 = %.4f$', R2), 'FontSize', 12, 'Color', 'k','Interpreter', 'latex','FontName', 'Arial');
%text(280, 260, sprintf('\\it{R^2} = %.4f', R2), 'FontSize', 12, 'Color', 'k', 'Interpreter', 'tex', 'FontName', 'Arial');

% Customize the plot
set(gcf,'unit','centimeters','position',[10,5,16,12]);
set(gca,'unit','centimeters','position',[2.5,2,12,9]);

set(gca, 'Box', 'off', ...
     'TickDir', 'in', 'TickLength', [0 0], ...
     'XMinorTick', 'off', 'YMinorTick', 'off', ...
     'XTick', 0:150:450, ...
     'YTick', 0:100:300, ...
     'XLim', [0 450], ...
     'YLim', [0 300]);
set(gca, 'FontSize', 16, 'FontName', 'Arial');

% Adding labels to the axes
xlabel('Experimental Contact Area (μm^2)');
ylabel('Simulated Contact Area (μm^2)');
axis square;

filename1 = sprintf('all.png');
print(gcf, filename1, '-dpng', '-r600');
