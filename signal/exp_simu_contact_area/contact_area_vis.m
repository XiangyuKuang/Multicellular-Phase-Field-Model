clear
load('velocity_noise_contact_area.mat','contact_area_sort','CellName')
SimNum = length(contact_area_sort);
CellNum = length(CellName);
mean_contact_area = zeros(CellNum);
for simidx = 1:SimNum
  mean_contact_area = mean_contact_area+contact_area_sort{simidx};
end
mean_contact_area = mean_contact_area/SimNum;

dot_size = [0.5,8];%range of dot sizes
I = mean_contact_area>0;
s = rescale(mean_contact_area(I),dot_size(1),dot_size(2));
figure;hold on
[i,j] = find(I);
for n = 1:length(i)
  if i(n)>j(n)
    plot([i(n),j(n)],[j(n),i(n)],'ko',MarkerSize=s(n),MarkerFaceColor='k')
  end
end

set(gca,'ydir','reverse','box','off','xtick',1:CellNum,'ytick',1:CellNum,'xticklabels',CellName,'yticklabels',CellName)
axis('equal',[0,CellNum+1,0,CellNum+1])

%% legend
legend_sizes = linspace(min(mean_contact_area(I)),max(mean_contact_area(I)),4);
legend_marker_sizes = rescale(legend_sizes,dot_size(1),dot_size(2));
legend_labels = arrayfun(@(x) sprintf('%.1f', x), legend_sizes, 'UniformOutput', false);

x_pos = 0; 
y_pos = (0:length(legend_sizes) - 1) * 1.5; 
figure;hold on
for k = 1:length(legend_sizes)
    plot(x_pos, y_pos(k), 'ko', 'MarkerSize', legend_marker_sizes(k), 'MarkerFaceColor', 'k');
    text(x_pos + 0.1, y_pos(k), legend_labels{k}, 'VerticalAlignment', 'middle','FontSize',15);
end
axis off