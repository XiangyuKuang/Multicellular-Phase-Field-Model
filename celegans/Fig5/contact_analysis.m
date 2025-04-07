% 1) move *_contact.mat of 100 simulations to one folder
%   noise_sigma0d2 for sn0d2_[1-100]/*_contact.mat
%   noise_v1 for vn1_[1-100]/*_contact.mat
% 2) change current folder to noise_sigma0d2 or noise_v1 folder

clear
load('../../CellParameters.mat','CellParameters')
% SimName = 'sn0d2';%noise_sigma0d2
SimName = 'vn1';%noise_v1
SimNum = 100;
%%
load('sn0_1_contact.mat')% baseline contact,
contact0 = contact; 
% load(sprintf("%s_%d_contact.mat",SimName,1))
load(sprintf("%s_%d_trajectory.mat",SimName,1))
CellIdx = trajectory{2};
[CellIdx,I] = sort(CellIdx);
CellName = CellParameters.Name(CellIdx);
CellNum = length(CellIdx);
t = trajectory{1}{CellNum+1};

nonrobust_contact = zeros(CellNum);
robust_contact = zeros(CellNum);
contact_stack = cell(SimNum,1);
contact_T = cell(SimNum,1);
Tn = length(t);
contact_diff = zeros(SimNum,Tn);
for simidx = 1:SimNum
  load(sprintf("%s_%d_contact.mat",SimName,simidx))
  
  contact_i = zeros(CellNum,CellNum);
  for tn = 1:Tn
    contact_i = contact_i+contact{1}{tn};
    contact_diff(simidx, tn) = sum(abs(contact{1}{tn} - contact0{1}{tn}),'all')/2;
  end
  contact_T{simidx} = contact{1}{tn}(I,I);
%   contact_i=contact_i(I,I);
%   robust_contact_i = contact_i==length(t)
  contact_stack{simidx} = contact_i(I,I);
  nonrobust_contact = nonrobust_contact+(contact_i>0 & contact_i<Tn);
  robust_contact = robust_contact+(contact_i==Tn);
end
robust_contact = robust_contact(I,I);
nonrobust_contact = nonrobust_contact(I,I);
%% Fig5a 2nd column
figure
for simidx = 1:SimNum
  semilogx(t-t(1),contact_diff(simidx,:)+0.1,'color',[0.5,0.5,0.8],'linewidth',0.8);hold on
end
semilogx(t-t(1),mean(contact_diff)+0.1,'k')
axis square
xticks([10,100,1000,10000])
ylabel('Contact Variation \eta_{\rmC}')
xlabel('{\it in silico} Time')
%% Fig5a 3rd column average area (final state)
load(sprintf('%s_contact_area.mat',SimName),'contact_area')
mean_contact_area = zeros(CellNum);
for simidx = 1:SimNum
  mean_contact_area = mean_contact_area+contact_area{simidx};
end
mean_contact_area = mean_contact_area(I,I)/SimNum;
Iv = mean_contact_area>0;
dot_size = [2,8];
s = rescale(mean_contact_area(Iv),dot_size(1),dot_size(2));
cx = rescale(robust_contact(Iv),0,1);
figure;hold on
[i,j] = find(Iv);
for n = 1:length(i)
  c = cx(n)*[1,0,0]+(1-cx(n))*[0,0,1];
  plot([i(n),j(n)],[j(n),i(n)],'o',MarkerSize=s(n),linewidth=0.5,MarkerFaceColor=c,MarkerEdgeColor=c)
end
% for n = 1:length(i)
%   if i(n)>j(n)
%     % robust_contact == SimNum  k
%     % robust_contact<SimNum & robust_contact>0  gray
%     % nonrobust_contact>0  w
%     if robust_contact(i(n),j(n)) == SimNum 
%       plot([i(n),j(n)],[j(n),i(n)],'ko',MarkerSize=s(n),MarkerFaceColor='k',linewidth=0.5)
%     elseif robust_contact(i(n),j(n)) >0
%       plot([i(n),j(n)],[j(n),i(n)],'ko',MarkerSize=s(n),MarkerFaceColor=[0.5,0.5,0.5],linewidth=0.5)
%     elseif nonrobust_contact(i(n),j(n)) >0
%       plot([i(n),j(n)],[j(n),i(n)],'ko',MarkerSize=s(n),linewidth=0.5)
%     end
%   end
% end
% imagesc(mean_duration)
set(gca,'ydir','reverse','box','off','xtick',1:CellNum,'ytick',1:CellNum,'xticklabels',CellName,'yticklabels',CellName)%,'
axis('square',[0,CellNum+1,0,CellNum+1])
cx = linspace(0,1,256)';
colormap(gca,cx*[1,0,0]+(1-cx)*[0,0,1])
c=colorbar(gca);
c.Ticks = [0,0.5,1];
c.TickLabels = [0,50,100];

legend_sizes = linspace(min(mean_contact_area(Iv)),max(mean_contact_area(Iv)),4);
legend_marker_sizes = rescale(legend_sizes,dot_size(1),dot_size(2));
legend_labels = arrayfun(@(x) sprintf(' %.1f', x), legend_sizes, 'UniformOutput', false);
x_pos = 0; 
y_pos = (0:length(legend_sizes) - 1) * 1.5; 
figure('unit','pixel','Position',[100,100,150,200]);hold on
for k = 1:length(legend_sizes)
    plot(x_pos, y_pos(k), 'ko', 'MarkerSize', legend_marker_sizes(k), 'MarkerFaceColor', 'k');
    text(x_pos + 0.1, y_pos(k), legend_labels{k}, 'VerticalAlignment', 'middle','FontSize',15);
end
axis off

figure('unit','pixel','Position',[100,100,200,150]);hold on
a=0.015;
ms_max = max(legend_marker_sizes);
ms_min = min(legend_marker_sizes);
s_max =  max(mean_contact_area(Iv));
f = @(s,ms_max,ms_min,a)(s-ms_min)/(ms_max-ms_min)*4.5;
if s_max>120
  ymax = y_pos(end);
  patch([-1.5,ymax,ymax],[0,0,a*ms_max],'k')
else
  ymax = f(120/s_max*ms_max,ms_max,ms_min);
  patch([-1.5,ymax,ymax],[0,0,a*ms_max*120/s_max],'k')
end
ylim([0,1])
set(gca,'tickdir','out','ycolor','none','TickLength',[0.02,0.02],'XTickLabelRotation',0)
xticks(f(linspace(0,120,4)/s_max*ms_max,ms_max,ms_min))
xticklabels(linspace(0,120,4))

%% Fig5a 4th column Area (final and inital)-Conservativeness 
figure
ia = mean_contact_area>0;% & mean_duration<1001;
plot(mean_contact_area(ia),robust_contact(ia),'k.','markersize',10)
axis('square')%,[-50,1050,-10,150]
ylabel('Conservativeness',FontSize=20)
xlabel('Mean Contact Area (μm^2)',FontSize=20)
ylim([-5,105])

load('../sn0_1/IC_contact_area.mat','contact_area_sort_IC')
% remove div interface
i = [2     1     4     3     6     5     8     7];
j = [1     2     3     4     5     6     7     8];
for n = 1:length(i)
  contact_area_sort_IC(i(n),j(n))=0;
end

ia = contact_area_sort_IC>0;% & mean_duration<1001;
plot(contact_area_sort_IC(ia),robust_contact(ia),'k.','markersize',10)
ylim([-5,105])
ylabel('Conservativeness',FontSize=20)
xlabel('Initial Contact Area (μm^2)',FontSize=20)
axis('square')