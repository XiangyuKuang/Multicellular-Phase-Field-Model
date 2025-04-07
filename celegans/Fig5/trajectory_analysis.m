% 1) move *_trajectory.mat of 100 simulations to one folder
%   noise_sigma0d2 for sn0d2_[1-100]/*_trajectory.mat
%   noise_v1 for vn1_[1-100]/*_trajectory.mat
% 2) change current folder to noise_sigma0d2 or noise_v1 folder

clear
CellNum = 12;
% load baseline trajectory
load('12_motif\8adhto12dx0d25dt0d1c18sigma0d5-0_3_trajectory.mat');%29930:10:end
b1 = trajectory{1,2}(1:CellNum); % baseline trajectory
lb1 = length(b1{1}); % baseline trajectory length
t1 = trajectory{1,2}{CellNum+1};
t1 = t1'-t1(1);
idx1e4 = find(t1==1e4);
%%
f1=figure;
f1.Position = [573 530 560 450];
% axes('Position',[0.19,0.23,0.79,0.75]);
%%
% diff between noise simulation trajectory and baseline trajectory
% phi
% SimNum = 10;
% SimName = 'phin5';

% velocity
SimNum = 100;
SimName = 'vn1';

% sigma
% SimNum = 100;
% SimName = 'sn0d2';

idx1st = 1;% start index of simulation
dist_mean1 = zeros(SimNum,lb1);
dist_end = zeros(SimNum,1);
for n = 1:SimNum
  load(sprintf('%s_%d_trajectory.mat',SimName,n+idx1st-1));
  xcell1 = trajectory{1,1}(1:CellNum); % noise sim trajectory
  l1 = min(length(xcell1{1}),lb1); % trajectory length
  dist_cell1 = zeros(CellNum,l1);
  for m = 1:CellNum
    dist_cell1(m,:) = vecnorm((xcell1{m}(1:l1,:)-b1{m}(1:l1,:))',2);
  end
  dist_mean1(n,:) = mean(dist_cell1);
  dist_end(n) = dist_mean1(n,idx1e4);

  loglog(t1,dist_mean1(n,:),'color',[0.5,0.5,0.8],'linewidth',0.8);hold on;%c(n,:)semilogx
%   semilogx(t1,dist_mean1(n,:),'color',[0.5,0.5,0.8],'linewidth',0.8);hold on;%c(n,:)
end
%%
dist_allmean1 = mean(dist_mean1);
t3 = t1(2:end);% (t1(2:end)+t1(1:end-1))'/2;
dist_allmean1_t = diff(dist_allmean1)/10;
dist_allmean1_tau = dist_allmean1_t.*t3*log(10);

th = dist_allmean1_tau(1)*5;
[~,ith] = find(abs(dist_allmean1_tau-th)<1e-3,1);
tth = t3(ith);

loglog(t1,dist_allmean1,'k')
set(gca,'fontsize',18)
ylabel('Position Variation \eta_{\rmP} (\mum)')
xlabel('{\it in silico} Time')
xticks([10,100,1000,10000])
% ylim([0.0001,0.005])
% xlim([1,10000])
box off
axis square
%%
[~,I]=sort(dist_end);
fprintf('%d,', I(end-10:end)')