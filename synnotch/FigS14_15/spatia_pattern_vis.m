% revision 
% supplement figures for synnotch spatial pattern
% based on spatial_pattern_comp_dev.mlx
clear
enable_save_rphi = false; %1
enable_save_thetaphi = true; %2
enable_save_prop = false; %3
enable_save1R = false; %1R
enable_legend = false;
ct_name = 'ct5';
% ct_name = 'ct3';
% ct_name = 'ct2';
fs = 25;% theta-phi plot fontsize 
%% color map
cmap = zeros(8,3);
c_t = turbo(8);
c_t(1,:) = c_t(1,:)*2.5;
%% asym
% load('400PS/CellParameters_200200_1.mat','CellParameters')
% cmap(1:3,:) = [[0.5,1,0.5];[1,0,0];[0,0,1]];
% cmap(4:8,:) = c_t([1,3,5,6,8],:);
% repeat_list = [6,7,8,9];
% 
% switch ct_name
%   case 'ct5'
%     [ct_stack,celltype_struct] = load_cell_type('400PS/asymmetric_t3_t6_t9.mat');
%     ts = [30000,60000,90000];
%     ct_i = 6;
%     group1 = [2,5];
%     group1_str = arrayfun(@(x)sprintf('Cell Type %d',x),group1,'UniformOutput',false);
%     group2 = [1,3,4]; %1221
% %     group1 = [2,3]; %0103
% %     group2 = [1,4,5];
%     base_ct = 1;
%   case 'ct3'
%     [ct_stack,celltype_struct] = load_cell_type('400PS/asymmetric_400_one_signal_t90000.mat');
%     ts = 90000;
%     ct_i = 1;
%     group1 = [2,3];
%     group2 = 1;
%     base_ct = 1;
%   case 'ct2'
%     group1 = 2;
%     group2 = 1;
%     base_ct = 1;
% end
% 
% save_name_list = cell(length(repeat_list)+1,1);
% sim_name_list = cell(length(repeat_list)+1,1);
% trajectory_name_list = cell(length(repeat_list)+1,1);
% contact_name_list = cell(length(repeat_list)+1,1);
% 
% save_name_list{1} = '400PS';
% sim_name_list{1} = '400PS/NewmainPSdt1d25CT1_200200_2';
% trajectory_name_list{1} = '400PS/NewmainPSdt1d25CT1_200200_2_trajectory.mat';
% contact_name_list{1} = '400PS/asymmetric_400_contact.mat';
% for n = 1:length(repeat_list)
%   repeat = repeat_list(n);
%   save_name_list{n+1} = sprintf('400PS_%d',repeat);
%   sim_name_list{n+1} = sprintf('400PS_%1$d/NewmainPSdt1d25CT1_200200_%1$d',repeat);
%   trajectory_name_list{n+1} = sprintf('400PS_%1$d/NewmainPSdt1d25CT1_200200_%1$d_trajectory.mat',repeat);
%   contact_name_list{n+1} = sprintf('400PS_%1$d/NewmainPSdt1d25CT1_200200_%1$d_contact.mat',repeat);
% end

%% sym
load('400L2/CellParameters_120280_1.mat','CellParameters')
cmap(1:3,:) = [[0.5,1,0.5];[0,0,1];[1,0,0]]; %L2
cmap(4:8,:) = c_t([1,3,5,6,8],:);

repeat_list = 2:5;
switch ct_name
  case 'ct5'
    [ct_stack,celltype_struct] = load_cell_type('400L2/symmetric_t3_t6_t9.mat');
    ts = [30000,60000,90000];
    ct_i = 365;
    group1 = [2,3];
    group2 = [1,4,5];
    base_ct = 1;
  case 'ct3'
    [ct_stack,celltype_struct] = load_cell_type('400L2/symmetric_400_one_signal_t90000.mat');
    ts = 90000;
    ct_i = 2;
    group1 = [2,3];
    group2 = 1;
    base_ct = 1;
  case 'ct2'
    group1 = 1;
    group2 = 2;
    base_ct=1;
end

save_name_list = cell(length(repeat_list)+1,1);
sim_name_list = cell(length(repeat_list)+1,1);
trajectory_name_list = cell(length(repeat_list)+1,1);
contact_name_list = cell(length(repeat_list)+1,1);

save_name_list{1} = '400L2';
sim_name_list{1} = '400L2/NewmainL2dt1d25CT1_120280_1';
trajectory_name_list{1} = '400L2/NewmainL2dt1d25CT1_120280_1_trajectory.mat';
contact_name_list{1} = '400L2/symmetric_400_contact.mat';
for n = 1:length(repeat_list)
  repeat = repeat_list(n);
  save_name_list{n+1} = sprintf('400L2_%d',repeat);
  sim_name_list{n+1} = sprintf('400L2_%1$d/NewmainL2dt1d25CT1_120280_%1$d',repeat);
  trajectory_name_list{n+1} = sprintf('400L2_%1$d/NewmainL2dt1d25CT1_120280_%1$d_trajectory.mat',repeat);
  contact_name_list{n+1} = sprintf('400L2_%1$d/NewmainL2dt1d25CT1_120280_%1$d_contact.mat',repeat);
end

%% analysis all
for k = 1:length(trajectory_name_list)
  close all
  trajectory_name = trajectory_name_list{k};
  contact_name = contact_name_list{k};
  save_name = save_name_list{k};
  load(trajectory_name,'trajectory')
  load(contact_name,'contact')

  %% analysis
  CellIdx  = trajectory{2};
  t = trajectory{1}{end};
  x_t = trajectory{1}(1:end-1);
  CellNum = length(CellIdx);
  T = 90000;
  tn = find(t==T);
  x_T = zeros(CellNum,3);
  for n = 1:CellNum
    x_T(n,:) = x_t{n}(tn,:);
  end

  % sphere coords
  center =  mean(x_T, 1);
  relative_pos = x_T - center;
  r = vecnorm(relative_pos, 2, 2);  % 径向距离
  theta = acos(relative_pos(:,3) ./ r);  % 极角
  phi = atan2(relative_pos(:,2), relative_pos(:,1));  % 方位角

  % signal program
  ct = CellParameters.CellType;
  if ~strcmp(ct_name,'ct2')
    signal_n = size(celltype_struct,2);
    tc = cell2mat(contact(:,2)); % contact map time point
    ct_t = zeros(CellNum,signal_n+1);
    ct_t(:,1) = ct; % celltype at each signal time step
    for n = 1:signal_n
      contact_t = contact{tc == ts(n),1};
      sender = celltype_struct{ct_i,n}.idx1;
      receiver = celltype_struct{ct_i,n}.idx2;
      CellIdx_r = CellIdx(ct==receiver);
      I = any(contact_t(ct == sender,ct == receiver),1);
      ctm = max(ct);
      ct(CellIdx_r(I)) = ctm+1;
      ct_t(:,n+1) = ct;
    end
  end
  uct = unique(ct);
  nct = length(uct);

  %% r-phi
  % all
  figure('Position',[100 100 500 500]);hold on
  % polar grid
  rmax = max(abs(r))*0.9;
  c = 0.*ones(1,3);
  lw = 1;
  plot([0,0],[-rmax,rmax],'Color',c,'LineWidth',lw)
  plot([-rmax,rmax],[0,0],'Color',c,'LineWidth',lw)
  rticks = rmax*[1/3,2/3];
  tr = linspace(0,2*pi,400);
  for rr = rticks
    plot(rr*cos(tr), rr*sin(tr), 'Color',c, 'LineStyle','--','LineWidth',lw);
  end

  for n = 1:nct
    I = ct == uct(n);
    r_type = r(I);
    phi_type = phi(I);
    x = r_type.*cos(phi_type);
    y = r_type.*sin(phi_type);
    plot(x,y,'.','Color',cmap(uct(n),:),MarkerSize=20)
  end
  axis equal off
  if enable_save_rphi
    print(sprintf('%s%s_1.png',save_name,ct_name),'-dpng','-r800')
  end
  % group1
  c1 = 0.8*ones(1,3);
  figure('Position',[100 100 500 500]);hold on
  % polar grid
  rmax = max(abs(r))*0.9;
  lw = 1;
  plot([0,0],[-rmax,rmax],'Color',c,'LineWidth',lw)
  plot([-rmax,rmax],[0,0],'Color',c,'LineWidth',lw)
  rticks = rmax*[1/3,2/3];
  tr = linspace(0,2*pi,400);
  for rr = rticks
    plot(rr*cos(tr), rr*sin(tr), 'Color',c, 'LineStyle','--','LineWidth',lw);
  end

  Ivis = zeros(CellNum,1);
  for n = group1
    I = ct == uct(n);
    Ivis = Ivis+I;
  end
  r_type = r(~Ivis);
  phi_type = phi(~Ivis);
  x = r_type.*cos(phi_type);
  y = r_type.*sin(phi_type);
  plot(x,y,'.','Color',c1,MarkerSize=10)

  for n = group1
    I = ct == uct(n);
    r_type = r(I);
    phi_type = phi(I);
    x = r_type.*cos(phi_type);
    y = r_type.*sin(phi_type);
    plot(x,y,'.','Color',cmap(uct(n),:),MarkerSize=20)
  end

  axis equal off
  if enable_save_rphi
    % print(sprintf('%s%s_1(1.png',save_name,ct_name),'-dpng','-r800')
    exportgraphics(gca,sprintf('%s%s_1(1.png',save_name,ct_name),"Resolution",800)
  end
  % group2
  figure('Position',[100 100 500 500]);hold on
  % polar grid
  rmax = max(abs(r))*0.9;
  lw = 1;
  plot([0,0],[-rmax,rmax],'Color',c,'LineWidth',lw)
  plot([-rmax,rmax],[0,0],'Color',c,'LineWidth',lw)
  rticks = rmax*[1/3,2/3];
  tr = linspace(0,2*pi,400);
  for rr = rticks
    plot(rr*cos(tr), rr*sin(tr), 'Color',c, 'LineStyle','--','LineWidth',lw);
  end

  Ivis = zeros(CellNum,1);
  for n = group2
    I = ct == uct(n);
    Ivis = Ivis+I;
  end
  r_type = r(~Ivis);
  phi_type = phi(~Ivis);
  x = r_type.*cos(phi_type);
  y = r_type.*sin(phi_type);
  plot(x,y,'.','Color',c1,MarkerSize=10)

  for n = group2
    I = ct == uct(n);
    r_type = r(I);
    phi_type = phi(I);
    x = r_type.*cos(phi_type);
    y = r_type.*sin(phi_type);
    plot(x,y,'.','Color',cmap(uct(n),:),MarkerSize=20)
  end

  axis equal off
  if enable_save_rphi
    exportgraphics(gca,sprintf('%s%s_1(2.png',save_name,ct_name),"Resolution",800)
    % print(sprintf('%s%s_1(2.png',save_name,ct_name),'-dpng','-r800')
  end

  %% theta-phi
  % all
  figure;hold on
  for n = 1:nct
    I = ct == uct(n);
    r_type = r(I);
    phi_type = phi(I);
    theta_type = theta(I);
    plot(phi_type,theta_type,'.','Color',cmap(uct(n),:),MarkerSize=20)
  end
  xlabel('\eta');
  ylabel('\theta');
  yticks([0,pi/2,pi])
  yticklabels({'0','\pi/2','\pi'})
  xticks([-pi,-pi/2,0,pi/2,pi])
  xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
  offset = 0.06;
  axis('equal',[-pi-offset,pi+offset,0-offset,pi+offset])
  set(gca,'tickdir','out')
  if enable_save_thetaphi
    print(sprintf('%s%s_2.png',save_name,ct_name),'-dpng','-r800')
  end
  % group1
  figure;hold on
  Ivis = zeros(CellNum,1);
  for n = group1
    I = ct == uct(n);
    Ivis = Ivis+I;
  end
  phi_type = phi(~Ivis);
  theta_type = theta(~Ivis);
  plot(phi_type,theta_type,'.','Color',c1,MarkerSize=10)
  for n = group1
    I = ct == uct(n);
    phi_type = phi(I);
    theta_type = theta(I);
    plot(phi_type,theta_type,'.','Color',cmap(uct(n),:),MarkerSize=20)
  end

  xlabel('\psi');
  ylabel('\theta');
  yticks([0,pi/2,pi])
  yticklabels({'0','\pi/2','\pi'})
  xticks([-pi,-pi/2,0,pi/2,pi])
  xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
  offset = 0.06;
  axis('equal',[-pi-offset,pi+offset,0-offset,pi+offset])
  set(gca,'tickdir','out','fontsize',fs)
  if enable_save_thetaphi
    % print(sprintf('%s%s_2(1.png',save_name,ct_name),'-dpng','-r800')
    exportgraphics(gca,sprintf('%s%s_2(1.png',save_name,ct_name),"Resolution",800)
  end
  %legend
  if enable_legend
    figure('Color', 'w', 'Position', [100, 100, 200, 100]); hold on
    for n = group1
       plot(NaN, NaN, '.', 'Color', cmap(uct(n),:), 'MarkerSize', 20);
    end
    legend(group1_str, 'Location', 'northwestoutside', 'FontSize', 15,'box','off');
    axis off; 
  end
  % group2
  figure;hold on
  Ivis = zeros(CellNum,1);
  for n = group2
    I = ct == uct(n);
    Ivis = Ivis+I;
  end
  phi_type = phi(~Ivis);
  theta_type = theta(~Ivis);
  plot(phi_type,theta_type,'.','Color',c1,MarkerSize=10)
  for n = group2
    I = ct == uct(n);
    phi_type = phi(I);
    theta_type = theta(I);
    plot(phi_type,theta_type,'.','Color',cmap(uct(n),:),MarkerSize=20)
  end

  xlabel('\psi');
  ylabel('\theta');
  yticks([0,pi/2,pi])
  yticklabels({'0','\pi/2','\pi'})
  xticks([-pi,-pi/2,0,pi/2,pi])
  xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
  offset = 0.06;
  axis('equal',[-pi-offset,pi+offset,0-offset,pi+offset])
  set(gca,'tickdir','out','fontsize',fs)
  if enable_save_thetaphi
    exportgraphics(gca,sprintf('%s%s_2(2.png',save_name,ct_name),"Resolution",800)
    % print(sprintf('%s%s_2(2.png',save_name,ct_name),'-dpng','-r800')
  end

  %% r-proportion
  nb = 12;
  r_edges = linspace(0, max(r), nb + 1);
  r_centers = (r_edges(1:end-1) + r_edges(2:end)) / 2;
%   display(max(r))
  count_matrix = zeros(nb, nct);
  ct = reshape(ct,[CellNum,1]);
  for i = 1:nb
    % 找到在当前径向区间内的细胞
    in_bin = (r >= r_edges(i)) & (r < r_edges(i+1));
    % 统计每种类型的数量
    for j = 1:nct
      count_matrix(i, j) = sum(in_bin & (ct == uct(j)));
    end
  end
  proportion_matrix = count_matrix ./ sum(count_matrix, 2);
  proportion_matrix(isnan(proportion_matrix)) = 0;
  I = any(proportion_matrix'>0);
  figure('Position',[100 100 500 500])
  h = bar(r_centers(I), proportion_matrix(I,:), 'stacked', 'BarWidth', 1, 'EdgeColor', 'none');
  for j = 1:nct
    h(j).FaceColor = cmap(uct(j), :);
    h(j).EdgeColor = [0.5,0.5,0.5];  % 去掉边框使图更清晰
  end
  ylabel('Cell Type Proportion')
  xlabel('Distance from Center (Grid Step)')
  axis square tight
  % legend_labels = arrayfun(@(x) sprintf('%d', x), uct, 'UniformOutput', false);
  % legend(h,legend_labels, Location = 'eastoutside',Box = 'off')
  if enable_save_prop
    print(sprintf('%s%s_3.png',save_name,ct_name),'-dpng','-r800')
  end

  %% rotated r-phi
  % using the symmetry axis as coordiate axis
  I = ct == base_ct;
  inertia_tensor = relative_pos(I,:)'*relative_pos(I,:);
  [eigvecs, eigvals] = eig(inertia_tensor);
  [sorted_eigvals, sort_idx] = sort(diag(eigvals),'descend');%
  sorted_eigvecs = eigvecs(:, sort_idx);
  v = sorted_eigvecs(:, 2);
  xyv = [v(1); v(2); 0];
  xyv_norm = norm(xyv);
  if xyv_norm<1e-5
    R = eye(3);
    fprintf('对称轴已经垂直于xy平面\n');
  else
    xyv = xyv/xyv_norm;
    rotation_axis = [-v(2);v(1);0];% 旋转轴：xy投影方向的垂直向量（在xy平面内）
    rotation_axis  = rotation_axis /norm(rotation_axis );
    rotation_angle = -atan2(v(3), xyv_norm);% 需要旋转的角度
    K = [0, -rotation_axis(3), rotation_axis(2);
         rotation_axis(3), 0, -rotation_axis(1);
         -rotation_axis(2), rotation_axis(1), 0];
    R = eye(3) + sin(rotation_angle) * K + (1 - cos(rotation_angle)) * (K * K);
  end
%   R = [sorted_eigvecs(:, 1), sorted_eigvecs(:, 3), sorted_eigvecs(:, 2)]';
  
  relative_posR = (R*relative_pos')';
  rR = vecnorm(relative_posR, 2, 2);  % 径向距离
  thetaR = acos(relative_posR(:,3) ./ r);  % 极角
  phiR = atan2(relative_posR(:,2), relative_posR(:,1));  % 方位角

  figure('Position',[100 100 500 500]);hold on
  % polar grid
  rmax = max(abs(r))*0.9;
  c = 0.*ones(1,3);
  lw = 1;
  plot([0,0],[-rmax,rmax],'Color',c,'LineWidth',lw)
  plot([-rmax,rmax],[0,0],'Color',c,'LineWidth',lw)
  rticks = rmax*[1/3,2/3];
  tr = linspace(0,2*pi,400);
  for rr = rticks
    plot(rr*cos(tr), rr*sin(tr), 'Color',c, 'LineStyle','--','LineWidth',lw);
  end

  for n = 1:nct
    I = ct == uct(n);
    r_type = rR(I);
    phi_type = phiR(I);
    x = r_type.*cos(phi_type);
    y = r_type.*sin(phi_type);
    plot(x,y,'.','Color',cmap(uct(n),:),MarkerSize=20)
  end
  axis equal off

  % group1
  figure('Position',[100 100 500 500]);hold on
  % polar grid
  rmax = max(abs(r))*0.9;
  lw = 1;
  plot([0,0],[-rmax,rmax],'Color',c,'LineWidth',lw)
  plot([-rmax,rmax],[0,0],'Color',c,'LineWidth',lw)
  rticks = rmax*[1/3,2/3];
  tr = linspace(0,2*pi,400);
  for rr = rticks
    plot(rr*cos(tr), rr*sin(tr), 'Color',c, 'LineStyle','--','LineWidth',lw);
  end

  Ivis = zeros(CellNum,1);
  for n = group1
    I = ct == uct(n);
    Ivis = Ivis+I;
  end
  r_type = rR(~Ivis);
  phi_type = phiR(~Ivis);
  x = r_type.*cos(phi_type);
  y = r_type.*sin(phi_type);
  plot(x,y,'.','Color',c1,MarkerSize=10)

  for n = group1
    I = ct == uct(n);
    r_type = rR(I);
    phi_type = phiR(I);
    x = r_type.*cos(phi_type);
    y = r_type.*sin(phi_type);
    plot(x,y,'.','Color',cmap(uct(n),:),MarkerSize=20)
  end

  axis equal off
  if enable_save1R
    % print(sprintf('%s%s_1(1.png',save_name,ct_name),'-dpng','-r800')
    exportgraphics(gca,sprintf('%s%s_1R(1.png',save_name,ct_name),"Resolution",800)
  end

  % group2
  figure('Position',[100 100 500 500]);hold on
  % polar grid
  rmax = max(abs(r))*0.9;
%   lw = 1;
  plot([0,0],[-rmax,rmax],'Color',c,'LineWidth',lw)
  plot([-rmax,rmax],[0,0],'Color',c,'LineWidth',lw)
  rticks = rmax*[1/3,2/3];
  tr = linspace(0,2*pi,400);
  for rr = rticks
    plot(rr*cos(tr), rr*sin(tr), 'Color',c, 'LineStyle','--','LineWidth',lw);
  end

  Ivis = zeros(CellNum,1);
  for n = group2
    I = ct == uct(n);
    Ivis = Ivis+I;
  end
  r_type = rR(~Ivis);
  phi_type = phiR(~Ivis);
  x = r_type.*cos(phi_type);
  y = r_type.*sin(phi_type);
  plot(x,y,'.','Color',c1,MarkerSize=10)

  for n = group2
    I = ct == uct(n);
    r_type = rR(I);
    phi_type = phiR(I);
    x = r_type.*cos(phi_type);
    y = r_type.*sin(phi_type);
    plot(x,y,'.','Color',cmap(uct(n),:),MarkerSize=20)
  end

  axis equal off
  if enable_save1R
    exportgraphics(gca,sprintf('%s%s_1R(2.png',save_name,ct_name),"Resolution",800)
    % print(sprintf('%s%s_1(2.png',save_name,ct_name),'-dpng','-r800')
  end
end
close all