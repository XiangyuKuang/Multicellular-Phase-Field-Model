clear
dirname = '400L2';
save_dirname = dirname;%'temp';
cellparaname = 'CellParameters_120280_1.mat';
load(sprintf('%s/%s',dirname,cellparaname),'CellParameters')
c = zeros(8,3);
c_t = turbo(8);
c_t(1,:) = c_t(1,:)*2.5;
c(1:3,:) = [[0.5,1,0.5];[0,0,1];[1,0,0]];
c(4:8,:) = c_t([1,3,5,6,8],:);
cell_type_number = 2; % cell type number in Fig7 

%% load cell type
if cell_type_number == 2
  M = 0;
  signal = {CellParameters.CellType'};
  skip_nosignal = 0;
  signal_n = zeros(1,10);
elseif cell_type_number == 3 
  load(sprintf('%s/%s.mat',dirname,'symmetric_400_one_signal_t90000'),'symmetric_400_one_signal_t90000')
  M = [1,2];
  signal = cell(2,size(symmetric_400_one_signal_t90000,2));
  for m = 1:2
    for n = 1:size(symmetric_400_one_signal_t90000,2)
      signal{m,n} = symmetric_400_one_signal_t90000{m,n}.cell_types;
    end
  end
  skip_nosignal = 1;
  signal_n = [0,0,0,0,0,0,0,0,1,1];% signal - timepoint
elseif cell_type_number == 5 
  load(sprintf('%s/%s.mat',dirname,'symmetric_t3_t6_t9'),'five_cell_type')
  M = [1,2,3,4,5,6,93,162,363,365,371,405,419,578, 647, 696];
  signal = cell(length(M),size(five_cell_type,2));
  for m = 1:length(M)
    for n = 1:size(five_cell_type,2)
      signal{m,n} = five_cell_type{M(m),n}.cell_types;
    end
  end
  skip_nosignal = 1;
  signal_n = [0,0,1,1,1,2,2,2,3,3];% signal - timepoint
end

%% mask
% mask = [];
mask = ones([150,150,150]);
mask(76:150,1:75,76:150) = 0;

%% visualization
synnotch_signal_vis(dirname,save_dirname,signal,M,signal_n,CellParameters,c,mask,skip_nosignal)
if vis_IC
  load('IC/N400r5dx0d5_1_0.mat','phi_save','CellIdx')
  synnotch_signal0_vis(phi_save,CellIdx,save_dirname,signal,M,CellParameters,c,mask)
end
%% colorbar
figure
x = [0,1,1,0];
y = [0,0,1,1];
for n = 1:size(c,1)
  patch(x+n,y,c(n,:),'EdgeColor','w','linewidth',2)
end
axis equal
axis off
print(gcf,sprintf('colorbar.png'),'-dpng','-r800')

%% functions
function synnotch_signal_vis(load_dirname,save_dirname,signal,signal_m,signal_n,CellParameters,cmap,mask,skip_nosignal)
% savename = dirname;
[time, idx, FileList, ~]=GetTime(load_dirname);
dr = 0.5;
viewVec = [-54.24,27.96];%[-28.95,12.05];
for n = 1:length(signal_n)
  if skip_nosignal
    if signal_n(n) == 0
      continue
    end
  end
  
  load(sprintf('%s/%s',load_dirname,FileList{idx(n)}),'phi_cropped','CellIdx','box')% load phi
  %   phi_save = cropped_phi_recovery(phi_cropped,box,[150,150,150]);
  if ~isempty(mask)
    phi_cropped = mask_phi(phi_cropped,box,mask);
    I = cellfun(@(x)any(x>0.34,'all'),phi_cropped);
    phi_cropped = phi_cropped(I);
    CellIdx = CellIdx(I);
    box = box(I,:);
  end
  for m = 1:size(signal,1)
    CellParametersVis = CellParameters;
    if signal_n(n) >0
      CellParametersVis.CellType = signal{m,signal_n(n)}';
    end
    CellParametersVis.Color = cmap(CellParametersVis.CellType,:);
    figure('unit','pixel','position',[0,0,630,540],'color',[1,1,1],'toolbar','none');
    MultiphaseDisp(phi_cropped,CellIdx,dr,CellParametersVis,'LegendOff','view',viewVec,'AddBox','revydir','box',box);
    axis off
    camproj perspective
    print(gcf,sprintf('%s/%s_%d_%d.png',save_dirname,load_dirname,signal_m(m),time(n)),'-dpng','-r800')
    close all
  end
  clearvars phi_save
  
end
end

function synnotch_signal0_vis(phi_save,CellIdx,dirname,signal,signal_m,CellParameters,cmap,mask)
dr = 0.5;
figname = '%s/%s_%d_%d.png';
viewVec = [-54.24,27.96];
if ~isempty(mask)
  for n = 1:length(CellIdx)
    phi_save{n} = mask.*phi_save{n};
  end
  figname = '%s/%s_%d_%d(1.png';
end
for m = 1:size(signal,1)
  CellParametersVis = CellParameters;
  CellParametersVis.CellType = signal{m,1}';
  CellParametersVis.Color = cmap(CellParametersVis.CellType,:);
  figure('unit','pixel','position',[0,0,630,540],'color',[1,1,1],'toolbar','none');
  MultiphaseDisp(phi_save,CellIdx,dr,CellParametersVis,'LegendOff','view',viewVec,'AddBox','revydir');
  axis off
  camproj perspective
  print(gcf,sprintf(figname,dirname,dirname,signal_m(m),0),'-dpng','-r800')
  close all
end
end