clear
load('../CellParameters.mat','CellParameters')
signalname = 'real_program_signal_m';
% signalname = 'all_program_signal_m';
% signalname = 'no_signal_m';
savename = signalname;
load(sprintf('%s.mat',signalname))
signal = eval(signalname);
%%
% set color map
n_state = length(signal);
cell_types_u = cell(1,n_state);
for n = 1:n_state
  cell_types = signal{n}.cell_types;
  cell_types_u{n} = unique(cell_types);
end
% c = turbo(max(cell2mat(cell_types_u)));

% c = turbo(8);% for real_program_signal, use the same color as all_program_signal
% c(1,:) = c(1,:)*2.5;

c = [0.752941176	0	0
0.3	0.3	0.3
0.439215686	0.188235294	0.62745098
0.184313725	0.458823529	0.71372549
0.717647059	0.035294118	0.670588235
0.298039216	0.48627451	0.17254902
1	0.752941176	0
0    0    0.5];
%%
% 1-8
dirname='../Fig2_S3/1to8dx0d25dt0d1adh_1';
[~, idx, FileList, ~]=GetTime(dirname);
MatList = cellfun(@(x)sprintf('%s/%s',dirname,x),FileList(idx),'UniformOutput',0);
% 12 
MatList = cat(1,MatList,{'../Fig2_S3/8adhto12dx0d25dt0d1c18sigma0d5-0_3_44700.mat'});
%%
dx=0.25;
th=0.25;
for n = 7% 1:n_state
  load(MatList{n+1},'phi_save','CellIdx','SimParameters')
  cell_types = signal{n}.cell_types;% order of cell_types is the same as sorted CellIdx
  % order cell_type by CellIdx
  [~,Ib]=ismember(CellIdx,sort(CellIdx)); % Ib recover sorted CellIdx back to CellIdx: CellIdx == CellIdx_s(Ib)
  cellcolor = c(cell_types(Ib),:);
  figure
  set(gcf,'unit','pixels','position',[0,0,2816/3,960]);
  MultiphaseDisp(phi_save,CellIdx,dx,CellParameters,'threshold',th,'sigma',SimParameters.sigma,'TripleView','revydir','color',cellcolor);
  state = length(CellIdx);
  print(gcf,sprintf('%s_%d(3.png',savename,state),'-dpng','-r800')
  close all
end

figure
x = [0,1,1,0];
y = [0,0,1,1];
for n = 1:size(c,1)
  patch(x+n,y,c(n,:),'EdgeColor','w','linewidth',2)
end
axis equal
axis off