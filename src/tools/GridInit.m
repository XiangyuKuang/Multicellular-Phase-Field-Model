function [x,k,dr] = GridInit(L,N,EnableGPU)
dim = length(L);
x = cell(dim,1);
k = cell(dim,1);
dr = L./N;
% periodic, phi(0)=phi(L)
switch dim 
  case 3
    [x{1},x{2},x{3}] = meshgrid(0:dr(1):L(1)-dr(1),0:dr(2):L(2)-dr(2),0:dr(3):L(3)-dr(3));
%     [x{1},x{2},x{3}] = meshgrid(-L(1)/2:dr(1):L(1)/2-dr(1),-L(2)/2:dr(2):L(2)/2-dr(2),-L(3)/2:dr(3):L(3)/2-dr(3));
    [kx,ky,kz]=meshgrid(fftshift(-floor(N(1)/2):ceil(N(1)/2)-1),fftshift(-floor(N(2)/2):ceil(N(2)/2)-1),fftshift(-floor(N(3)/2):ceil(N(3)/2)-1));
    k{1} = 2*pi/L(1)*1i*kx;
    k{2} = 2*pi/L(2)*1i*ky;
    k{3} = 2*pi/L(3)*1i*kz;
  case 2
    [x{1},x{2}] = meshgrid(0:dr(1):L(1)-dr(1),0:dr(2):L(2)-dr(2));
%     [x{1},x{2}] = meshgrid(-L(1)/2:dr(1):L(1)/2-dr(1),-L(2)/2:dr(2):L(2)/2-dr(2));
    [kx,ky]=meshgrid(fftshift(-floor(N(1)/2):ceil(N(1)/2)-1),fftshift(-floor(N(2)/2):ceil(N(2)/2)-1));
    k{1} = 2*pi/L(1)*1i*kx;
    k{2} = 2*pi/L(2)*1i*ky;
  case 1
    x{1} = -L(1)/2:dr(1):L(1)/2-dr(1);
    k{1} = 2*pi/L(1)*1i*(fftshift(-floor(N(1)/2):ceil(N(1)/2)-1));
end
if nargin>2
  if EnableGPU
    x = cellfun(@gpuArray,x,'UniformOutput',0);
    k = cellfun(@gpuArray,k,'UniformOutput',0);
  end
end