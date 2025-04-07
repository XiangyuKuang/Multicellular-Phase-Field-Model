function [phi,box] = crop_phi(phi_save,th)
arguments
  phi_save cell
  th double = 1e-3
end
N = length(phi_save);
phi = cell(N,1);
box = zeros(N,6);
sz = size(phi_save{1});
for n = 1:N
  i = find(phi_save{n}>=th);
  [i1,i2,i3] = ind2sub(sz,i);
  box(n,:) = [min(i1),max(i1),min(i2),max(i2),min(i3),max(i3)];
  phi{n} = phi_save{n}(box(n,1):box(n,2),box(n,3):box(n,4),box(n,5):box(n,6));
end