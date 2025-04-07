SimNum = 100;
SimName = 'vn1';
% SimName = 'sn0d2';
dr = 0.25;
th = 0.44;
contact_area = cell(SimNum,1);
mem_area = cell(SimNum,1);
for SimIdx = 1:SimNum
  load(sprintf('%1$s/%1$s_%2$d/%1$s_%2$d_44700.mat',SimName,SimIdx),'phi_save','CellIdx');
  [contact_area{SimIdx}, mem_area{SimIdx}] = ContactArea(phi_save,dr,th);
end
save(sprintf('%1$s/%1$s_contact_area.mat',SimName),"mem_area","contact_area")