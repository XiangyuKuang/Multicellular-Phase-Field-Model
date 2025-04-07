% first change current folder to the folder storing simulation data
load('../../CellParameters.mat','CellParameters')
AutoImg('dx',0.25,'ImgType',3,'CellPara',CellParameters,'threshold',0.25,'save','close','fig')