function SlicePos = loadSlicePos(Path,SliceN)
% Path = 'G:\invivoReg\TestAni\SlicePos.txt';
D = readtable(Path); 
dd = table2cell(D(:,1));
gg = strcmp(dd,SliceN);

SlicePos = table2array(D(gg,2));