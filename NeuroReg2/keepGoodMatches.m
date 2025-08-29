function [TransTable,GoodMatch] = keepGoodMatches(TransTable,X,Y,DataSets,pt_list_vol,AngleTolerance,DistTolerance)
dataZ = DataSets.dataZ;
data_slice = DataSets.data_slice;

TransParameters = table2array(TransTable(:,2:end));
if nargin < 6
    AngleTolerance = 40;
    DistTolerance = 100;
end

GoodMatch = zeros(size(TransParameters,1),1);
fprintf('\n finding good matches using brain surface')
tic
for i = 1:size(TransParameters,1)%54
    t = TransParameters(i,4:6)';
     [~,R,~,~,~] = ...
        neuroReg.rotateCells(pt_list_vol,...
        TransParameters(i,1),TransParameters(i,2),TransParameters(i,3));
    M = [R',-R'*t]; % M: slice to volume. Default.
%     [~,b_plane,~] = neuroReg.cutVolume(dataZ,data_slice,M);
     b_plane = neuroReg.cutVolumeBorder(dataZ,data_slice,M);
    
    GoodMatch(i) = validateMatch(mean(b_plane(1:2)),min(b_plane(3:4)),X,Y,AngleTolerance,-TransParameters(i,1),DistTolerance);
end
toc
fprintf('\n Done')
TransTable.Intensity(logical(GoodMatch)) = TransTable.Intensity(logical(GoodMatch)) + 10;
% TransTable = [TransTable(logical(GoodMatch),:);TransTable(~logical(GoodMatch),:)];
[~,S] = sort(TransTable.Intensity,'descend');
TransTable = TransTable(S,:);
end
