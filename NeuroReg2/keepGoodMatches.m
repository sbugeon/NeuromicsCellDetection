function [TransTable,GoodMatch] = keepGoodMatches(TransTable,X,Y,DataSets,pt_list_vol,AngleTolerance,DistTolerance,AngleCorrection,output_points)
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
    if ~isempty(output_points) % using detected zstack brain surface
    M1 = [R,t]; % Volume to Slice
    surf_stack =  M1*[output_points,ones(size(output_points(:,1)))]';
    gg = abs(surf_stack(2,:))<20;
    surf_stack = surf_stack(:,gg);
    b_plane=[];
    else % using zstack boundaries
        b_plane = neuroReg.cutVolumeBorder(dataZ,data_slice,M);
        surf_stack =[];
    end
    if isempty(b_plane) & isempty(surf_stack) % case where stack surface is not in match 
        b_plane = neuroReg.cutVolumeBorder(dataZ,data_slice,M);
    end
    GoodMatch(i) = validateMatch(b_plane,X,Y,AngleTolerance,-TransParameters(i,1),DistTolerance,AngleCorrection,surf_stack);
end
toc
fprintf('\n Done')
TransTable.Intensity(logical(GoodMatch)) = TransTable.Intensity(logical(GoodMatch)) + 10;
% TransTable = [TransTable(logical(GoodMatch),:);TransTable(~logical(GoodMatch),:)];
[~,S] = sort(TransTable.Intensity,'descend');
TransTable = TransTable(S,:);
%%

%     [data_out,b_plane,~] = neuroReg.cutVolume(dataZ,data_slice,M);
%     x1 = min(b_plane(1,:));
%     x2 = max(b_plane(1,:));
%     y1 = min(b_plane(2,:));
%     y2 = max(b_plane(2,:));
%     [~,ix1] = min(abs(data_slice.x - x1));
%     [~,ix2] = min(abs(data_slice.x - x2));
%     [~,iy1] = min(abs(data_slice.y - y1));
%     [~,iy2] = min(abs(data_slice.y - y2));
%     data_cut_now = data_slice;
%     data_cut_now.value = 0*data_cut_now.value;
%     data_cut_now.value(ix1:ix2,iy1:iy2) = data_out.value;
%     
%     gg = abs(surf_stack(2,:))<20;
%     
%     figure()
%     clf
%     imshow(data_cut_now.value,[])
%     hold on
%     scatter(surf_stack(3,gg),surf_stack(1,gg))
%     pause