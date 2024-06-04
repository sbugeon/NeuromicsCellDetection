
function [ScoreStack , ScoreSlice] = make_icp(TransTable,pt_list_slice, pt_list_vol,data_slice,dataZ_mid,d)
% ICP based on the current transform matrix.
% obj.icp plots and print the icp result.
% obj.icp(d) plots and print the icp result, with the maximum
% neigbhor range d.

if nargin<5
    d = 15;
end
MaxCellMatch = 40; % maximum of cells
TransParameters = table2array(TransTable(:,2:end));

ScoreStack = NaN(size(TransParameters,1),1);
ScoreSlice = NaN(size(TransParameters,1),1);
 fprintf('\n Getting Peaks  \n')
for i = 1: size(TransParameters,1)
    %% Get M2
    [~,R,~,~,~] = ...
        neuroReg.rotateCells(pt_list_vol,...
        TransParameters(i,1),TransParameters(i,2),TransParameters(i,3));
    t = TransParameters(i,4:6)';
    M = [R',-R'*t]; % M: slice to volume. Default.
    M_23 = M(:,[1 3 4]);
    M1 = [R,t];
    
    [~,~,MyNeighb,~] = neuroReg.PointCloudRegister2(pt_list_vol',pt_list_slice',M_23',d,[]);
    try
    I = neuroReg.findCutStack(dataZ_mid,data_slice,M1);
    catch
        ScoreStack(i) = 0;
        ScoreSlice(i) = 0;
        continue;
    end
    pt_list_vol_rotated = M1*[pt_list_vol;ones(size(pt_list_vol(1,:)))];
   
    if ~isempty(I.vert)
        xv = I.vert(:,1);
        yv = I.vert(:,2);
        k = boundary(xv,yv,0.001);
        xq  = pt_list_slice(1,:)';
        yq  = pt_list_slice(2,:)';
        
        pt_slice_in = inpolygon(xq,yq,xv(k),yv(k));
        pt_vol_in = pt_list_vol_rotated(2,:)>= -d &   pt_list_vol_rotated(2,:)<= d;
        
        % Score for Cells in Zstack matched
        ScoreStack(i) = (length(unique(MyNeighb))/sum(pt_vol_in)) + (length(unique(MyNeighb))^2/MaxCellMatch^2);
        % Score for Cells in Slice matched
        ScoreSlice(i) =  length(unique(MyNeighb))/sum(pt_slice_in)+ (length(unique(MyNeighb))^2/MaxCellMatch^2);
        
       
    else
        ScoreStack(i) = 0;
        ScoreSlice(i) = 0;
    end
    
end
