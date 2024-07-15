function [v_temp1,v_temp2,pt_list_vol,pt_list_slice] = getOverlap(h,TransTable,DataSets,pt_list_vol0,pt_list_slice0,Option)
% plotTransform visualize the registration from dataZ to data_slice.
% h specifies the Figure handle to plot in.
% TransTable is a 1-by-7 table. It can be one row from the output of
% rotationCorr3. It records:
% Intensity, Angle, Translation
% [{intensity},{alpha,beta,gamma},{tx,ty,yz}]
% dataZ is the volume data.
% data_slice is the particular slice.

TransTable = table2array(TransTable);
TransParameters = TransTable(1,2:end);
Integ = Option.Integ;

Subsampling = 0.4;
dataZ = DataSets.dataZ;
dataZ = neuroReg.subsample_data(dataZ, Subsampling);
data_slice = DataSets.data_slice;
data_slice.value = double(data_slice.value);
data_slice=neuroReg.subsample_data(data_slice, Subsampling);

%%
[~,R,~,~,~] = ...
    neuroReg.rotateCells(pt_list_vol0,...
    TransParameters(1),TransParameters(2),TransParameters(3));
t = TransParameters(4:6)';
M = [R',-R'*t]; % M: slice to volume. Default.
    M1 = [R,t]; % Volume to Slice

%% Plot the slice from volume
[data_cut1,b_plane,~] = neuroReg.cutVolume(dataZ,data_slice,M,Integ);
%%
x1 = min(b_plane(1,:));
x2 = max(b_plane(1,:));
y1 = min(b_plane(2,:));
y2 = max(b_plane(2,:));
[~,ix1] = min(abs(data_slice.x - x1));
[~,ix2] = min(abs(data_slice.x - x2));
[~,iy1] = min(abs(data_slice.y - y1));
[~,iy2] = min(abs(data_slice.y - y2));
data_slice_now.x = data_slice.x(ix1:ix2);
data_slice_now.y = data_slice.y(iy1:iy2);
data_slice_now.value = data_slice.value(ix1:ix2,iy1:iy2);
%%
pt_list_slice(:,1) = pt_list_slice0(:,1) - ix1;
pt_list_slice(:,2) = pt_list_slice0(:,2) - iy1;

pt_list_vol = M1*[pt_list_vol0;ones(size(pt_list_vol0(1,:)))];
% pf = abs(pt_list_vol_rotated(2,:))<Integ;
% pt_list_vol = pt_list_vol_rotated(:,pf);

pt_list_vol(1,:) = pt_list_vol(1,:) - ix1;
pt_list_vol(3,:) = pt_list_vol(3,:) - iy1;
%%

v_temp1 = data_cut1.value;
v_temp1 = intensity_normalize(v_temp1);
v_temp1 = flipud(v_temp1');

v_temp2 = data_slice_now.value;
v_temp2 = intensity_normalize(v_temp2);
v_temp2 = flipud(v_temp2');

end
function v_out = intensity_normalize(v_in)
v_in = v_in - min(v_in(:));
v_in(isnan(v_in))=0;
v_out = v_in/max(v_in(:));
v_out = double(v_out);
end
