function plotTransform(h,TransTable,DataSets,pt_list_vol,pt_list_slice,Option,pathway,M_icp_s2v,M_icp_v2s)
% plotTransform visualize the registration from dataZ to data_slice.
% h specifies the Figure handle to plot in.
% TransTable is a 1-by-7 table. It can be one row from the output of
% rotationCorr3. It records:
% Intensity, Angle, Translation
% [{intensity},{alpha,beta,gamma},{tx,ty,yz}]
% dataZ is the volume data.
% data_slice is the particular slice.

TransTable = table2array(TransTable);
Corr = TransTable(1,1);
TransParameters = TransTable(1,2:end);
StepX = Option.StepX;
CellRadius = Option.CellRadius;
Integ = Option.Integ;
dataZ = DataSets.dataZ;
data_slice = DataSets.data_slice;
data_slice_bw_low = DataSets.data_slice_bw_low;
icp_flag = 0;
if nargin == 10
    icp_flag = 1;
end
%% Visualization for transformation

if icp_flag==1 % Plot icp result
    M = M_icp_s2v; % M: slice to volume. Default.
    M1 = M_icp_v2s; % Volume to Slice
else
    [~,R,~,~,~] = ...
        neuroReg.rotateCells(pt_list_vol,...
        TransParameters(1),TransParameters(2),TransParameters(3));
    t = TransParameters(4:6)';
    M = [R',-R'*t]; % M: slice to volume. Default.
    M1 = [R,t]; % Volume to Slice
end
figure(h);
clf(h);
%% Plot the slice from volume
subplot(2,4,3);
[data_cut1,b_plane,~] = neuroReg.cutVolume(dataZ,data_slice,M,Integ);
neuroReg.plotData2(data_cut1);
title(['Cut from Volume, Thickness =',num2str(Integ),'um']);
%% Plot detected cells
figure(h);
subplot(2,4,8);
xs_full = pt_list_slice(1,:);
ys_full = pt_list_slice(2,:);
xb = b_plane(1,:);
yb = b_plane(2,:);
x_lim = [xb(1),xb(3)];
y_lim = [yb(1),yb(2)];
pt_list_slice_now = pt_list_slice(:,inpolygon(xs_full,ys_full,xb,yb));
pt_list_vol_rotated = M1*[pt_list_vol;ones(size(pt_list_vol(1,:)))];
pf = abs(pt_list_vol_rotated(2,:))<Integ;
pt_list_vol_now = pt_list_vol_rotated(:,pf);
xs = pt_list_slice_now(1,:);
ys = pt_list_slice_now(2,:);
xv = pt_list_vol_now(1,:);
yv = pt_list_vol_now(3,:);
dv = pt_list_vol_now(2,:);
pt_now_area = ones(1,sum(pf==1)).*exp(-(dv/CellRadius).^2/2);
hold on;
scatter(xs,ys,72,'k+');
scatter(xv,yv,72*pt_now_area,'ro');
axis equal;
xlim([b_plane(1,1),b_plane(1,3)]);
ylim([b_plane(2,1),b_plane(2,2)]);
hold off;
title('Detected Cells. + = Sclice, o = Volume');

%% Plot cut from reconstruction
figure(h);
subplot(2,4,7);
data_cut_r = neuroReg.renderCell2(x_lim,y_lim,Option.StepX,[xv;yv],pt_now_area,Option);
nf = isnan(data_cut_r.value);
data_cut_r.value = data_cut_r.value - mean(data_cut_r.value(~nf))*Option.MagicNumber;
data_cut_r.value(nf) = 0;
neuroReg.plotData2(data_cut_r);
title('Reconstruction')

%% Plot the geometry
figure(h);
subplot(2,4,1);
[xx0,zz0] = ndgrid(data_slice.x,data_slice.y);
yy0 = zeros(size(xx0));
surf(xx0,yy0,zz0,double(data_slice.value));
shading interp;
axis equal;
axis vis3d;
hold on;
vts = getCube([dataZ.x(1),dataZ.y(1),dataZ.z(1)],...
    [dataZ.x(end)-dataZ.x(1),dataZ.y(end)-dataZ.y(1),dataZ.z(end)-dataZ.z(1)]);
vts1 = M1*[reshape(vts,[3,24]);ones(1,24)];
vts1 = reshape(vts1,[3 4 6]);
plot3(b_plane(1,:),zeros(size(b_plane(1,:))),b_plane(2,:),'r');
drawCube(vts1,'r');
hold off;
title('Geometry');
%% Plot the slice from the ex-vivo slice
figure(h);
subplot(2,4,2);
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
neuroReg.plotData2(data_slice_now);
title('Slice');

%% Plot the correlation function XCC
figure(h);
subplot(2,4,5);
Option2 = Option;
Option2.Res0 = 0.003;
data_cut_now = data_cut1;
nf = isnan(data_cut_now.value);
data_cut_now.value = data_cut_now.value - mean(data_cut_now.value(~nf));
data_cut_now.value(nf) = 0;
data_cut_now = neuroReg.downSample(data_cut_now,StepX,[],StepX);
data_corr_now = data_slice_bw_low;
data_corr_now.value = filter2(data_cut_now.value,data_slice_bw_low.value);
neuroReg.plotData2(data_corr_now);
hold on;
plot(b_plane(1,:),b_plane(2,:),'r');
hold off;
title(['XCC, max: ',num2str(Corr)]);

%% Plot XCC from reconstruction
figure(h);
subplot(2,4,6);
data_corr_now_r = data_corr_now;
data_corr_now_r.value = filter2(data_cut_r.value,data_slice_bw_low.value);
neuroReg.plotData2(data_corr_now_r);
hold on;
plot(b_plane(1,:),b_plane(2,:),'r');
hold off;
title('XCC from reconstruction')

%% Overlap plot
figure(h);
h_temp=subplot(2,4,4);

v_temp1 = data_cut1.value;
v_temp1 = intensity_normalize(v_temp1);
v_temp1 = flipud(v_temp1');
v_temp2 = data_slice_now.value;
v_temp2 = intensity_normalize(v_temp2);
v_temp2 = flipud(v_temp2');

IM_V(:,:,1) = v_temp2;
IM_V(:,:,2) = v_temp1;
IM_V(:,:,3) = zeros(size(v_temp2,1),size(v_temp2,2));
imshow(IM_V,'Parent',h_temp);

title('Overlap');
%% Save subplots as individual tiff
if ~isempty(pathway)
    for i = 1:8
        figure(h);
        hax = subplot(2,4,i);
        hfig = figure(i);
        hax_new = copyobj(hax, hfig);
        set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
        saveas(hfig, fullfile(pathway,['MatchFound_', num2str(i) , '.tif']))
    end
end
end
function v_out = intensity_normalize(v_in)
v_in = v_in - min(v_in(:));
v_in(isnan(v_in))=0;
% v_out(v_out<0)=0;
% v_out = v_out/max(v_out(:));
v_out = v_in/max(v_in(:));
end
function vts = getCube ( origin, size )
x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1)+size(1)/2;
y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2)+size(2)/2;
z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3)+size(3)/2;
vts(1,:,:) = x;
vts(2,:,:) = y;
vts(3,:,:) = z;
end

function drawCube(vts,c)
for i=1:6
    hold on
    plot3(vts(1,:,i),vts(2,:,i),vts(3,:,i),c);
    hold off
end
h=patch(vts(1,:,6),vts(2,:,6),vts(3,:,6),'y');
alpha(h,0.3);
end


