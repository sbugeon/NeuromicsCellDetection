function [pt_list, pt_area] = gaussfiltdetectCells3(res, Res0, data, SizeLimitPx, stepX, stepY, stepZ,Option)

tic
fprintf('Finalizing result...\n');
res2 = imhmin(res, Res0*0.4);
im_bw = single(res2>Res0);
 im_bw = imclearborder(im_bw); % connectivity parameter???
CC = bwconncomp(im_bw,18);
stats = regionprops(CC,'Area','Centroid','Image','BoundingBox');
pt_list = [nan nan nan]';
pt_area = [];
k=1;
for i = 1:length(stats)
    if stats(i).Area>SizeLimitPx(1)
        if stats(i).Area<SizeLimitPx(2)
            pt_list(:,k) = stats(i).Centroid';
            pt_area(k) = stats(i).Area;
            k=k+1;
        else
            bw = stats(i).Image;
            [pt_this,pt_area_this] = neuroReg.splitCells(bw,stepX,stats(i).Area,0);
            n = length(pt_area_this);
            if n>=1
                for j = 1:n
                    pt_list(:,k) = pt_this(:,j) + stats(i).BoundingBox(1:3)';
                    pt_area(k) = pt_area_this(j);
                    k = k+1;
                end
            end
        end
        
    end
end
pt_area = pt_area*(stepX*stepY*stepZ);
if strcmp(Option.Detect3Mode,'Green')
    % If it Green mode, the area will not be accurate.
    % To fix this, a unified cell volume is assigned.
    pt_area = ones(size(pt_area))*mean(Option.SizeLimit);
end
temp = pt_list(1,:);
pt_list(1,:) = pt_list(2,:);
pt_list(2,:) = temp;
pt_list(1,:) = pt_list(1,:)*stepX + min(data.x) - stepX;
pt_list(2,:) = pt_list(2,:) * stepY + min(data.y) - stepY;
pt_list(3,:) = pt_list(3,:) * stepZ + min(data.z) - stepZ;
fprintf('Result done.\n Number of Cells = %d\n',length(pt_area));
toc
%% Visualized feedback
assignin('base','data_temp_749',data);
neuroReg.PlotSlices('data_temp_749','y',pt_list,pt_area,7,[]);

fprintf('Plot done \n')
end