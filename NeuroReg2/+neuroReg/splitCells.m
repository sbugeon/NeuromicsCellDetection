function [pt,pt_area] = splitCells(bw,step,Size_cell, ShowGraph)
min_dist = 7/step; %%%%%%%%%%%%%%%%%%%

D = bwdist(~bw);
D = -D;
D(~bw) = Inf;

mask = imextendedmin(D,0.3);
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = bw;
bw3(Ld2 == 0) = 0;
Ld2(~bw3) = 0;
stats = regionprops(Ld2,'Area','Centroid');
n = length(stats);
SizeLimitPx = 0.05*sum(bw(:))/n; %%%%%%%%%%%%%%%%%%%
pt_area = [];
pt = [nan,nan,nan]';
k=1;

for i = 1:n
    if stats(i).Area>SizeLimitPx % Leave out small parts
        pt(:,k) = stats(i).Centroid';
        pt_area(k) = stats(i).Area;
        k = k+1;
    end
end


% Visualization for testing
if ShowGraph
    figure(10098);cla;
    IS = patch(isosurface(bw,0.5), 'Clipping','off');
    axis equal, title(num2str(Size_cell))
    xlabel x, ylabel y, zlabel z,
    view(3), camlight, lighting gouraud
    % ax = gca;
    % ax.Clipping = 'off';   % get the current axis
    % IS.Clipping = 'off';    % turn clipping off
    IS.FaceColor = 'blue';
    IS.EdgeColor = 'none';
    
    figure(5),cla
    for i = 1 : 16
        h(i) = subplot(4,4,i);
    end
    for i = 1 : n
        G=(Ld2==i);
        
        h(i) = subplot(4,4,i);
        cla(h(i))
        isosurface(G)
        xlabel x, ylabel y, zlabel z
        view(3), camlight, lighting gouraud
        axis equal
    end
    
end

if k>= 3
    X=pt';
    Z = linkage(X,'ward');
    %     dendrogram(Z,'ColorThreshold',min_dist)
    c = cluster(Z, 'Cutoff', min_dist, 'criterion', 'distance');
    [idx,C] = kmeans(X,max(c));
    
    pt = C';
    pt_area_old = pt_area; pt_area = [];
    for i = 1 : size(pt,2)
        pt_area(i) = sum(pt_area_old(idx == i));
    end
    
    if ShowGraph
        figure(69);
        cla, isosurface(bw,0.5), axis equal, title('BW')
        xlabel x, ylabel y, zlabel z
        view(3), camlight, lighting gouraud, alpha(0.2)
        hold on
        scatter3(C(:,1),C(:,2),C(:,3),'b' )
        hold on
        scatter3(X(:,1),X(:,2),X(:,3), '.r' )
    end
end

if ShowGraph
for i = 1 : 16
    cla(h(i));
end
end