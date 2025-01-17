function [M, Error, MyNeighb, UseMe] = PointCloudRegister2(y0, x0, M0, DistScale, Options)
% [M, Error] = PointCloudRegister(y, x, M0, DistScale, Options)
% 
% Perform point cloud registration to map points x onto points y by
% iterative closest point: repeatedly finding the best y for each x, 
% and doing linear regression to find the M that maps best maps x to y
%
% inputs:
% M0 is initial transformation matrix, default identity
% DistScale: any x whose nearest neighbor is further than this won't count,
% default inf
%
% Options: what type of fit. For now ignored, the only option is a general linear
% model where x gets an extra column of ones and M is 2x3.

MaxIter = 100;
Interactive = 0;

[nP, nD] = size(x0);
x = [x0, ones(nP, 1)];

if nargin<3 || isempty(M0)
    M0 = [eye(nD) ; zeros(1, nD)];
elseif min(size(M0))==1
    % vector input means that is a shift
    M0 = [eye(nD) ; M0(:)'];
end
M = M0;

if nargin<4 || isempty(DistScale)
    DistScale = inf;
end

% make kd tree - default options!
k0 = KDTreeSearcher(y0);

% find well isolated points as those whose second neighbor is far
%[~, d2] = k0.knnsearch(y0, 'k', 2);
% if isfinite(DistScale)
%     y = y0(d2(:,2)>DistScale*2,:);
% else
%     y=y0;
% end
y=y0;

k = KDTreeSearcher(y);

Neighbor = zeros(nP,1);
for i=1:MaxIter
    LastNeighbor = Neighbor;
    xM = x*M;
    [Neighbor, Dist] = k.knnsearch(xM);
    UseMe = Dist<DistScale;
    MyNeighb = Neighbor(UseMe);
    M = x(UseMe,:)\y(MyNeighb,:);
    
    Error = sqrt(mean(Dist(UseMe).^2));
    
    if Interactive
        figure(2938764);
        fprintf('Iteration %d: %d matches, mean error %f\n', i, sum(UseMe), Error);
        clf; hold on
         plot(y(:,1), y(:,2), 'r.');
         plot(xM(:,1), xM(:,2), 'g.');
        plot([xM(UseMe,1) y(MyNeighb,1)]', [xM(UseMe,2) y(MyNeighb,2)]', 'k-');
        
        
        drawnow;
        pause
    end
    if isequal(LastNeighbor, Neighbor); break; end
end

