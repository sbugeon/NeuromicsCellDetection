function GoodMatch = validateMatch(b_plane,X,Y,AngleTol,Alpha,DistTol,AngleCorrection,surf_stack)


[X,ia] = sort(X);
Y = Y(ia);
% Y = smooth(Y);
xq = [X(1):100:X(end)];
vq = interp1(X,Y,xq);

if isempty(surf_stack) 
[pgx,pgy] = minboundparallelogram(b_plane(:,1),b_plane(:,2));
pgxU = pgx(1:4);
pgyU = pgy(1:4);
[~,S] = sort(pgyU);
x0 = mean(pgx(S(1:2)));
y0 =  mean(pgy(S(1:2)));

% find 2 closest points
% tic
D = pdist([x0,xq; y0,vq]');Z = squareform(D);
[S,Use] = sort(Z(2:end,1));Use = Use(1:2);

% tic
% find angle of brain surface using 2 cloest points
opp = abs(diff(vq(Use)));
adj = abs(diff(xq(Use)));
A = rad2deg(atan(opp/adj));

GoodMatch = abs(A + Alpha - AngleCorrection)<AngleTol ;
else
    
   
    % recover angle and distance of stack surface 
    pt1 = surf_stack([1,3],:)';
    % keep only center of stack
    M = mean(pt1); M = M(1);
    pt1 = pt1(abs(pt1(:,1)-M)<50,:);
    
    D = pdist([pt1(:,1)',xq; pt1(:,2)',vq]');Z = squareform(D);
    [~,Use] = sort(mean(Z(size(pt1,1)+1:end,:),2));Use = Use(1:2);
    pt2 = [xq(Use);vq(Use)]';
    Angle = getAngleBetweenPointClouds(pt1, pt2);
    S = getMeanDistanceBetweenLines([pt1,ones(size(pt1,1),1)], [pt2,ones(size(pt2,1),1)]);
    
    GoodMatch = abs(Angle)<AngleTol ;
end

GoodMatch  =  GoodMatch & S(1)<DistTol;


%%
% % if GoodMatch
% figure(2)
% clf
% scatter(X,Y)
% hold on
% 
% scatter(xq,vq ,'+r')
% scatter(xq(Use),vq(Use) ,'+k')
% % scatter(x0,y0,'^r')
% scatter(surf_stack(1,:),surf_stack(3,:),'sk')
% scatter(pt1(:,1),pt1(:,2),'sg')
% % plot([x0 xq(Use(1))],[y0 vq(Use(1))])
% % plot([x0 xq(Use(2))],[y0 vq(Use(2))])
% plot([-1000 ; X ; max(X)+1000; -1000],[ -1000 ; Y ; -1000; -1000],'r')
% title([ num2str(GoodMatch), ' ',num2str(Angle), ' ',num2str(S(1))])
% axis equal
% pause
% end