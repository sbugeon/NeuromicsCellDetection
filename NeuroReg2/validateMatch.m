function GoodMatch = validateMatch(x0,y0,X,Y,AngleTol,Alpha)



% interpolate drawn boundaries
% tic
xq = [X(1):100:X(end)];
vq = interp1(X,Y,xq);
% toc


% find 2 closest points
% tic
D = pdist([x0,xq; y0,vq]');Z = squareform(D);
[S,Use] = sort(Z(2:end,1));Use = Use(1:2);
% toc

% tic
% find angle of brain surface using 2 cloest points
opp = abs(diff(vq(Use)));
adj = abs(diff(xq(Use)));
A = rad2deg(atan(opp/adj));

GoodMatch = abs(A + Alpha)<AngleTol;

MatchIn = inpolygon(x0,y0,[-1000 ; X ; max(X)+1000],[ -1000 ; Y ; -1000]) ;

GoodMatch  = GoodMatch & ~MatchIn & S(1)<100;
% toc
%%
% figure(1)
% clf
% scatter(X,Y)
% hold on
% scatter(x0,y0,'or')
% scatter(xq,vq ,'+k')
% plot([x0 xq(Use(1))],[y0 vq(Use(1))])
% plot([x0 xq(Use(2))],[y0 vq(Use(2))])
% plot([-1000 ; X ; max(X)+1000; -1000],[ -1000 ; Y ; -1000; -1000],'r')
% title(num2str(GoodMatch))
%     pause