function TransTable = keepGoodMatches(TransTable,X,Y)

TransParameters = table2array(TransTable(:,2:end));
AngleTolerance = 10;
GoodMatch = zeros(size(TransParameters,1),1);
for i = 1:size(TransParameters,1)
    t = TransParameters(i,4:6)';
    
    GoodMatch(i) = validateMatch(t(1),t(3),X,Y,AngleTolerance,-TransParameters(i,1));

end

TransTable = TransTable(logical(GoodMatch),:);
end
