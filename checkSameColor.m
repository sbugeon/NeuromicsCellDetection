function Same = checkSameColor(ZstackN,SliceN,Pairs)
Same=[];
for i=1:size(Pairs,1)
    Same(i) = contains(ZstackN,Pairs(i,1)) & contains(SliceN,Pairs(i,2));
end
Same = any(Same);