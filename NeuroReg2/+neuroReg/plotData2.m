function plotData2(data)
% plotSlice(data) visualize a 2-D data set
% data.value = fliplr(data.value);
imagesc(data.x,data.y,data.value');
colormap(gray)
axis equal;
axis tight;
shading flat;
set(gca,'Ydir','normal')
end

