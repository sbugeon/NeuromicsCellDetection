function [I,R] = data2img(data)
% This function deals with the image storage convention.
% imshow(I) will result in the same image as pcolor(data.value')
% R specifies picture position in the real world.

I = flipud(data.value');
R = imref2d(size(I),[data.x(1),data.x(end)],-[data.y(end),data.y(1)]);

end

