function I = findCutStack(dataZ_mid,data_slice,M1)

% from slice, stack and transformation, finds the boundary of the
% intersection of the two datasets

% depends upon Matt J (2020). Analyze N-dimensional Polyhedra in terms of 
% Vertices or (In)Equalities 
%(https://www.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-polyhedra-in-terms-of-vertices-or-in-equalities), MATLAB Central File Exchange. Retrieved January 10, 2020.

data = dataZ_mid;
step = 20;
data.x = data.x(1:step:end);
data.y = data.y(1:step:end);
data.z = data.z(1:step:end);
[X,Y,Z] = meshgrid(data.x,data.y,data.z);


c=cat(3,X,Y,Z);
d=reshape(c,[],3);
d_rotated = M1*[d,ones(size(d(:,1)))]';
d_rotated = d_rotated';
d_rotated  = [d_rotated(:,1) d_rotated(:,3) d_rotated(:,2)];

dataS.x = data_slice.x(1:step:end);
dataS.y = data_slice.y(1:step:end);

[X,Y] = meshgrid(dataS.x,dataS.y);


cS=cat(2,X,Y);
dS=reshape(cS,[],2);
dS = [dS zeros(length(dS),1)];

 I = intersectionHull('vert',dS,'vert',d_rotated); 
 

end

