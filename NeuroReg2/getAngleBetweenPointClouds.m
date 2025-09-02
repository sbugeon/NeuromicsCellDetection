function angle_degrees = getAngleBetweenPointClouds(point_cloud1, point_cloud2)
%GETANGLEBETWEENPOINTCLOUDS Calculates the acute angle between two 3D point clouds.
%
%   This function determines the best-fit line for each point cloud using
%   Principal Component Analysis (PCA) and calculates the acute angle
%   between these two lines.
%
%   Args:
%       point_cloud1: An N x 3 matrix of points [x, y, z].
%       point_cloud2: An M x 3 matrix of points [x, y, z].
%
%   Returns:
%       angle_degrees: The acute angle (from 0 to 90 degrees) between the
%                      best-fit lines of the two clouds.

    % 1. Find the first principal component (direction vector) for each cloud.
    % The pca function automatically centers the data. The first column of
    % the 'coeff' matrix is the direction of highest variance (the line's vector).
    [coeff1, ~] = pca(point_cloud1);
    dir_vec1 = coeff1(:, 1);

    [coeff2, ~] = pca(point_cloud2);
    dir_vec2 = coeff2(:, 1);

    % 2. Calculate the acute angle using the dot product of the unit vectors.
    % The absolute value of the dot product gives the cosine of the acute angle.
    dot_product = abs(dir_vec1' * dir_vec2);

    % 3. Clamp value to handle potential floating-point inaccuracies near 1.0.
    dot_product = min(1.0, dot_product);

    % 4. Calculate the angle in radians and convert to degrees.
    angle_degrees = rad2deg(acos(dot_product));

end