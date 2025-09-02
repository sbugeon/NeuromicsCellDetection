function mean_dist = getMeanDistanceBetweenLines(point_cloud1, point_cloud2)
%GETMEANDISTANCESHORTERTOLONGER Calculates a one-way mean distance between two lines.
%
%   (Corrected Version)
%   This function identifies the shorter point cloud (by number of points)
%   and calculates the mean distance from its points to the best-fit line
%   of the longer point cloud.
%
%   Args:
%       point_cloud1: An N x 3 matrix of points [x, y, z].
%       point_cloud2: An M x 3 matrix of points [x, y, z].
%
%   Returns:
%       mean_dist: A single scalar value representing the one-way mean distance.

    % 1. Determine which point cloud is shorter and which is longer.
    if size(point_cloud1, 1) < size(point_cloud2, 1)
        shorter_cloud = point_cloud1;
        longer_cloud = point_cloud2;
    else
        shorter_cloud = point_cloud2;
        longer_cloud = point_cloud1;
    end

    % 2. Get the properties of the LONGER line.
    [dir_vec_longer, centroid_longer] = get_line_properties(longer_cloud);

    % 3. Calculate vectors from each point in the SHORTER cloud to the
    %    centroid of the LONGER line.
    vecs_to_line = shorter_cloud - centroid_longer;

    % 4. *** FIX IS HERE ***
    % Replicate the direction vector to match the size of vecs_to_line
    % before taking the row-wise cross product.
    dir_vec_mat = repmat(dir_vec_longer', size(shorter_cloud, 1), 1);
    
    % Calculate the norm of the cross product for each point to get distances.
    dists = vecnorm(cross(vecs_to_line, dir_vec_mat, 2), 2, 2);

    % 5. The result is the mean of these calculated distances.
    mean_dist = mean(dists);

end

function [direction_vector, centroid] = get_line_properties(points)
    % Helper function to find the line's direction and a point on the line.
    centroid = mean(points, 1);
    [coeff, ~] = pca(points);
    direction_vector = coeff(:, 1); % PCA returns a column vector
end