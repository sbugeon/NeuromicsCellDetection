function output_points = getbrainsurfaceZstack(data,THRESHOLD)


dd = data.value;
SMOOTH_FACTOR = 50;


% 1. Define the column indices we will process.
j_indices = 1:size(dd, 2);
% 2. Work on a subset of the data corresponding to the j_indices.
% This is more memory-efficient than operating on the whole 'dd' matrix.
dd_subset = dd(:, j_indices, :);
% 3. Find all values above the threshold in a single vectorized operation.
% This creates a 3D logical (true/false) matrix.
is_above_threshold = (dd_subset > THRESHOLD);
% 4. Find the index of the *first* occurrence of 'true' along the 3rd dimension.
% The max function, when used on a logical array, finds the first '1' (true).
% The second output argument of max is the index of the maximum value.
[~, surface_raw] = max(is_above_threshold, [], 3);
% 5. Handle the case where no value in a vector exceeded the threshold.
% In this scenario, max returns an index of 1, which is incorrect.
% We find where no threshold was crossed and set those indices to NaN.
no_match_found = ~any(is_above_threshold, 3);
surface_raw(no_match_found) = NaN;
% 6. Pre-allocate the final output matrix.
surface0 = NaN(size(surface_raw, 1), size(surface_raw, 2));
% 7. Apply smoothing to each column. A loop here is fine and efficient,
% as it operates on columns and avoids the slow, deeply nested structure.
for k = 1:size(surface_raw, 2)
    surface0(:, k) = smooth(surface_raw(:, k), SMOOTH_FACTOR);
end

% --- Plane Fitting and Grid Sampling Script ---
% --- Configuration ---
% Define the number of points you want in your final output grid.
GRID_POINTS_X = 20; % Number of points along the X-axis of the new grid
GRID_POINTS_Y = 20; % Number of points along the Y-axis of the new grid

% --- 1. Prepare Data for Plane Fitting ---
% First, create the X and Y coordinate grids corresponding to 'surface0'.
y_coords = 1:size(surface0, 1); 
x_coords = j_indices;
[X, Y] = meshgrid(x_coords, y_coords);

% The equation of a plane is z = a*x + b*y + c. We need to find a, b, and c.
% We need to solve the system of linear equations in a least-squares sense.
% Extract the valid (non-NaN) data points into column vectors.
valid_indices = ~isnan(surface0);
x_data = X(valid_indices);
y_data = Y(valid_indices);
z_data = surface0(valid_indices);

% Construct the matrix for the linear system.
% The columns are [x, y, 1] to solve for the coefficients [a; b; c].
A = [x_data, y_data, ones(size(x_data))];

% --- 2. Solve for the Plane Coefficients ---
disp('Fitting plane to the surface data...');

% Use MATLAB's backslash operator (mldivide) for a robust and efficient
% least-squares solution. This finds the coefficients that best fit the data.
% plane_coeffs will be a 3x1 vector: [a; b; c]
plane_coeffs = A \ z_data;
a = plane_coeffs(1);
b = plane_coeffs(2);
c = plane_coeffs(3);

fprintf('Best-fit plane equation: Z = (%.4f) * X + (%.4f) * Y + (%.4f)\n', a, b, c);

% --- 3. Generate the Sample Grid on the Fitted Plane ---

disp('Generating sample points on the fitted plane...');

% Create a new, evenly spaced grid of X and Y coordinates.
% This grid will cover the same range as the original data.
x_plane_vector = linspace(min(x_data), max(x_data), GRID_POINTS_X);
y_plane_vector = linspace(min(y_data), max(y_data), GRID_POINTS_Y);

% Use meshgrid to create the sampling grid.
[X_plane, Y_plane] = meshgrid(x_plane_vector, y_plane_vector);

% Calculate the Z values for every point on this new grid using the plane equation.
Z_plane = a * X_plane + b * Y_plane + c;

% --- 4. Format the Final Output ---

% Reshape the X, Y, and Z grid matrices into a single N x 3 matrix of 3D points.
% Each row represents a point [x, y, z].
output_points = [X_plane(:), Y_plane(:), Z_plane(:)];

fprintf('Generated %d sample points forming a grid on the plane.\n', size(output_points, 1));
disp('The variable "output_points" contains the [X, Y, Z] coordinates.');
disp('Preview of the first 5 points:');
disp(output_points(1:5, :));

% --- 5. (Optional) Visualize the Result ---
disp('Creating visualization...');
figure;
hold on; % Allow multiple plots in the same figure

% Plot the original surface data as a mesh
h_orig = surf(X, Y, surface0);
h_orig.FaceAlpha = 0.7; % Make it slightly transparent
h_orig.EdgeColor = 'none'; % Hide the mesh lines for clarity
title('Original Surface and Fitted Plane');

% Plot the new, fitted plane on top of the original data
h_plane = surf(X_plane, Y_plane, Z_plane);
h_plane.FaceColor = 'r'; % Make the plane red
h_plane.EdgeColor = [0.7 0 0]; % Darker red edges

% Plot the sampled grid points as black dots on the plane
% We add a small offset in Z to make them visible on top of the plane surface
plot3(output_points(:,1), output_points(:,2), output_points(:,3) + 0.1, 'k.', 'MarkerSize', 15);

xlabel('X-axis (Column Index)');
ylabel('Y-axis (Row Index)');
zlabel('Z-axis (K Index)');
legend('Original Surface', 'Fitted Plane', 'Sampled Grid Points');
grid on;
axis tight;
view(3); % Set to 3D view
hold off;
