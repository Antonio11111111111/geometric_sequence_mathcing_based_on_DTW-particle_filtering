function [R, t] = ICCP_Find_Rigid_Transform(P_source, P_target)
% Name: ICCP_Find_Rigid_Transform
% Description: Finds the optimal 2D rigid transformation (Rotation R,
%              Translation t) that maps a source point set (P_source) to a
%              target point set (P_target) by minimizing the squared error.
%              Uses the SVD method.
% Parameters:
%   P_source: (N x 2) Source point set.
%   P_target: (N x 2) Target point set.
% Output:
%   R: (2 x 2) Optimal rotation matrix.
%   t: (2 x 1) Optimal translation vector.
% Example:
%   [R, t] = ICCP_Find_Rigid_Transform(P_ins, Y_contours);

% 1. Calculate centroids
centroid_source = mean(P_source, 1);
centroid_target = mean(P_target, 1);

% 2. Center the point clouds
P_source_centered = P_source - centroid_source;
P_target_centered = P_target - centroid_target;

% 3. Calculate the covariance matrix H
H = P_source_centered' * P_target_centered;

% 4. Compute SVD
[U, ~, V] = svd(H);

% 5. Calculate rotation matrix R
R = V * U';

% 6. Handle reflection case (det(R) = -1)
if det(R) < 0
    %fprintf('Reflection detected, correcting...\n');
    V(:, 2) = -V(:, 2);
    R = V * U';
end

% 7. Calculate translation vector t
t = centroid_target' - R * centroid_source';

end
