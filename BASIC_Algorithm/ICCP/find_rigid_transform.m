% =========================================================================
% File 4:  find_rigid_transform.m
% =========================================================================

% find_rigid_transform.m
%
% Given two point sets, P (INS indicated) and Y (closest contour),
% calculates the optimal rigid transform (R, t) that minimizes
% d = sum || y_i - (R*p_i + t) ||^2.

function [R, t] = find_rigid_transform(P, Y)
    
    N = size(P, 1);
    % Check if point set is empty
    if N == 0
        R = eye(2); t = zeros(2,1);
        return; % If empty, return identity transform and exit
    end

    % 1. Calculate centroids of point sets
    p_bar = mean(P, 1)'; 
    y_bar = mean(Y, 1)'; 

    % 2. Calculate relative coordinates
    delta_P = P' - p_bar; 
    delta_Y = Y' - y_bar; 

    % 3. Construct the 2x2 matrix M
    M = delta_Y * delta_P'; 

    % 4. Calculate rotation matrix R (using SVD)
    [U, ~, V] = svd(M);
    R = U * V'; 
    
    % Ensure R is a pure rotation matrix (prevent reflection)
    if det(R) < 0
        V(:, end) = -V(:, end);
        R = U * V';
    end

    % 5. Calculate translation vector t
    t = y_bar - R * p_bar; 

end  % This is the single, required 'end' for the function