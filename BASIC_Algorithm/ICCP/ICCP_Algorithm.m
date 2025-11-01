function [P_corrected, R_final, t_final, errors, k] = ICCP_Algorithm(P_ins, M_measured, X_map, Y_map, Z_map)
% Name: ICCP_Algorithm
% Description: Implements the classic Iterative Closest Contour Point (ICCP) 
%              algorithm for a 2D rigid transformation (rotation + translation).
% Parameters:
%   P_ins:      (N x 2) INS indicated track points (source points).
%   M_measured: (N x 1) Corresponding geomagnetic measurements.
%   X_map:      (matrix) X-coordinates of the map grid.
%   Y_map:      (matrix) Y-coordinates of the map grid.
%   Z_map:      (matrix) Z-values (magnetic intensity) of the map grid.
% Output:
%   P_corrected: (N x 2) Matched/corrected track.
%   R_final:     (2 x 2) Final accumulated rotation matrix.
%   t_final:     (2 x 1) Final accumulated translation vector.
%   errors:      (k x 1) RMSE history for each iteration.
%   k:           (scalar) Number of iterations performed.
% Example:
%   [P_corr, R, t, err, iters] = ICCP_Algorithm(P_ins, M_m, X_m, Y_m, Z_m);

% Iteration parameters
max_iter = 20;      % Max number of iterations
tol = 1e-6;         % Convergence tolerance
P_current = P_ins;  % Track at current iteration
errors = zeros(max_iter, 1);

% Accumulated transform parameters
R_total = eye(2);
t_total = zeros(2, 1);

disp('Iteration started:');
fprintf('Iter |      RMSE\n');
fprintf('---------------------\n');

% Loop until convergence or max iterations
for k = 1:max_iter
    
    % Step 1: Find closest contour points (Matching)
    % Y is the set of target points on the map contours
    Y = ICCP_Find_Closest_Contour_Points(P_current, M_measured, X_map, Y_map, Z_map);
    
    % Step 2: Solve for optimal rigid transform (Alignment)
    % Solves for R, t such that Y = R * P_current + t
    [R, t] = ICCP_Find_Rigid_Transform(P_current, Y);
    
    % Step 3: Apply transformation
    % Update the track points: P_new = R*P_current + t
    P_new = (R * P_current' + t)';
    
    % Accumulate transformation
    R_total = R * R_total;
    t_total = R * t_total + t;
    
    % Step 4: Check for convergence
    % Calculate RMSE between P_new and P_current
    diff = P_new - P_current;
    rmse = sqrt(mean(sum(diff.^2, 2)));
    errors(k) = rmse;
    
    fprintf('%4d | %13.8f\n', k, rmse);
    
    P_current = P_new;
    
    if rmse < tol
        disp('Algorithm converged.');
        break; % Convergence tolerance met
    end
end

if k == max_iter
    disp('Reached max iterations.');
end

errors = errors(1:k); % Trim error log to effective iterations
P_corrected = P_current;
R_final = R_total;
t_final = t_total;

end
