% =========================================================================
% File 3:  iccp_algorithm.m
% =========================================================================

% iccp_algorithm.m
%
% Implements the classic ICCP algorithm (rigid transform only).
%
% Inputs:
%   P_ins:      INS indicated track points (N x 2)
%   M_measured: Corresponding geomagnetic measurements (N x 1)
%   X_map, Y_map, Z_map: Geomagnetic reference map
%
% Outputs:
%   P_corrected: Matched/corrected track (N x 2)
%   R_final:     Final rotation matrix
%   t_final:     Final translation vector
%   errors:      RMSE history for each iteration

function [P_corrected, R_final, t_final, errors] = iccp_algorithm(P_ins, M_measured, X_map, Y_map, Z_map)

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
    Y = find_closest_contour_points(P_current, M_measured, X_map, Y_map, Z_map);
    
    % Step 2: Solve for optimal rigid transform (Alignment)
    % Solve only for rotation R and translation t
    [R, t] = find_rigid_transform(P_current, Y);
    
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