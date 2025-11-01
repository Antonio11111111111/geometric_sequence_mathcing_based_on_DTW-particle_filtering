% =========================================================================
% File 1:  Main_ICCP_Geomagnetic.m
% (This is the main file you need to run)
% =========================================================================
% Main_ICCP_Geomagnetic.m
%
% Demonstration of the classic ICCP algorithm for geomagnetic matching.
% (Performs rigid transformation only: rotation and translation)
clear; clc; close all;
disp('Setting up simulation environment...');

% 1. Set simulation parameters and get data
% Note: Simulation data still includes a 1.1x scale error 
% to observe the limitations of the classic ICCP algorithm.
[X_map, Y_map, Z_map, P_truth, P_ins, M_measured] = ICCP_Setup_Simulation();

disp('Running classic ICCP algorithm...');
% 2. Run the classic ICCP matching algorithm (solving only for R and t)
[P_corrected, R_final, t_final, errors, iter_count] = ICCP_Algorithm(P_ins, M_measured, X_map, Y_map, Z_map);

disp('Algorithm finished.');
disp('Final rotation matrix (R):');
disp(R_final);
disp('Final translation vector (t):');
disp(t_final);

% 3. Visualize results
disp('Generating visualization plots...');

% === Figure 1: Trajectory Matching Result (Enhanced) ===
figure('Name', 'Classic ICCP Geomagnetic Matching Result', 'NumberTitle', 'off', 'Color', 'w');
set(gca, 'FontSize', 11, 'FontName', 'Arial');
hold on;

% Plot geomagnetic reference map (filled contours)
contourf(X_map, Y_map, Z_map, 20, 'LineStyle', 'none');
colormap(gca, 'parula');
cb = colorbar;
cb.Label.String = 'Magnetic Field (nT)';

% Plot trajectories
h1 = plot(P_truth(:,1), P_truth(:,2), 'k-', 'LineWidth', 2.5);
h2 = plot(P_ins(:,1), P_ins(:,2), 'r--', 'LineWidth', 2.5);
h3 = plot(P_corrected(:,1), P_corrected(:,2), 'g-o', ...
    'LineWidth', 2, ...
    'MarkerSize', 4, ...
    'MarkerFaceColor', 'g', ...
    'MarkerEdgeColor', 'k');

title('Classic ICCP Geomagnetic Matching (Rigid Transform)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('East Position / m');
ylabel('North Position / m');
legend([h1, h2, h3], ...
       'Ground Truth', ...
       'INS Indicated Track (Initial)', ...
       'ICCP Matched Track (Corrected)', ...
       'Location', 'northwest');
axis equal;
grid on;
box on;
hold off;

% === Figure 2: RMSE Convergence (Enhanced) ===
figure('Name', 'ICCP Algorithm Convergence', 'NumberTitle', 'off', 'Color', 'w');
set(gca, 'FontSize', 11, 'FontName', 'Arial');
hold on;

plot(1:iter_count, errors, 'b-x', 'LineWidth', 1.5, 'MarkerSize', 6);

title('Matching Error Convergence Curve', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Iteration');
ylabel('Root Mean Square Error (RMSE)');
grid on;
box on;
if iter_count > 1
    xlim([1, iter_count]);
end
legend('RMSE', 'Location', 'northeast');
hold off;

disp('Visualization complete.');
