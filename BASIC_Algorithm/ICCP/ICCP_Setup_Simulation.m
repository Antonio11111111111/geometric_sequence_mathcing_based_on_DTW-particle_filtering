function [X_map, Y_map, Z_map, P_truth, P_ins, M_measured] = ICCP_Setup_Simulation()
% Name: ICCP_Setup_Simulation
% Description: Generates a simulated geomagnetic map, ground truth
%              trajectory, and erroneous INS track with measurements.
% Parameters:
%   None
% Output:
%   X_map: (matrix) X-coordinates of the map grid.
%   Y_map: (matrix) Y-coordinates of the map grid.
%   Z_map: (matrix) Z-values (magnetic intensity) of the map grid.
%   P_truth: (N x 2) Ground truth trajectory points.
%   P_ins: (N x 2) Erroneous INS-indicated trajectory points.
%   M_measured: (N x 1) Geomagnetic measurements along the truth track.
% Example:
%   [X_map, Y_map, Z_map, P_truth, P_ins, M_measured] = ICCP_Setup_Simulation();

% 1. Generate geomagnetic map (using peaks function to simulate complex field)
[X_map, Y_map] = meshgrid(linspace(-25, 25, 200), linspace(-25, 25, 200));
Z_map = peaks(X_map, Y_map) * 100 + 50000; % Simulate field around 50000 nT

% 2. Generate ground truth trajectory
num_points = 100; % Number of track points
t_path = linspace(0, 1.5*pi, num_points);
P_truth = [-1 + 1.5*cos(t_path); 0.5 + 1.5*sin(t_path)]';
P_truth = P_truth + [linspace(0, 0.5, num_points)', linspace(0, -0.5, num_points)'];

% 3. Get true geomagnetic measurements (with slight noise)
% Sample on the geomagnetic map using interp2
M_measured = interp2(X_map, Y_map, Z_map, P_truth(:,1), P_truth(:,2));
M_measured = M_measured + randn(num_points, 1) * 0.5; % Measurement noise

% 4. Generate INS indicated track with errors
% Apply a known affine transform (scale + rotation + translation)
theta = 15 * (pi/180); % Rotate 15 degrees
R_truth = [cos(theta), -sin(theta); sin(theta), cos(theta)];
s_truth = 1.1;  % Scale factor 110% (affine error)
t_truth = [0.3; -0.4]; % Translation
ins_noise = 0.02; % INS track noise

% Apply transform: P_ins = s * R * P_truth + t
P_ins = (s_truth * R_truth * P_truth' + t_truth)';
P_ins = P_ins + randn(size(P_ins)) * ins_noise;

end
