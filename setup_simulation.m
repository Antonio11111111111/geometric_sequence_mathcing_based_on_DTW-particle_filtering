% =========================================================================
% File 2:  setup_simulation.m
% =========================================================================

% setup_simulation.m
%
% Generate simulation data:
% 1. Geomagnetic reference map (Z_map)
% 2. Ground truth trajectory (P_truth)
% 3. Geomagnetic measurements corresponding to the true track (M_measured)
% 4. INS indicated track with affine transform + noise (P_ins)

function [X_map, Y_map, Z_map, P_truth, P_ins, M_measured] = setup_simulation()

% 1. Generate geomagnetic map (using peaks function to simulate complex field)
[X_map, Y_map] = meshgrid(linspace(-3, 3, 200), linspace(-3, 3, 200));
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