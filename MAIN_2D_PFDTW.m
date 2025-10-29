%==========================================================================
% SCRIPT: MAIN_2D_PFDTW
%
% DESCRIPTION:
%   Main simulation script for a 2D Adaptive Particle Filter (APF)
%   aided by PDR and geomagnetic matching using DTW.
%   The particle count (M) is adapted based on the minimum DTW distance.
%
% DEPENDENCIES:
%   - Signal Processing Toolbox (for 'dtw' function)
%   - Particle_Filter_DTW_Step_2D.m
%   - Geometric_Map_Generator.m
%   - get_next_step_random.m
%==========================================================================

clc;
clear;
close all;

%% 1. SIMULATION PARAMETERS
%==========================================================================

% --- Simulation Control ---
M_init = 500;            % Initial particle count
NUM_STEPS = 100;          % Number of simulation steps
SEQUENCE_LEN = 50;        % Sequence length for DTW

% --- Map Dimensions ---
MAP_X_LEN = 50;           % Map width (X)
MAP_Y_LEN = 50;           % Map height (Y)

% --- Noise Parameters ---
SENSOR_NOISE_STD = 0.5;   % Std dev of sensor noise (added to live sequence)
DTW_NOISE_STD = 30;      % Std dev for weight calculation (converts DTW dist to weight)
process_noise.step_std = 0.5;          % Process noise: std dev for step length
process_noise.theta_std = deg2rad(10); % Process noise: std dev for heading

% --- Particle Re-injection Noise ---
pos_std = 5.0;           % Std dev for re-injected particle position
ang_std = 0.2;            % Std dev for re-injected particle angle

% --- APF (Adaptive Particle Filter) Parameters ---
% APF.M_min = 500;          % Minimum allowed particle count
% APF.M_max = 100000;      % Maximum allowed particle count
% APF.DTW_THRESH_HIGH = 22; % DTW distance above which M is increased
% APF.DTW_THRESH_LOW  = 10.0; % DTW distance below which M is decreased
KLD.n_min = 500;          % 最小粒子数
KLD.n_max = 100000;       % 最大粒子数 (安全上限)
KLD.bin_size_xy = 0.5;    % [重要] X/Y 分箱大小 (米)
KLD.epsilon = 0.05;       % [重要] KLD 误差界 (e.g., 0.05)
KLD.delta = 0.01;         % [重要] KLD 概率界 (1-0.01 = 99% 置信度)
KLD.DTW_THRESH_HIGH = 23.9; % 加回这一行

%%
sigma = 0.05;
eps = 0.5;

%% 2. MAP GENERATION
%==========================================================================
fprintf('Generating geomagnetic map...\n');

[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
Mag = imgaussfilt(Mag_raw, 3.0); % Smooth the map

% Store map in a struct
geo_map.X_grid = X;
geo_map.Y_grid = Y;
geo_map.Mag_map = Mag;

fprintf('Map generation complete.\n');


%% 3. INITIALIZATION
%==========================================================================
fprintf('Initializing simulation...\n');

% --- Initial True State ---
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; % [x, y, theta]

% --- Initial Particle Set ---
INIT_POS_STD = 5.0; % Std dev for initial particle positions
INIT_ANG_STD = 0.5; % Std dev for initial particle angles

particles = zeros(M_init, 3);
particles(:, 1) = true_state(1) + randn(M_init, 1) * INIT_POS_STD;
particles(:, 2) = true_state(2) + randn(M_init, 1) * INIT_POS_STD;
particles(:, 3) = true_state(3) + randn(M_init, 1) * INIT_ANG_STD;

% --- History Logs (for plotting) ---
full_true_path_history = zeros(NUM_STEPS, 3); % True [x, y, theta]
full_pdr_step_history = zeros(NUM_STEPS, 2);  % True [step_len, d_theta]
true_path_history = zeros(NUM_STEPS, 2);      % True [x, y]
estimated_path_history = zeros(NUM_STEPS, 2); % Estimated [x, y]

% --- Store Initial State (t=1) ---
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);
estimated_path_history(1, :) = true_path_history(1, :);

% --- Initialize Adaptive M Variables ---
current_M = M_init;                     % M_t (current particle count)
M_history = zeros(NUM_STEPS, 1);        % Log of M_t over time
M_history(1) = current_M;


%% 4. RUN SIMULATION
%==========================================================================
fprintf('Running %d simulation steps (DTW-APF)...\n', NUM_STEPS);
h_waitbar = waitbar(0, 'Running DTW-APF 2D Particle Filter...');

for t = 2:NUM_STEPS
   
    % --- 4a. Simulate True Motion (PDR) ---
    % Get the next true state and the PDR step that led to it
    [true_state, pdr_step] = Get_Next_Step_2D(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    
    % Log the true motion
    full_pdr_step_history(t, :) = pdr_step;
    full_true_path_history(t, :) = true_state;
    
    % --- 4b. Prepare Function Inputs ---
    % Get the history window for this step
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    
    % Input 1: PDR History (U_t)
    pdr_history_for_function = zeros(SEQUENCE_LEN, 2);
    actual_len = end_idx - start_idx + 1;
    pdr_history_for_function(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
    
    % Input 2: "Live" Sensor Sequence (Y_t)
    live_sequence = zeros(1, SEQUENCE_LEN);
    path_segment = zeros(SEQUENCE_LEN, 3);
    path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
    
    for k = 1:SEQUENCE_LEN
        pos_x = path_segment(k, 1);
        pos_y = path_segment(k, 2);
        
        if pos_x ~= 0 || pos_y ~= 0
            % Sample magnetic field from map: h(X_t)
            live_sequence(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                                       pos_x, pos_y, 'linear', 0);
        end
    end
    % Add sensor noise: + e_t
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD;
    
    % --- 4c. Call the 2D Particle Filter Step ---
    if t == 2 && ~exist('dtw', 'file')
        error('Function "dtw" not found. This script requires the Signal Processing Toolbox.');
    end
    
    % Get new particles, best guess, and the minimum DTW distance
    [particles_out, best_guess_state, dist] = Particle_Filter_DTW_Step_2D(particles, live_sequence, ...
                                    pdr_history_for_function, geo_map, process_noise, DTW_NOISE_STD);
    
    
    % --- 4d. APF Control: Re-injection and Adaptation ---
    
    % Particle Re-injection (prevents particle depletion)
    M_current_step = size(particles_out, 1);
    N_reset = round(M_current_step * 0.05); % 5% of particles
    indices_to_reset = randperm(M_current_step, N_reset);
    
    % Reset selected particles to be near the best guess
    particles_out(indices_to_reset, 1) = best_guess_state(1) + randn(N_reset, 1) * pos_std;
    particles_out(indices_to_reset, 2) = best_guess_state(2) + randn(N_reset, 1) * pos_std;
    particles_out(indices_to_reset, 3) = best_guess_state(3) + randn(N_reset, 1) * ang_std;

    % M-Adaptation (adjust particle count)
    % M_new = Adapt_Particle_Count_DTW(current_M, dist, ...
    %             APF.DTW_THRESH_LOW, APF.DTW_THRESH_HIGH, ...
    %             APF.M_min, APF.M_max);
    % M_new = Adapt_Particle_Count_DTW_KLD(particles_out, KLD.bin_size_xy, ...
    %                             MAP_X_LEN, MAP_Y_LEN, ...
    %                             KLD.epsilon, KLD.delta, ...
    %                             KLD.n_min, KLD.n_max)    % Apply the new particle count
    % 传入 d_min 和 DTW_THRESH_HIGH
M_new = Adapt_Particle_Count_DTW_KLD(particles_out, dist, KLD.bin_size_xy, ...
                                MAP_X_LEN, MAP_Y_LEN, ...
                                KLD.epsilon, KLD.delta, ...
                                KLD.n_min, KLD.n_max, ...
                                KLD.DTW_THRESH_HIGH);
    if M_new ~= current_M
        particles_next = Adjust_Particle_Set(particles_out, M_new);
    else
        particles_next = particles_out; % No change
    end
    
    % Update state for the next iteration
    particles = particles_next;
    current_M = size(particles, 1);
                                    
    % --- 4e. Store Results ---
    true_path_history(t, :) = true_state(1:2);
    estimated_path_history(t, :) = best_guess_state(1:2);
    M_history(t) = current_M;
    
    % --- 4f. Update Waitbar ---
    waitbar(t/NUM_STEPS, h_waitbar, sprintf('DTW-APF (M=%d, dist=%.1f)', current_M, dist));
end

close(h_waitbar);
fprintf('Simulation complete.\n');


%% 5. PLOTTING
%==========================================================================
figure('Position', [100, 100, 1600, 500]);

% --- Plot 1: 2D Path ---
subplot(1, 3, 1); 
imagesc(geo_map.X_grid(1,:), geo_map.Y_grid(:,1), geo_map.Mag_map);
hold on;
axis xy; 
colormap('jet');
colorbar;
title('2D Path Tracking (DTW-APF) - Random Path'); 
xlabel('X Position');
ylabel('Y Position');
plot(particles(:, 1), particles(:, 2), 'k.', 'MarkerSize', 2, 'DisplayName', 'Final Particle Cloud');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', 'True Path');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'DTW-APF Estimate'); 
plot(true_path_history(1, 1), true_path_history(1, 2), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'True Start');
plot(estimated_path_history(1, 1), estimated_path_history(1, 2), 'gs', 'MarkerSize', 12, 'DisplayName', 'Est. Start');
plot(true_path_history(end, 1), true_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'True End');
plot(estimated_path_history(end, 1), estimated_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'Est. End');
legend('show', 'Location', 'best');
axis equal; 
axis([0 MAP_X_LEN 0 MAP_Y_LEN]);

% --- Plot 2: Error ---
subplot(1, 3, 2); 
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('Localization Error (Euclidean Distance)');
xlabel('Time Step');
ylabel('Error (in meters/units)');
grid on;
ylim([0, Inf]);

% --- Plot 3: Adaptive Particle Count (M_t) ---
subplot(1, 3, 3);
plot(1:NUM_STEPS, M_history, 'b-', 'LineWidth', 2);
title('Adaptive Particle Count (M_t)');
xlabel('Time Step');
ylabel('Particle Count (M)');
grid on;
line([1, NUM_STEPS], [KLD.M_min, KLD.M_min], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'M_min');
line([1, NUM_STEPS], [KLD.M_max, KLD.M_max], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'M_max');
legend('M_t', 'Bounds', 'Location', 'best');
ylim([APF.M_min*0.8, APF.M_max*1.2]);




