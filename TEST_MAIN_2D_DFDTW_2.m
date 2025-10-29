% =========================================================================
% Main Script: 2D Particle Filter (PF) + DTW Simulation
% (Adaptive M, Forced Diversification + Robustness Analysis)
% =========================================================================
clc;
clear;
close all;

%% 1. Simulation Parameters
% --- Particle Count Adaptation (Elvira et al. approach) ---
INITIAL_N_PARTICLES = 1000;      % Increased initial particles
M_MIN = 100;                     % Minimum particles
M_MAX = 10000;                   % Maximum particles
K_SAMPLES = 15;                  % (K) Virtual observation samples
WINDOW_SIZE = 50;                % (W) Statistical window size
P_LOW = 0.9;                     % (p_l) Low p-value threshold (Increase M)
P_HIGH = 0.99;                   % (p_h) High p-value threshold (Decrease M) - Made less aggressive
ADAPT_M_FLAG = true;             % Flag to enable/disable M adaptation

% --- Forced Diversification (Jittering) ---
JITTER_POS_STD = 0.05;           % Position jitter std dev
JITTER_ANG_STD = deg2rad(0.5);   % Angle jitter std dev

% --- Simulation & Map ---
NUM_STEPS = 100;                 % Simulation steps
SEQUENCE_LEN = 50;               % Sequence length for DTW
MAP_X_LEN = 50;
MAP_Y_LEN = 50;

% --- Noise Parameters ---
SENSOR_NOISE_STD = 0.5;          % Sensor noise (for live_sequence & virtual obs)
DTW_NOISE_STD = 25.0;            % DTW weighting parameter (Initial value)
process_noise.step_std = 0.35;   % Particle process noise (step) > real (0.2)
process_noise.theta_std = deg2rad(12.5); % Particle process noise (angle) > real (10.0 deg)

%% 2. Generate 2D Geomagnetic Map
fprintf('Generating non-symmetric smooth geomagnetic map...\n');
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]);
Mag = imgaussfilt(Mag_raw, 3.0); % Gaussian smoothing
geo_map.X_grid = X;
geo_map.Y_grid = Y;
geo_map.Mag_map = Mag;
fprintf('Map generated. Size: %d x %d\n', MAP_X_LEN, MAP_Y_LEN);

%% 3. Initialization
fprintf('Initializing simulation...\n');
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; % Initial true state [x, y, theta]

% --- Adaptive M variables ---
M_t = INITIAL_N_PARTICLES;  % Current particle count
history_A_K = nan(WINDOW_SIZE, 1); % History window for A_K statistic

% --- Particles (Gaussian Ball Initialization) ---
INIT_POS_STD = 2.0;  % Initial position std dev
INIT_ANG_STD = 0.5;  % Initial angle std dev (rad)
particles = zeros(M_t, 3);
particles(:, 1) = true_state(1) + randn(M_t, 1) * INIT_POS_STD;
particles(:, 2) = true_state(2) + randn(M_t, 1) * INIT_POS_STD;
particles(:, 3) = true_state(3) + randn(M_t, 1) * INIT_ANG_STD;

% --- History Arrays ---
full_true_path_history = zeros(NUM_STEPS, 3);
full_pdr_step_history = zeros(NUM_STEPS, 2);
true_path_history = zeros(NUM_STEPS, 2);
estimated_path_history = zeros(NUM_STEPS, 2);
full_M_history = zeros(NUM_STEPS, 1);

% --- Robustness Analysis variables ---
t_robust = 50; % Time step for analysis
particles_state_t_minus_1 = []; % To store particle state before the robust step

% --- Store Initial State (t=1) ---
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);
estimated_path_history(1, :) = true_path_history(1, :); % Initialize estimate at true start
full_M_history(1) = M_t;

%% 4. Run Simulation
fprintf('Running %d simulation steps...\n', NUM_STEPS);
h_waitbar = waitbar(0, 'Running 2D Particle Filter (Adaptive M, Jittering)...');

for t = 2:NUM_STEPS

    % --- 4a. Simulate True Motion (PDR) ---
    [true_state, pdr_step] = get_next_step_random(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step;
    full_true_path_history(t, :) = true_state;

    % --- 4b. Prepare Filter Inputs (PDR History & Live Sequence) ---
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    pdr_history_for_function = zeros(SEQUENCE_LEN, 2);
    actual_len = end_idx - start_idx + 1;
    pdr_history_for_function(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);

    live_sequence = zeros(1, SEQUENCE_LEN);
    path_segment = zeros(SEQUENCE_LEN, 3);
    path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
    for k = 1:SEQUENCE_LEN
        pos_x = path_segment(k, 1); pos_y = path_segment(k, 2);
        if pos_x == 0 && pos_y == 0, live_sequence(k) = 0;
        else, live_sequence(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, pos_x, pos_y, 'linear', 0); end
    end
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD;

    % --- 4c. Adaptive M Evaluation (Predict Step for Statistic) ---
    if ADAPT_M_FLAG
        y_t_real = live_sequence(end); % Current real observation value
        pdr_step_t = pdr_history_for_function(end, :); % Current PDR step
        propagated_particles_for_eval = propagate_particles_manual(particles, pdr_step_t, process_noise, MAP_X_LEN, MAP_Y_LEN);
        
        % --- 4d. Generate Statistic A_K,M,t ---
        a_k_m_t = Generate_A_K_Statistic(propagated_particles_for_eval, y_t_real, geo_map, SENSOR_NOISE_STD, K_SAMPLES);
        
        % --- 4e. Update Statistic History Window ---
        history_A_K(mod(t-2, WINDOW_SIZE) + 1) = a_k_m_t;
    end

    % --- [Modified] Store State for Robustness Analysis ---
    if t == t_robust - 1
        particles_state_t_minus_1 = particles; % Save state *before* step t_robust is processed
    end
    % --- End Modified ---

    % --- 4f. Call Main Particle Filter Step ---
    [particles_out, best_guess, ~] = Particle_Filter_DTW_Step_2D(particles, live_sequence, ...
                                    pdr_history_for_function, geo_map, process_noise, DTW_NOISE_STD, ...
                                    JITTER_POS_STD, JITTER_ANG_STD);

    % --- 4g. Store Results ---
    true_path_history(t, :) = true_state(1:2);
    estimated_path_history(t, :) = best_guess(1:2);

    % --- 4h. Adapt Particle Count M_t ---
    if ADAPT_M_FLAG && (t > WINDOW_SIZE) && (mod(t, WINDOW_SIZE) == 0) % Evaluate every W steps
        M_t = Adapt_Particle_Count(size(particles,1), history_A_K, K_SAMPLES, WINDOW_SIZE, P_LOW, P_HIGH, M_MIN, M_MAX);
        particles = Adjust_Particle_Set(particles_out, M_t); % Adjust set size for next iteration
    else
        particles = particles_out; % Use output if not adapting
    end

    % --- 4i. Record M ---
    full_M_history(t) = size(particles, 1);

    waitbar(t/NUM_STEPS, h_waitbar);
end
close(h_waitbar);
fprintf('Simulation finished.\n');
fprintf('Average particle count (steps 2-%d): %.0f\n', NUM_STEPS, mean(full_M_history(2:end)));

%% 5. Robustness Analysis (Varying DTW_NOISE_STD at t_robust)
fprintf('Running robustness analysis at t = %d...\n', t_robust);
nominal_estimate_robust = estimated_path_history(t_robust, :); % Estimate from main run
param_name = 'DTW\_NOISE\_STD'; % Underscore escaped for title
nominal_dtw = DTW_NOISE_STD; % Get the initial value used
param_values = nominal_dtw * [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0];
robustness_errors = zeros(length(param_values), 1);

if isempty(particles_state_t_minus_1)
    warning('State before t_robust not captured. Skipping robustness analysis.');
else
    % Extract inputs for step t_robust
    start_idx_robust = max(1, t_robust - SEQUENCE_LEN + 1);
    end_idx_robust = t_robust;
    pdr_history_robust = zeros(SEQUENCE_LEN, 2);
    actual_len_robust = end_idx_robust - start_idx_robust + 1;
    pdr_history_robust(end-actual_len_robust+1:end, :) = full_pdr_step_history(start_idx_robust:end_idx_robust, :);

    live_sequence_robust = zeros(1, SEQUENCE_LEN);
    path_segment_robust = zeros(SEQUENCE_LEN, 3);
    path_segment_robust(end-actual_len_robust+1:end, :) = full_true_path_history(start_idx_robust:end_idx_robust, :);
    for k = 1:SEQUENCE_LEN
        pos_x_r = path_segment_robust(k, 1); pos_y_r = path_segment_robust(k, 2);
        if pos_x_r == 0 && pos_y_r == 0, live_sequence_robust(k) = 0;
        else, live_sequence_robust(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, pos_x_r, pos_y_r, 'linear', 0); end
    end
    live_sequence_robust = live_sequence_robust + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD; % Re-add noise for consistency (optional)

    for i = 1:length(param_values)
        current_dtw_noise = param_values(i);
        fprintf('  Testing %s = %.2f\n', param_name, current_dtw_noise);

        % Run ONE step of the filter with the modified parameter
        [~, best_guess_robust, ~] = Particle_Filter_DTW_Step_2D( ...
                                        particles_state_t_minus_1, ... % Use saved state
                                        live_sequence_robust, ...
                                        pdr_history_robust, ...
                                        geo_map, ...
                                        process_noise, ... % Keep process noise nominal
                                        current_dtw_noise, ... % Use varied DTW noise
                                        JITTER_POS_STD, ...
                                        JITTER_ANG_STD);

        % Calculate deviation from nominal estimate at t_robust
        robustness_errors(i) = norm(best_guess_robust(1:2) - nominal_estimate_robust);
    end
end
fprintf('Robustness analysis finished.\n');

%% 6. Plotting (2x2 Layout)
figure('Position', [100, 100, 1200, 1000]); % Adjusted figure size

% --- Plot 1: 2D Path Tracking ---
subplot(2, 2, 1);
imagesc(geo_map.X_grid(1,:), geo_map.Y_grid(:,1), geo_map.Mag_map);
hold on; axis xy; colormap('jet'); colorbar;
title('2D Path Tracking (Adaptive M + DTW + Jitter)');
xlabel('X Position'); ylabel('Y Position');
plot(particles(:, 1), particles(:, 2), 'k.', 'MarkerSize', 1, 'DisplayName', 'Final Particles'); % Reduced marker size
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2, 'DisplayName', 'True Path');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1, 'DisplayName', 'PF Estimate');
plot(true_path_history(1, 1), true_path_history(1, 2), 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start (True)');
plot(estimated_path_history(1, 1), estimated_path_history(1, 2), 'gs', 'MarkerSize', 10, 'DisplayName', 'Start (Est)');
plot(true_path_history(end, 1), true_path_history(end, 2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'End (True)');
plot(estimated_path_history(end, 1), estimated_path_history(end, 2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'End (Est)');
legend('Location', 'best'); axis equal; axis([0 MAP_X_LEN 0 MAP_Y_LEN]);

% --- Plot 2: Robustness Analysis ---
subplot(2, 2, 2);
if isempty(particles_state_t_minus_1)
    text(0.5, 0.5, 'Robustness analysis skipped.', 'HorizontalAlignment', 'center');
    title(sprintf('Robustness Analysis (Skipped)'));
else
    plot(param_values, robustness_errors, 'm-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
    hold on;
    plot(nominal_dtw, 0, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'k'); % Mark nominal value (error=0 by definition)
    title(sprintf('Robustness to %s at t=%d', param_name, t_robust));
    xlabel(sprintf('Parameter Value (%s)', param_name));
    ylabel(sprintf('Deviation from Nominal Estimate at t=%d', t_robust));
    legend('Deviation', 'Nominal Value', 'Location', 'best');
    grid on;
end

% --- Plot 3: Estimation Error ---
subplot(2, 2, 3);
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('Estimation Error (Euclidean Distance)'); xlabel('Time Step'); ylabel('Error');
grid on; ylim([0, Inf]);

% --- Plot 4: Particle Count M_t ---
subplot(2, 2, 4);
plot(1:NUM_STEPS, full_M_history, 'b-', 'LineWidth', 2);
hold on;
line([1, NUM_STEPS], [M_MIN, M_MIN], 'Color', 'r', 'LineStyle', '--');
line([1, NUM_STEPS], [M_MAX, M_MAX], 'Color', 'r', 'LineStyle', '--');
title('Adaptive Particle Count (M_t)'); xlabel('Time Step'); ylabel('Particle Count');
legend('M_t', 'Min/Max Bounds', 'Location', 'best'); grid on; ylim([0, M_MAX * 1.1]);

% =========================================================================
% Local Function Definitions (Copied from previous steps)
% =========================================================================
function geometric_map = Geometric_Map_Generator(D, map_size)
% [From user code]
switch D, case 1, if numel(map_size)<1, error('Size vector too short for D=1'); end, map_length=map_size(1); geometric_map=rand(1,map_length); case 2, if numel(map_size)<2, error('Size vector too short for D=2'); end, map_rows=map_size(1); map_cols=map_size(2); geometric_map=rand(map_rows,map_cols); otherwise, error('Invalid dimension D'); end
end
% =========================================================================
function [new_state, pdr_step] = get_next_step_random(old_state, MAP_X_LEN, MAP_Y_LEN)
% [From user code, with original noise values]
MEAN_STEP_LEN=1.0; STEP_LEN_STD=0.2; MEAN_TURN_RAD=0.0; TURN_STD_RAD=deg2rad(10.0); %<-- Kept original higher noise here
random_step_len=MEAN_STEP_LEN+randn()*STEP_LEN_STD; random_d_theta=MEAN_TURN_RAD+randn()*TURN_STD_RAD;
if random_step_len<0, random_step_len=0; end; pdr_step=[random_step_len, random_d_theta];
old_x=old_state(1); old_y=old_state(2); old_theta=old_state(3);
new_theta=old_theta+random_d_theta; new_x=old_x+sin(new_theta)*random_step_len; new_y=old_y+cos(new_theta)*random_step_len;
new_x=max(1,min(MAP_X_LEN,new_x)); new_y=max(1,min(MAP_Y_LEN,new_y)); new_state=[new_x,new_y,new_theta];
end
% =========================================================================
function [particles_out, best_guess, d_min] = Particle_Filter_DTW_Step_2D(particles_in, live_sequence, ...
                                            pdr_history, geo_map, process_noise, DTW_NOISE_STD, ...
                                            JITTER_POS_STD, JITTER_ANG_STD)
% [User code + Modifications + Jitter]
    M = size(particles_in, 1);
    L = size(pdr_history, 1);
    K = length(live_sequence); % Added K for sequence length

    % --- 1. Propagation (Prediction) ---
    last_pdr_step = pdr_history(end, :);
    step_noise = randn(M, 1) * process_noise.step_std;
    theta_noise = randn(M, 1) * process_noise.theta_std;
    steps = last_pdr_step(1) + step_noise;
    steps(steps<0) = 0; % Ensure steps are non-negative
    d_thetas = last_pdr_step(2) + theta_noise;
    old_x = particles_in(:, 1); old_y = particles_in(:, 2); old_theta = particles_in(:, 3);
    new_theta = old_theta + d_thetas;
    new_x = old_x + sin(new_theta) .* steps;
    new_y = old_y + cos(new_theta) .* steps;
    new_x = max(1, min(geo_map.X_grid(1,end), new_x));
    new_y = max(1, min(geo_map.Y_grid(end,1), new_y));
    propagated_particles = [new_x, new_y, new_theta];

    % --- 2. Weighting ---
    weights = zeros(M, 1);
    dtw_variance = DTW_NOISE_STD^2;
    all_distances = zeros(M, 1);

    for m = 1:M
        % --- 2a. Path Reconstruction ---
        particle_path = zeros(K, 2); % Use K (sequence length)
        current_pos = propagated_particles(m, 1:2); current_theta = propagated_particles(m, 3);
        
        % Reconstruct path segment matching live_sequence length
        temp_path_store = zeros(K, 3); % Store [x, y, theta] temporarily
        temp_path_store(end,:) = [current_pos, current_theta];

        for k_rev = K-1:-1:1
             pdr_idx = size(pdr_history,1) - (K-1-k_rev) -1; % Map k_rev to pdr_history index (adjusting for end offset)
             if pdr_idx < 1, continue; end % Handle start phase where history is shorter
             
             pdr_step_len = pdr_history(pdr_idx+1, 1); % Step that led *to* k_rev+1
             pdr_d_theta  = pdr_history(pdr_idx+1, 2); % Angle change that led *to* k_rev+1
             
             prev_theta = temp_path_store(k_rev+1, 3) - pdr_d_theta;
             curr_x = temp_path_store(k_rev+1, 1);
             curr_y = temp_path_store(k_rev+1, 2);
             curr_theta = temp_path_store(k_rev+1, 3); % Angle *at* k_rev+1
             
             prev_x = curr_x - sin(curr_theta) * pdr_step_len; % Use angle *at* end of step
             prev_y = curr_y - cos(curr_theta) * pdr_step_len;

             temp_path_store(k_rev, :) = [prev_x, prev_y, prev_theta];
        end
        particle_path = temp_path_store(:, 1:2); % Extract reconstructed x, y path

        % --- 2b. Generate Map Sequence ---
        map_sequence = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                              particle_path(:, 1), particle_path(:, 2), 'linear', 0);
        map_sequence(isnan(map_sequence)) = 0; % Handle NaNs from interp2 if path goes OOB during reconstruction

        % --- 2c. Calculate DTW Distance and Convert to Weight ---
         if size(live_sequence, 2) ~= K || size(map_sequence, 1) ~= K
             warning('Sequence length mismatch: live=%d, map=%d', size(live_sequence,2), size(map_sequence,1));
             distance = Inf; % Assign high distance if lengths mismatch
         else
             if exist('dtw', 'file') % Check if Signal Processing Toolbox is available
                 distance = dtw(live_sequence(:), map_sequence(:));
             else
                 distance = dtw_cost(live_sequence(:), map_sequence(:)); % Use fallback
             end
         end
        all_distances(m) = distance;
        weights(m) = exp(-distance^2 / (2 * dtw_variance));
    end
    % --- 2d. Normalize Weights ---
    sum_weights = sum(weights);
    if sum_weights > 1e-15, weights = weights / sum_weights; else weights = ones(M, 1) / M; end

    % --- 3. Resampling (Systematic) ---
    resampled_particles = zeros(M, 3);
    cdf = cumsum(weights);
    r_0 = rand() / M;
    idx_j = 1;
    for m = 1:M
        U = r_0 + (m-1)/M;
        while idx_j <= M && U > cdf(idx_j) % Added boundary check
             idx_j = idx_j + 1;
        end
        if idx_j > M, idx_j = M; end % Ensure index does not exceed M
        resampled_particles(m, :) = propagated_particles(idx_j, :);
    end

    % --- 3b. Forced Diversification (Jittering) ---
    jitter_pos = randn(M, 2) * JITTER_POS_STD;
    jitter_ang = randn(M, 1) * JITTER_ANG_STD;
    particles_out = resampled_particles;
    particles_out(:, 1:2) = particles_out(:, 1:2) + jitter_pos;
    particles_out(:, 3) = particles_out(:, 3) + jitter_ang;
    particles_out(:, 1) = max(1, min(geo_map.X_grid(1,end), particles_out(:, 1)));
    particles_out(:, 2) = max(1, min(geo_map.Y_grid(end,1), particles_out(:, 2)));

    % --- 4. Estimation ---
    best_guess(1:2) = mean(particles_out(:, 1:2), 1);
    best_guess(3) = atan2(mean(sin(particles_out(:, 3))), mean(cos(particles_out(:, 3)))); % Circular mean for angle
    d_min = min(all_distances); % Minimum DTW distance in this step
end
% =========================================================================
function [propagated_particles] = propagate_particles_manual(particles_in, pdr_step, process_noise, MAP_X_LEN, MAP_Y_LEN)
% [New function for Adaptive M statistic]
M=size(particles_in,1); step_noise=randn(M,1)*process_noise.step_std; theta_noise=randn(M,1)*process_noise.theta_std;
steps=pdr_step(1)+step_noise; steps(steps<0)=0; d_thetas=pdr_step(2)+theta_noise;
old_x=particles_in(:,1); old_y=particles_in(:,2); old_theta=particles_in(:,3);
new_theta=old_theta+d_thetas; new_x=old_x+sin(new_theta).*steps; new_y=old_y+cos(new_theta).*steps;
new_x=max(1,min(MAP_X_LEN,new_x)); new_y=max(1,min(MAP_Y_LEN,new_y)); propagated_particles=[new_x,new_y,new_theta];
end
% =========================================================================
function [a_k_m_t] = Generate_A_K_Statistic(propagated_particles, y_t_real, geo_map, SENSOR_NOISE_STD, K)
% [New function for Adaptive M statistic] (Algorithm 2)
M_t=size(propagated_particles,1); if M_t==0, a_k_m_t=0; return; end % Handle empty case
indices=randi(M_t,K,1); % Draw K samples with replacement
predicted_means=interp2(geo_map.X_grid,geo_map.Y_grid,geo_map.Mag_map,propagated_particles(indices,1),propagated_particles(indices,2),'linear',0);
predicted_means(isnan(predicted_means)) = 0; % Handle NaNs
fictitious_obs=predicted_means+randn(K,1)*SENSOR_NOISE_STD; % Generate K virtual observations
a_k_m_t=sum(fictitious_obs<y_t_real); % Count how many are less than the real observation
end
% =========================================================================
function [new_M] = Adapt_Particle_Count(current_M, history_A_K, K, W, p_low, p_high, M_min, M_max)
% [New function for Adaptive M] (Algorithm 3, modified adaptation logic)
valid_history = history_A_K(~isnan(history_A_K)); % Use only valid entries in window
if length(valid_history) < W / 2, new_M = current_M; return; end % Don't adapt if window not filled enough
W_actual = length(valid_history); % Actual number of samples in window
O_j=zeros(K+1,1); for j=0:K, O_j(j+1)=sum(valid_history==j); end;
E_j=W_actual/(K+1); % Expected count based on actual window size
chi2_stat=sum((O_j-E_j).^2/(E_j+1e-9)); % Chi-squared statistic
p_value=1-chi2cdf(chi2_stat,K); % P-value
if p_value<=p_low % If p-value is low (observed distribution deviates significantly from uniform) -> Increase M
    new_M=min(round(current_M*1.5), M_max); % Increase by 50%
    fprintf(' (t=%d, p=%.3f <= %.2f, W=%d) Increasing M -> %d\n', evalin('base','t'), p_value, p_low, W_actual, new_M);
elseif p_value>=p_high % If p-value is high (observed distribution is very close to uniform) -> Decrease M
    new_M=max(round(current_M/1.5), M_min); % Decrease by 33%
    fprintf(' (t=%d, p=%.3f >= %.2f, W=%d) Decreasing M -> %d\n', evalin('base','t'), p_value, p_high, W_actual, new_M);
else % Otherwise, keep M the same
    new_M=current_M;
    % fprintf(' (t=%d, p=%.3f, W=%d) Keeping M = %d\n', evalin('base','t'), p_value, W_actual, new_M);
end
end
% =========================================================================
function [new_particle_set] = Adjust_Particle_Set(old_particle_set, new_M)
% [New function for Adaptive M] (Algorithm 4)
current_M=size(old_particle_set,1); if new_M==current_M, new_particle_set=old_particle_set; return; end
if new_M>current_M % Increase particles
    num_to_add=new_M-current_M;
    indices_to_add=randi(current_M,num_to_add,1); % Sample with replacement from existing
    new_particles=old_particle_set(indices_to_add,:);
    % Add small jitter to new particles to avoid exact duplicates
    jitter_pos=randn(num_to_add,2)*0.01; % Very small jitter
    jitter_ang=randn(num_to_add,1)*0.001;
    new_particles(:,1:2)=new_particles(:,1:2)+jitter_pos;
    new_particles(:,3)=new_particles(:,3)+jitter_ang;
    new_particle_set=[old_particle_set; new_particles];
else % Decrease particles
    indices_to_keep=randperm(current_M,new_M); % Randomly select which to keep
    new_particle_set=old_particle_set(indices_to_keep,:);
end
end
% =========================================================================
function cost = dtw_cost(s1, s2)
% [New fallback DTW function]
n=length(s1); m=length(s2); if n==0 || m==0, cost=Inf; return; end % Handle empty input
D=zeros(n+1,m+1); D(:,1)=Inf; D(1,:)=Inf; D(1,1)=0; s1=s1(:); s2=s2(:);
for i=2:n+1, for j=2:m+1, cost_val=abs(s1(i-1)-s2(j-1)); D(i,j)=cost_val+min([D(i-1,j),D(i,j-1),D(i-1,j-1)]); end; end; cost=D(n+1,m+1);
end
% =========================================================================