%==========================================================================
% SCRIPT: MAIN_2D_PF_KLD_v2
%
% DESCRIPTION:
%   Main simulation script for a 2D Adaptive Particle Filter (APF)
%   aided by PDR and geomagnetic matching using DTW.
%   - Particle count (n) is adapted using KLD-SAMPLING.
%   - Particle re-injection is used to combat particle impoverishment.
%
% DEPENDENCIES:
%   - Signal Processing Toolbox (for 'dtw' function)
%   - calculate_n_kld.m
%   - propagate_one_particle.m
%   - weight_one_particle.m
%   - get_bin_for_particle.m
%   - get_next_step_random.m
%==========================================================================
clc;
clear;
close all;
%% 1. SIMULATION PARAMETERS
%==========================================================================
M_init = 500;            
NUM_STEPS = 100;          
SEQUENCE_LEN = 50;        
MAP_X_LEN = 50;           
MAP_Y_LEN = 50;           
SENSOR_NOISE_STD = 0.5;   
DTW_NOISE_STD = 10; % [修改] 使用一个更有区分度的值
process_noise.step_std = 0.5;          
process_noise.theta_std = deg2rad(10); 

% --- KLD-Sampling Parameters ---
KLD.epsilon = 0.05;       % 允许的最大近似误差
KLD.delta = 0.01;         % 置信度: (1-delta) = 99%
KLD.n_min = 500;          % 最小粒子数
KLD.n_max = 100000;       % 最大粒子数 (安全上限)

% --- KLD Binning Parameters ---
BIN.xy_size = 1.0;        % [修改] 增大箱子尺寸 (米)
BIN.theta_size = deg2rad(15); % [修改] 增大箱子尺寸 (弧度)

% --- Re-injection Parameters ---
RESET.percentage = 0.05;  % 重注入 5% 的粒子
RESET.pos_std = 3.0;      % 位置抖动标准差
RESET.ang_std = deg2rad(5); % 角度抖动标准差

%% 2. MAP GENERATION
%==========================================================================
fprintf('Generating geomagnetic map...\n');
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
Mag = imgaussfilt(Mag_raw, 3.0); 
geo_map.X_grid = X;
geo_map.Y_grid = Y;
geo_map.Mag_map = Mag;
fprintf('Map generation complete.\n');
%% 3. INITIALIZATION
%==========================================================================
fprintf('Initializing simulation...\n');
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; 
INIT_POS_STD = 5.0; 
INIT_ANG_STD = 0.5; 

% [修改] 粒子现在是结构体数组 (Struct Array)
particles(1:M_init) = struct('x', 0, 'y', 0, 'theta', 0, 'weight', 1/M_init);
for i = 1:M_init
    particles(i).x = true_state(1) + randn() * INIT_POS_STD;
    particles(i).y = true_state(2) + randn() * INIT_POS_STD;
    particles(i).theta = true_state(3) + randn() * INIT_ANG_STD;
end

full_true_path_history = zeros(NUM_STEPS, 3); 
full_pdr_step_history = zeros(NUM_STEPS, 2);  
true_path_history = zeros(NUM_STEPS, 2);      
estimated_path_history = zeros(NUM_STEPS, 2); 
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);
estimated_path_history(1, :) = true_path_history(1, :);
M_history = zeros(NUM_STEPS, 1);        
M_history(1) = M_init;
dist_history = zeros(NUM_STEPS, 1);
%% 4. RUN SIMULATION
%==========================================================================
fprintf('Running %d simulation steps (KLD-PF-DTW + Re-injection)...\n', NUM_STEPS);
h_waitbar = waitbar(0, 'Running KLD-PF 2D Particle Filter...');

for t = 2:NUM_STEPS
   
    % --- 4a. Simulate True Motion (PDR) ---
    [true_state, pdr_step] = Get_Next_Step_2D(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step;
    full_true_path_history(t, :) = true_state;
    
    % --- 4b. Prepare Function Inputs ---
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    pdr_history_for_function = zeros(SEQUENCE_LEN, 2);
    actual_len = end_idx - start_idx + 1;
    pdr_history_for_function(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
    
    live_sequence = zeros(1, SEQUENCE_LEN);
    path_segment = zeros(SEQUENCE_LEN, 3);
    path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
    
    for k_seq = 1:SEQUENCE_LEN
        pos_x = path_segment(k_seq, 1);
        pos_y = path_segment(k_seq, 2);
        if pos_x ~= 0 || pos_y ~= 0
            live_sequence(k_seq) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                                       pos_x, pos_y, 'linear', 0);
        end
    end
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD;
    
    % --- 4c. KLD-Sampling 循环 ---
    
    n = 0;
    k = 0;
    n_chi = Inf;
    occupied_bins = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    
    S_t_minus_1 = particles; 
    weights_prev = [S_t_minus_1.weight];
    
    sum_weights_prev = sum(weights_prev);
    if sum_weights_prev < 1e-15
        weights_prev = ones(1, length(S_t_minus_1)) / length(S_t_minus_1);
    else
        weights_prev = weights_prev / sum_weights_prev;
    end
    cdf_prev = cumsum(weights_prev);
    
    S_t = struct('x', {}, 'y', {}, 'theta', {}, 'weight', {}); 
    total_weight_alpha = 0; 
    d_min = Inf;            
    
    last_pdr_step_for_loop = full_pdr_step_history(t, :);

    while (n < n_chi || n < KLD.n_min) && n < KLD.n_max
        
        n = n + 1; 
        
        idx_j = 1 + sum(cdf_prev < rand()); 
        p_old = S_t_minus_1(idx_j);
        
        p_new = propagate_one_particle(p_old, last_pdr_step_for_loop, process_noise, geo_map);
        
        [weight, dtw_dist] = weight_one_particle(p_new, live_sequence, ...
            pdr_history_for_function, geo_map, DTW_NOISE_STD);
        
        if dtw_dist < d_min
            d_min = dtw_dist; 
        end
        
        S_t(n).x = p_new.x;
        S_t(n).y = p_new.y;
        S_t(n).theta = p_new.theta;
        S_t(n).weight = weight;
        
        total_weight_alpha = total_weight_alpha + weight; 

        bin_id = get_bin_for_particle(p_new, BIN.xy_size, BIN.theta_size);
        if ~isKey(occupied_bins, bin_id)
            occupied_bins(bin_id) = true;
            k = occupied_bins.Count;
        end
        
        if k > 1 && n >= KLD.n_min
            n_chi = calculate_n_kld(k, KLD.epsilon, KLD.delta);
        end
    end 
    
    % --- 4d. 归一化, 估计, 重注入 ---
    
    % --- 4d. 归一化, 估计, 重注入 ---
    
    n_out = n;
    
    if total_weight_alpha > 1e-15
        normalized_weights = [S_t.weight] / total_weight_alpha;
    else
        normalized_weights = ones(1, n_out) / n_out;
    end
    for i = 1:n_out
        S_t(i).weight = normalized_weights(i);
    end
    
    % [修改] 估计最佳猜测 (现在使用加权平均)
    mean_x = sum([S_t.x] .* [S_t.weight]);
    mean_y = sum([S_t.y] .* [S_t.weight]);
    mean_sin_th = sum(sin([S_t.theta]) .* [S_t.weight]);
    mean_cos_th = sum(cos([S_t.theta]) .* [S_t.weight]);
    best_guess_state = [mean_x, mean_y, atan2(mean_sin_th, mean_cos_th)];
    
    
    % *** [关键修改：全局重注入] ***
    N_reset = round(n_out * RESET.percentage);
    indices_to_reset = randperm(n_out, N_reset);
    indices_to_keep = setdiff(1:n_out, indices_to_reset);
    
    % 1. 重置选中的粒子 (随机撒满全图)
    for idx = indices_to_reset
        S_t(idx).x = rand() * MAP_X_LEN; % 随机 X
        S_t(idx).y = rand() * MAP_Y_LEN; % 随机 Y
        S_t(idx).theta = rand() * 2 * pi;  % 随机角度
        S_t(idx).weight = 1.0 / n_out; % 赋予平均权重
    end
    
    % 2. 重新缩放未重置的粒子权重，使总和为 1
    weight_sum_kept = sum([S_t(indices_to_keep).weight]);
    weight_sum_reset = N_reset / n_out;
    
    if weight_sum_kept > 1e-15
        scale_factor = (1.0 - weight_sum_reset) / weight_sum_kept;
        for idx = indices_to_keep
            S_t(idx).weight = S_t(idx).weight * scale_factor;
        end
    else
        % 如果所有保留的粒子权重都为0 (不太可能但需处理)
        for i = 1:n_out
            S_t(i).weight = 1.0 / n_out;
        end
    end
    % *** [修改结束] ***
    
    particles_out = S_t;
                                    
    % --- 4e. Store Results ---
    % ... (后续代码不变) ...
                                    
    % --- 4e. Store Results ---
    particles = particles_out; 
    current_M = n_out;
    
    true_path_history(t, :) = true_state(1:2);
    estimated_path_history(t, :) = best_guess_state(1:2);
    M_history(t) = current_M;
    dist_history(t) = d_min;
    
    % --- 4f. Update Waitbar ---
    waitbar(t/NUM_STEPS, h_waitbar, sprintf('KLD-PF (M=%d, dist=%.1f)', current_M, d_min));
end
close(h_waitbar);
fprintf('Simulation complete.\n');
%% 5. PLOTTING
%==========================================================================
figure('Position', [100, 100, 1800, 500]);
% --- Plot 1: 2D Path ---
subplot(1, 4, 1); 
imagesc(geo_map.X_grid(1,:), geo_map.Y_grid(:,1), geo_map.Mag_map);
hold on;
axis xy; 
colormap('jet');
colorbar;
title('2D Path Tracking (KLD-PF + Re-injection)'); 
xlabel('X Position');
ylabel('Y Position');
plot([particles.x], [particles.y], 'k.', 'MarkerSize', 2, 'DisplayName', 'Final Particle Cloud');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', 'True Path');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'KLD-PF Estimate'); 
plot(true_path_history(1, 1), true_path_history(1, 2), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'True Start');
plot(estimated_path_history(1, 1), estimated_path_history(1, 2), 'gs', 'MarkerSize', 12, 'DisplayName', 'Est. Start');
plot(true_path_history(end, 1), true_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'True End');
plot(estimated_path_history(end, 1), estimated_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'Est. End');
legend('show', 'Location', 'best');
axis equal; 
axis([0 MAP_X_LEN 0 MAP_Y_LEN]);
% --- Plot 2: Error ---
subplot(1, 4, 2); 
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('Localization Error (Euclidean Distance)');
xlabel('Time Step');
ylabel('Error (in meters/units)');
grid on;
ylim([0, Inf]);
% --- Plot 3: Adaptive Particle Count (M_t) ---
subplot(1, 4, 3);
plot(1:NUM_STEPS, M_history, 'b-', 'LineWidth', 2);
title('KLD-Sampling Particle Count (M_t)');
xlabel('Time Step');
ylabel('Particle Count (M)');
grid on;
line([1, NUM_STEPS], [KLD.n_min, KLD.n_min], 'Color', 'red', 'LineStyle', '--');
line([1, NUM_STEPS], [KLD.n_max, KLD.n_max], 'Color', 'red', 'LineStyle', '--');
legend('M_t (KLD)', 'Bounds', 'Location', 'best');
ylim([KLD.n_min*0.8, KLD.n_max*1.2]);
% --- Plot 4: Min DTW Distance ---
subplot(1, 4, 4);
plot(1:NUM_STEPS, dist_history, 'm-', 'LineWidth', 1.5);
title('Minimum DTW Distance');
xlabel('Time Step');
ylabel('DTW Distance');
grid on;
ylim([0, Inf]);
%% ========================================================================
%  FILE: calculate_n_kld.m
%  ========================================================================
function n_chi = calculate_n_kld(k, epsilon, delta)
    if k <= 1
        n_chi = Inf;
        return;
    end
    z_1_minus_delta = norminv(1 - delta, 0, 1);
    k_minus_1 = double(k - 1);
    term_1 = 2.0 / (9.0 * k_minus_1);
    term_2 = sqrt(term_1) * z_1_minus_delta;
    inner_term_cubed = (1.0 - term_1 + term_2)^3;
    chi_squared = k_minus_1 * inner_term_cubed;
    n_chi = chi_squared / (2.0 * epsilon);
    n_chi = ceil(n_chi);
end

%% ========================================================================
%  FILE: propagate_one_particle.m
%  ========================================================================
function p_out = propagate_one_particle(p_in, last_pdr_step, process_noise, geo_map)
    step_noise = randn() * process_noise.step_std;
    theta_noise = randn() * process_noise.theta_std;
    
    steps = last_pdr_step(1) + step_noise;
    d_theta = last_pdr_step(2) + theta_noise;
    
    old_x = p_in.x;
    old_y = p_in.y;
    old_theta = p_in.theta;
    
    new_theta = old_theta + d_theta;
    new_x = old_x + sin(new_theta) * steps;
    new_y = old_y + cos(new_theta) * steps;
    
    new_x = max(1, min(geo_map.X_grid(1,end), new_x));
    new_y = max(1, min(geo_map.Y_grid(end,1), new_y));
    
    p_out.x = new_x;
    p_out.y = new_y;
    p_out.theta = new_theta;
end

%% ========================================================================
%  FILE: weight_one_particle.m
%  ========================================================================
function [weight, dtw_dist] = weight_one_particle(p_propagated, live_sequence, pdr_history, geo_map, DTW_NOISE_STD)
    L = size(pdr_history, 1);
    
    particle_path = zeros(L, 2); 
    current_pos = [p_propagated.x, p_propagated.y];
    current_theta = p_propagated.theta;
    particle_path(end, :) = current_pos;
    
    for k = L-1:-1:1
        pdr_step_len = pdr_history(k+1, 1);
        pdr_d_theta = pdr_history(k+1, 2);
        prev_theta = current_theta - pdr_d_theta;
        prev_x = current_pos(1) - sin(current_theta) * pdr_step_len;
        prev_y = current_pos(2) - cos(current_theta) * pdr_step_len;
        particle_path(k, :) = [prev_x, prev_y];
        current_pos = [prev_x, prev_y];
        current_theta = prev_theta;
    end
    
    map_sequence = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                          particle_path(:, 1), particle_path(:, 2), 'linear', 0);
                      
    dtw_dist = dtw(live_sequence(:), map_sequence(:));
    
    dtw_variance = DTW_NOISE_STD^2;
    weight = exp(-dtw_dist^2 / (2 * dtw_variance));
end

%% ========================================================================
%  FILE: get_bin_for_particle.m
%  ========================================================================
function bin_id = get_bin_for_particle(particle, BIN_SIZE_XY, BIN_SIZE_THETA)
    x_bin = floor(particle.x / BIN_SIZE_XY);
    y_bin = floor(particle.y / BIN_SIZE_XY);
    
    theta_norm = mod(particle.theta, 2*pi); 
    theta_bin = floor(theta_norm / BIN_SIZE_THETA);
    
    bin_id = sprintf('%d_%d_%d', x_bin, y_bin, theta_bin);
end