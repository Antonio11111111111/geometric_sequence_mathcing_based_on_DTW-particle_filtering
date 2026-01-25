function [particles_out, best_guess, d_min] = Particle_Filter_DTW_Step_2D_CPFTest(particles_in_gpu, live_sequence_gpu, ...
                                            pdr_history_cpu, geo_map_gpu, process_noise, DTW_NOISE_STD)
%
% Function: Particle_Filter_DTW_Step_2D_CPFTest
%
% Summary: [GPU/DDTW + Standard Resampling VERSION]
%          执行一步2D粒子滤波。
%          此版本使用 DDTW (导数DTW) 进行权重计算 [cite: 280]，并使用
%          标准系统重采样 (SIR) 以确保稳定性。
%
% 流程:
%          1. Propagate (传播)
%          2. Weight (使用 DDTW)
%          3. Standard Systematic Resampling (标准重采样)
%
% -------------------------------------------------------------------------
    
    M = size(particles_in_gpu, 1);
    L = size(pdr_history_cpu, 1);
    
    % --- 1. Propagation (Prediction) [ON GPU] ---
    last_pdr_step_cpu = single(pdr_history_cpu(end, :)); 
    step_noise_gpu = randn(M, 1, 'single', 'gpuArray') * single(process_noise.step_std);
    theta_noise_gpu = randn(M, 1, 'single', 'gpuArray') * single(process_noise.theta_std);
    steps_gpu = last_pdr_step_cpu(1) + step_noise_gpu;
    d_thetas_gpu = last_pdr_step_cpu(2) + theta_noise_gpu;
    old_theta_gpu = particles_in_gpu(:, 3);
    new_theta_gpu = old_theta_gpu + d_thetas_gpu;
    new_x_gpu = particles_in_gpu(:, 1) + sin(new_theta_gpu) .* steps_gpu;
    new_y_gpu = particles_in_gpu(:, 2) + cos(new_theta_gpu) .* steps_gpu;
    
    % Boundary check on GPU
    x_min = geo_map_gpu.X_grid(1, 1);
    x_max = geo_map_gpu.X_grid(1, end);
    y_min = geo_map_gpu.Y_grid(1, 1);
    y_max = geo_map_gpu.Y_grid(end, 1);
    new_x_gpu = max(x_min, min(x_max, new_x_gpu));
    new_y_gpu = max(y_min, min(y_max, new_y_gpu));
    
    propagated_particles_gpu = [new_x_gpu, new_y_gpu, new_theta_gpu];
    
    % --- 2. Weighting ---
    % --- 2a. Path Reconstruction (Vectorized on GPU) ---
    all_particle_paths_x_gpu = zeros(M, L, 'single', 'gpuArray');
    all_particle_paths_y_gpu = zeros(M, L, 'single', 'gpuArray');
    all_particle_thetas_gpu = zeros(M, L, 'single', 'gpuArray');
    all_particle_paths_x_gpu(:, L) = propagated_particles_gpu(:, 1);
    all_particle_paths_y_gpu(:, L) = propagated_particles_gpu(:, 2);
    all_particle_thetas_gpu(:, L) = propagated_particles_gpu(:, 3);
    
    for k = L-1:-1:1
        pdr_step_len = single(pdr_history_cpu(k+1, 1));
        pdr_d_theta = single(pdr_history_cpu(k+1, 2));
        current_thetas_gpu = all_particle_thetas_gpu(:, k+1);
        current_X_gpu = all_particle_paths_x_gpu(:, k+1);
        current_Y_gpu = all_particle_paths_y_gpu(:, k+1);
        prev_thetas_gpu = current_thetas_gpu - pdr_d_theta;
        prev_X_gpu = current_X_gpu - sin(current_thetas_gpu) * pdr_step_len;
        prev_Y_gpu = current_Y_gpu - cos(current_thetas_gpu) * pdr_step_len;
        all_particle_thetas_gpu(:, k) = prev_thetas_gpu;
        all_particle_paths_x_gpu(:, k) = prev_X_gpu;
        all_particle_paths_y_gpu(:, k) = prev_Y_gpu;
    end
    
    % --- 2b. Generate Map Sequences (Vectorized on GPU) ---
    all_map_sequences_gpu = interp2(geo_map_gpu.X_grid, geo_map_gpu.Y_grid, geo_map_gpu.Mag_map, ...
                                  all_particle_paths_x_gpu, all_particle_paths_y_gpu, 'linear', 0);
                              
    % --- 2c. Calculate DTW & Weights (Parallel on CPU) ---
    all_map_sequences_cpu = gather(all_map_sequences_gpu);
    live_sequence_cpu = gather(live_sequence_gpu);
    live_sequence_cpu = live_sequence_cpu(:)'; 
    
    % --- [DDTW MODIFICATION] ---
    % 计算 "live" 序列的导数 (一阶差分)，用于DDTW [cite: 280]
    % 序列长度将变为 L-1
    live_sequence_deriv = diff(live_sequence_cpu);
    % --- [END MODIFICATION] ---
    
    weights_cpu = zeros(M, 1, 'single');
    all_distances_cpu = zeros(M, 1, 'single');
    dtw_variance = single(DTW_NOISE_STD^2); 
    
    parfor m = 1:M
        map_seq = all_map_sequences_cpu(m, :);

        % --- [DDTW MODIFICATION] ---
        % 计算 "map" 序列的导数
        map_seq_deriv = diff(map_seq);

        % 使用导数序列 (形状) 计算DTW [cite: 280]
        distance = dtw(live_sequence_deriv, map_seq_deriv);
        % --- [END MODIFICATION] ---

        all_distances_cpu(m) = distance;
        weights_cpu(m) = 1./(1+(distance./dtw_variance).^2);
    end
    % --- [NEW: Multi-Window Ensemble Settings] ---
    % 定义您的三个尺度 (L-1 是因为 diff() 减少了1)
    % L_full = length(live_sequence_deriv); % e.g., 99
    % L_med  = 50;  % 中尺度 (示例)
    % L_short = 20; % 短尺度 (示例)
    % 
    % % (确保 L_short 和 L_med 不大于 L_full)
    % L_med = min(L_med, L_full);
    % L_short = min(L_short, L_full);
    % 
    % % 提取 "live" 序列的末端
    % live_deriv_full = live_sequence_deriv;
    % live_deriv_med  = live_sequence_deriv(end-L_med+1:end);
    % live_deriv_short = live_sequence_deriv(end-L_short+1:end);
    % % --- [END NEW] ---
    % 
    % parfor m = 1:M
    %     map_seq_full = all_map_sequences_cpu(m, :);
    % 
    %     % --- [DDTW MODIFICATION] ---
    %     map_seq_deriv_full = diff(map_seq_full);
    % 
    %     % 1. 提取地图序列的末端
    %     map_deriv_med  = map_seq_deriv_full(end-L_med+1:end);
    %     map_deriv_short = map_seq_deriv_full(end-L_short+1:end);
    % 
    %     % 2. 在所有三个尺度上计算DTW
    %     dist_long  = dtw(live_deriv_full, map_seq_deriv_full);
    %     dist_med   = dtw(live_deriv_med, map_deriv_med);
    %     dist_short = dtw(live_deriv_short, map_deriv_short);
    % 
    %     % 3. 计算三个独立的权重 (使用您现有的Cauchy函数)
    %     w_long  = 1./(1+(dist_long ./dtw_variance).^2);
    %     w_med   = 1./(1+(dist_med  ./dtw_variance).^2);
    %     w_short = 1./(1+(dist_short./dtw_variance).^2);
    % 
    %     % 4. [关键] 集成权重：必须在所有尺度上都匹配良好
    %     weights_cpu(m) = w_long * w_med * w_short; 
    % 
    %     % 仍然只记录最长距离用于调试
    %     all_distances_cpu(m) = dist_long;
    %     % --- [END MODIFICATION] ---
    % end
    
    % --- 2d. Normalize Weights (ON GPU) ---
    weights_gpu = gpuArray(weights_cpu);
    all_distances_gpu = gpuArray(all_distances_cpu);
    
    sum_weights = sum(weights_gpu);
    if sum_weights > 1e-15
        weights_gpu = weights_gpu / sum_weights;
    else
        weights_gpu = ones(M, 1, 'single', 'gpuArray') / M;
    end
    
    %% --- 3. Standard Systematic Resampling [ON GPU] ---
    % (已恢复为标准SIR滤波器以确保稳定性)
    
    cdf_gpu = cumsum(weights_gpu);
    r_0_gpu = rand(1, 'single', 'gpuArray') / M;
    U_gpu = r_0_gpu + (gpuArray(single(transpose(1:M))) - 1) / M;
    
    [~, indices_gpu] = histc(U_gpu, [0; cdf_gpu]);
    indices_gpu = indices_gpu + 1; 
    indices_gpu = min(indices_gpu, M); 
    
    particles_out = propagated_particles_gpu(indices_gpu, :);
    
    % --- 3d. Final Boundary Check (on GPU) ---
    particles_out(:, 1) = max(x_min, min(x_max, particles_out(:, 1)));
    particles_out(:, 2) = max(y_min, min(y_max, particles_out(:, 2)));
    particles_out(:, 3) = mod(particles_out(:, 3) + pi, 2*pi) - pi;

    % --- 4. Estimation (ON GPU) ---
    mean_pos_gpu = mean(particles_out(:, 1:2), 1);
    mean_sin_gpu = mean(sin(particles_out(:, 3)));
    mean_cos_gpu = mean(cos(particles_out(:, 3)));
    mean_theta_gpu = atan2(mean_sin_gpu, mean_cos_gpu);
    
    best_guess = [mean_pos_gpu, mean_theta_gpu];
    d_min = min(all_distances_gpu);
end