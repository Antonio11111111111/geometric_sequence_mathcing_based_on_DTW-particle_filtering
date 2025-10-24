%% [修改] 函数签名现在返回第三个参数: d_min
function [particles_out, best_guess, d_min] = Particle_Filter_DTW_Step_2D(particles_in, live_sequence, ...
                                            pdr_history, geo_map, process_noise, DTW_NOISE_STD)
% 摘要: 执行一步粒子滤波 (传播、加权、重采样)
% 权重: 基于DTW（动态时间规整）
% 注意: 此函数需要 Signal Processing Toolbox (用于 'dtw' 函数)
    
    M = size(particles_in, 1);
    L = size(pdr_history, 1);
    
    % --- 1. 传播 (Prediction/Propagation) ---
    % ( ... 您的传播代码 ... )
    last_pdr_step = pdr_history(end, :);
    step_noise = randn(M, 1) * process_noise.step_std;
    theta_noise = randn(M, 1) * process_noise.theta_std;
    steps = last_pdr_step(1) + step_noise;
    d_thetas = last_pdr_step(2) + theta_noise;
    old_x = particles_in(:, 1);
    old_y = particles_in(:, 2);
    old_theta = particles_in(:, 3);
    new_theta = old_theta + d_thetas;
    new_x = old_x + sin(new_theta) .* steps;
    new_y = old_y + cos(new_theta) .* steps;
    new_x = max(1, min(geo_map.X_grid(1,end), new_x));
    new_y = max(1, min(geo_map.Y_grid(end,1), new_y));
    propagated_particles = [new_x, new_y, new_theta];

    % --- 2. 加权 (Weighting) ---
    weights = zeros(M, 1);
    dtw_variance = DTW_NOISE_STD^2;
    
    %% [修改] 新增一个数组来存储所有粒子的DTW距离
    all_distances = zeros(M, 1); 
    
    for m = 1:M
        % --- 2a. 重建粒子路径 (Path Reconstruction) ---
        % ( ... 您的路径重建代码 ... )
        particle_path = zeros(L, 2); 
        current_pos = propagated_particles(m, 1:2);
        current_theta = propagated_particles(m, 3);
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
        
        % --- 2b. 生成地图序列 ---
        map_sequence = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                              particle_path(:, 1), particle_path(:, 2), 'linear', 0);
                          
        % --- 2c. 计算 DTW 距离并转为权重 ---
        distance = dtw(live_sequence(:), map_sequence(:));
        
        %% [修改] 存储该粒子的距离
        all_distances(m) = distance;
        
        % 使用高斯核将距离转换为权重
        weights(m) = exp(-distance^2 / (2 * dtw_variance));
    end
    
    % 归一化权重
    % ( ... 您的归一化代码 ... )
    sum_weights = sum(weights);
    if sum_weights > 1e-15
        weights = weights / sum_weights;
    else
        weights = ones(M, 1) / M;
    end
    
    % --- 3. 重采样 (Resampling) ---
    % ( ... 您的重采样代码 ... )
    particles_out = zeros(M, 3);
    cdf = cumsum(weights);
    r_0 = rand() / M; 
    idx_j = 1;
    for m = 1:M
        U = r_0 + (m-1)/M;
        while U > cdf(idx_j)
            idx_j = idx_j + 1;
        end
        particles_out(m, :) = propagated_particles(idx_j, :);
    end
    
    % --- 4. 估计 (Estimation) ---
    % ( ... 您的估计代码 ... )
    best_guess = mean(particles_out, 1);
    best_guess(3) = atan2(mean(sin(particles_out(:, 3))), mean(cos(particles_out(:, 3))));
    
    %% [修改] 在函数末尾, 找出并返回最小的DTW距离
    d_min = min(all_distances);
end