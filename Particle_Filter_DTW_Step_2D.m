function [particles_out, best_guess] = Particle_Filter_DTW_Step_2D(particles_in, live_sequence, ...
                                            pdr_history, geo_map, process_noise, DTW_NOISE_STD)
% 摘要: 执行一步粒子滤波 (传播、加权、重采样)
% 权重: 基于DTW（动态时间规整）
% 注意: 此函数需要 Signal Processing Toolbox (用于 'dtw' 函数)
    
    M = size(particles_in, 1);
    L = size(pdr_history, 1);
    
    % --- 1. 传播 (Prediction/Propagation) ---
    % 我们使用 *最后一步* 的PDR来传播所有粒子
    last_pdr_step = pdr_history(end, :);
    
    % 添加 *过程噪声* (您的 "修改 2")
    step_noise = randn(M, 1) * process_noise.step_std;
    theta_noise = randn(M, 1) * process_noise.theta_std;
    
    steps = last_pdr_step(1) + step_noise;
    d_thetas = last_pdr_step(2) + theta_noise;
    
    % 矢量化传播
    old_x = particles_in(:, 1);
    old_y = particles_in(:, 2);
    old_theta = particles_in(:, 3);
    
    new_theta = old_theta + d_thetas;
    new_x = old_x + sin(new_theta) .* steps;
    new_y = old_y + cos(new_theta) .* steps;
    
    % 边界检查
    new_x = max(1, min(geo_map.X_grid(1,end), new_x));
    new_y = max(1, min(geo_map.Y_grid(end,1), new_y));
    
    propagated_particles = [new_x, new_y, new_theta];

    % --- 2. 加权 (Weighting) ---
    weights = zeros(M, 1);
    dtw_variance = DTW_NOISE_STD^2;
    
    for m = 1:M
        % --- 2a. 重建粒子路径 (Path Reconstruction) ---
        % 对于每个粒子, 我们必须 "倒带" 它的历史来构建它的路径
        particle_path = zeros(L, 2); % 存储 [x, y]
        
        % 路径的 *最后* 一点是我们刚刚 *传播* 到的点
        current_pos = propagated_particles(m, 1:2);
        current_theta = propagated_particles(m, 3);
        particle_path(end, :) = current_pos;
        
        % 循环 "倒带"
        for k = L-1:-1:1
            % 获取此步骤的 PDR (没有噪声!)
            pdr_step_len = pdr_history(k+1, 1);
            pdr_d_theta = pdr_history(k+1, 2);
            
            % 反转运动
            prev_theta = current_theta - pdr_d_theta;
            prev_x = current_pos(1) - sin(current_theta) * pdr_step_len;
            prev_y = current_pos(2) - cos(current_theta) * pdr_step_len;
            
            % 存储倒带的路径点
            particle_path(k, :) = [prev_x, prev_y];
            
            % 更新 "倒带" 状态
            current_pos = [prev_x, prev_y];
            current_theta = prev_theta;
        end
        
        % --- 2b. 生成地图序列 ---
        % 从地图中采样该粒子路径对应的地磁值
        map_sequence = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                              particle_path(:, 1), particle_path(:, 2), 'linear', 0);
                          
        % --- 2c. 计算 DTW 距离并转为权重 ---
        % (需要 Signal Processing Toolbox)
        distance = dtw(live_sequence(:), map_sequence(:));
        
        % 使用高斯核将距离转换为权重
        weights(m) = exp(-distance^2 / (2 * dtw_variance));
    end
    
    % 归一化权重
    sum_weights = sum(weights);
    if sum_weights > 1e-15
        weights = weights / sum_weights;
    else
        % 粒子全部失效 (权重为0) -> 重置
        weights = ones(M, 1) / M;
    end
    
    % --- 3. 重采样 (Resampling) ---
    % 使用标准 "Low Variance" / "Systematic" 重采样
    particles_out = zeros(M, 3);
    cdf = cumsum(weights);
    r_0 = rand() / M; % 初始随机数
    
    idx_j = 1;
    for m = 1:M
        U = r_0 + (m-1)/M;
        while U > cdf(idx_j)
            idx_j = idx_j + 1;
        end
        particles_out(m, :) = propagated_particles(idx_j, :);
    end
    
    % --- 4. 估计 (Estimation) ---
    % 最佳估计是重采样后粒子的均值
    best_guess = mean(particles_out, 1);
    
    % (处理角度均值，确保在 -pi 到 pi 之间)
    best_guess(3) = atan2(mean(sin(particles_out(:, 3))), mean(cos(particles_out(:, 3))));
end