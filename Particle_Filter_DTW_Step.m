function [new_particles, best_guess] = Particle_Filter_DTW_Step_ign(old_particles, live_sequence, ...
                                    geo_map, window_radius, process_noise, dtw_noise_std)
% 这就是你们方案的核心: 结合了DTW的粒子滤波
%
%   old_particles:   [Nx1] 上一步的粒子
%   live_sequence:   [1xK] "真实" 的地磁序列 (来自传感器)
%   geo_map:         [1xM] 完整的地磁地图
%   window_radius:   窗口半径 (例如 2)
%   process_noise:   运动噪声
%   dtw_noise_std:   DTW 距离的噪声方差 (用于高斯函数)
%
    N = length(old_particles);
    map_len = length(geo_map);
    window_len = 2*window_radius + 1;
    
    % --- 1. PREDICT (预测) ---
    % "抖动" 所有粒子 (使用简化的 -1, 0, +1 运动模型)
    noise = round(randn(N, 1) * process_noise);
    predicted_particles = old_particles + noise;
    
    % 边界检查: 确保粒子不会掉出 "可测量" 的地图范围
    predicted_particles(predicted_particles <= window_radius) = window_radius + 1;
    predicted_particles(predicted_particles >= (map_len - window_radius)) = map_len - window_radius - 1;

    % --- 2. UPDATE (更新权重) ---
    weights = zeros(N, 1); % 存储每个粒子的权重
    
    for i = 1:N
        % a. 获取第 i 个粒子的预测位置
        particle_loc = predicted_particles(i);
        
        % b. 从 "地图" 中提取这个粒子 "应该" 看到的序列
        map_sequence_indices = (particle_loc - window_radius) : (particle_loc + window_radius);
        map_sequence = geo_map(map_sequence_indices);
        

        dist = simple_global_dtw(live_sequence, map_sequence);
        

        weights(i) = exp( -(dist^2) / (2 * dtw_noise_std^2) );
    end
    

    if sum(weights) == 0
        weights = ones(N, 1); % 如果所有粒子都跑丢了，重置
    end
    weights = weights / sum(weights);

    % --- 3. RESAMPLE (重采样) ---
    % "杀死" 低权重的粒子, "复制" 高权重的粒子
    new_particles = zeros(N, 1);
    c = cumsum(weights);
    u = (0:N-1)/N + rand()/N; 
    
    i_ptr = 1;
    for j = 1:N
        while u(j) > c(i_ptr)
            i_ptr = i_ptr + 1;
        end
        new_particles(j) = predicted_particles(i_ptr);
    end

    % --- 4. ESTIMATE (估计) ---
    % 我们的最佳猜测是所有新粒子的平均位置
    best_guess = mean(new_particles);
end
