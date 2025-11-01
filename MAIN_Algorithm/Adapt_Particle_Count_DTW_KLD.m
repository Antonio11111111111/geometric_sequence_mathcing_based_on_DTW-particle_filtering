%% ========================================================================
%  FILE: Adapt_Particle_Count_DTW_KLD.m
%  (v5 - 修正了非整数 M_new 的错误)
%  ========================================================================
function M_new = Adapt_Particle_Count_DTW_KLD(particles_out, d_min, bin_size_xy, ...
                                    map_size_x, map_size_y, ...
                                    kld_epsilon, kld_delta, ...
                                    kld_n_min, kld_n_max, ...
                                    DTW_THRESH_HIGH)
    % ... (前面 1-3 步, 计算 k 的代码保持不变) ...
    current_M = size(particles_out, 1);
    
    % --- 1. 初始化 2D 分箱网格 ---
    num_bins_x = ceil(map_size_x / bin_size_xy);
    num_bins_y = ceil(map_size_y / bin_size_xy);
    bin_grid = false(num_bins_y, num_bins_x);
    k = 0; 
    
    % --- 2. 对所有粒子进行分箱 ---
    for i = 1:current_M
        x = particles_out(i, 1);
        y = particles_out(i, 2);
        bin_x = min(num_bins_x, max(1, floor(x / bin_size_xy) + 1));
        bin_y = min(num_bins_y, max(1, floor(y / bin_size_xy) + 1));
        if ~bin_grid(bin_y, bin_x)
            bin_grid(bin_y, bin_x) = true;
            k = k + 1;
        end
    end
    
    % --- 4. 调用 KLD 计算 (基础估算) ---
    M_kld = calculate_n_kld(k, kld_epsilon, kld_delta);
    
    % --- 5. 鲁棒性检查: "观测质量" ---
    if (d_min > DTW_THRESH_HIGH) && (M_kld < kld_n_max)
        % *** 发生锁死 ***
        
        % [修正] 确保 M_new 永远是整数
        % M_new = current_M * 1.5; % <-- 错误行
        M_new = ceil(current_M * 1.3); % <-- 正确行 (使用 ceil 或 round)
        
    else
        % *** 正常情况 ***
        M_new = M_kld;
    end

    % --- 6. 应用边界 ---
    M_new = max(kld_n_min, M_new);
    M_new = min(kld_n_max, M_new);
    
    % [安全措施] 确保最终输出是整数
    M_new = round(M_new);
end