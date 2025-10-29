%% ========================================================================
%  FILE: Adjust_Particle_Set.m
%  (一个鲁棒的粒子数调整函数)
%  ========================================================================
function particles_next = Adjust_Particle_Set(particles_curr, M_new)
    % 摘要:
    %   将当前的粒子集 'particles_curr' 调整为 'M_new' 个粒子。
    %   - 如果 M_new > M_curr, 随机复制现有粒子。
    %   - 如果 M_new < M_curr, 随机丢弃多余粒子。
    %
    % 输入:
    %   particles_curr: M_curr x 3 矩阵 [x, y, theta]
    %   M_new:          目标粒子数
    %
    % 输出:
    %   particles_next: M_new x 3 矩阵 [x, y, theta]

    % --- 1. [安全措施] 确保所有输入都是整数 ---
    M_new = round(M_new);
    M_curr = size(particles_curr, 1);

    % --- 2. 比较粒子数 ---
    if M_new == M_curr
        % 数量不变
        particles_next = particles_curr;
        return;
        
    elseif M_new > M_curr
        % --- 3. 增加粒子 ---
        N_add = M_new - M_curr; % N_add 现在保证是整数
        
        % 随机选择 N_add 个现有粒子进行复制
        % [修正] 这就是你出错的第 14 行
        indices = randi(M_curr, N_add, 1); 
        
        particles_to_add = particles_curr(indices, :);
        
        % (可选) 给新粒子增加一点点噪声, 避免完全重叠
        % noise_pos = (rand(N_add, 2) - 0.5) * 0.1; % +/- 5cm
        % noise_ang = (rand(N_add, 1) - 0.5) * deg2rad(1); % +/- 0.5 deg
        % particles_to_add(:, 1:2) = particles_to_add(:, 1:2) + noise_pos;
        % particles_to_add(:, 3) = particles_to_add(:, 3) + noise_ang;
        
        particles_next = [particles_curr; particles_to_add];
        
    else % M_new < M_curr
        % --- 4. 减少粒子 (丢弃) ---
        % M_new 此时 < M_curr
        
        % 随机选择 M_new 个粒子保留下来
        indices_to_keep = randperm(M_curr, M_new);
        
        particles_next = particles_curr(indices_to_keep, :);
    end
end