clc;
clear;
close all;
%% 1. 仿真参数 (Parameters)
% [H-APF] M_init 对应 PDF 中的 M_init 
M_init = 2000;            % 初始粒子数量 
NUM_STEPS = 100;          % 仿真步数 (T) 
SEQUENCE_LEN = 25;        % 序列长度 (L) [cite: 19, 23]
% 地图尺寸
MAP_X_LEN = 50;
MAP_Y_LEN = 50;
% 噪声参数
SENSOR_NOISE_STD = 0.5;   % 传感器噪声 (sigma_e) [cite: 15, 23]
% [修改 1: 大幅提高DTW噪声标准差]
DTW_NOISE_STD = 25;     
% [修改 2: 调大粒子过程噪声]
process_noise.step_std = 0.3;         % (大于真实的 0.2)
process_noise.theta_std = deg2rad(12.5); % (大于真实的 10.0 度)
% [H-APF] ----------------------------------------------------
% [H-APF] 添加 Algorithm 1 和 3 所需的 H-APF 参数 
fprintf('加载 H-APF 自适应参数...\n');
APF.M_min = 500;          % 最小粒子数 (M_min) 
APF.M_max = 10000;        % 最大粒子数 (M_max) 
APF.K = 10;               % 秩统计量参数 (K) [cite: 23, 24]
APF.W = 20;               % chi-2 检验的窗口大小 (W) 
APF.p_l = 0.05;           % p-value 下阈值 (p_l) 
APF.p_h = 0.95;           % p-value 上阈值 (p_h) 
% [H-APF] ----------------------------------------------------
%% 2. 生成2D地磁地图 (Map Generation)
fprintf('生成非对称平滑地磁地图...\n');
% 1. 生成坐标网格 (X, Y grids)
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
% 2. 调用你的 Generator 生成 "白噪声"
Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
% 3. [!!] 对白噪声地图进行高斯平滑 [!!]
Mag = imgaussfilt(Mag_raw, 3.0); 
% 4. 将所有内容打包到 geo_map 结构体中
geo_map.X_grid = X;
geo_map.Y_grid = Y;
geo_map.Mag_map = Mag;
fprintf('地图加载完成. 尺寸: %d x %d\n', MAP_X_LEN, MAP_Y_LEN);
%% 3. 初始化 (Initialization)
fprintf('初始化仿真...\n');
% 真实状态: [x, y, theta]
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; 
% 粒子: [Nx3] 矩阵, 每一行是 [x, y, theta]
INIT_POS_STD = 5.0;  % 初始位置的标准差 (5个单位)
INIT_ANG_STD = 0.5;  % 初始角度的标准差 (0.5 弧度)
particles = zeros(M_init, 3); % [H-APF] 使用 M_init 
% 1. 在真实 X (true_state(1)) 附近初始化
particles(:, 1) = true_state(1) + randn(M_init, 1) * INIT_POS_STD;
% 2. 在真实 Y (true_state(2)) 附近初始化
particles(:, 2) = true_state(2) + randn(M_init, 1) * INIT_POS_STD;
% 3. 在真实 Theta (true_state(3)) 附近初始化
particles(:, 3) = true_state(3) + randn(M_init, 1) * INIT_ANG_STD;
% 历史记录 (用于绘图)
full_true_path_history = zeros(NUM_STEPS, 3); % 存储 [x, y, theta]
full_pdr_step_history = zeros(NUM_STEPS, 2);  % 存储 [step_len, d_theta]
true_path_history = zeros(NUM_STEPS, 2);      % 存储 [x, y] 用于绘图
estimated_path_history = zeros(NUM_STEPS, 2); % 存储 [x, y] 用于绘图
% 存储初始状态 (t=1)
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);
estimated_path_history(1, :) = true_path_history(1, :);
% [H-APF] ----------------------------------------------------
% [H-APF] 初始化 Algorithm 1 所需的变量 
current_M = M_init;                     % M_t, 初始为 M_1 
S_buffer = nan(APF.W, 1);               % 循环缓冲区 S 
M_history = zeros(NUM_STEPS, 1);        % 存储 M_t 历史
M_history(1) = current_M;
% [H-APF] ----------------------------------------------------
%% 4. 运行仿真 (Run Simulation)
fprintf('运行 %d 步仿真 (H-APF + DTW)...\n', NUM_STEPS);
h_waitbar = waitbar(0, '运行 H-APF 2D 粒子滤波+DTW...');
for t = 2:NUM_STEPS
   
    % --- 4a. 模拟真实运动 (PDR) ---
    [true_state, pdr_step] = get_next_step_random(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    
    % 存储完整的历史记录 (用于主滤波器)
    full_pdr_step_history(t, :) = pdr_step;   % 这是 u_t [cite: 7]
    full_true_path_history(t, :) = true_state; % 这是 X_t
    
    % --- [H-APF] 4b. 获取瞬时真实观测 (M-控制器需要) ---
    % (对应 Alg 1, line 10: 获取 y_t^real)
    % (对应 Eq 14: y_t = h(X_t) + e_t)
    y_t_real = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                      true_state(1), true_state(2), 'linear', 0) + ...
                      randn() * SENSOR_NOISE_STD; % [cite: 14, 15]
                  
    % --- [H-APF] 4c. H-APF M-控制器: 预测与统计 ---
    % (对应 Alg 1, block [b, c, d]) 
    
    % [b. M-控制器: 预测] (Alg 1, lines 12-16) 
    % P_t_bar = f(P_{t-1}, u_t, v_t)
    P_t_bar = propagate_particles(particles, pdr_step, process_noise); % [cite: 13]
    
    % [c. M-控制器: 生成统计量] (Alg 1, line 18) [cite: 18]
    a_t = Generate_AK_Statistic(P_t_bar, y_t_real, geo_map, ...
                              SENSOR_NOISE_STD^2, APF.K); % [cite: 18]
    
    % [d. M-控制器: 更新历史] (Alg 1, line 21) [cite: 21]
    % (t=2 是第一次迭代, 存入索引 1)
    S_buffer_idx = mod(t - 2, APF.W) + 1;
    S_buffer(S_buffer_idx) = a_t; % [cite: 21]
    
    % --- 4d. H-APF 主滤波器: 状态估计 ---
    % (对应 Alg 1, block [e]) 
    
    % i. 获取PDR历史 (U_t) (Alg 1, line 26) 
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    pdr_history_for_function = zeros(SEQUENCE_LEN, 2);
    actual_len = end_idx - start_idx + 1;
    pdr_history_for_function(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
    
    % ii. 生成 "Live" 传感器序列 (Y_t) (Alg 1, line 25) [cite: 25]
    live_sequence = zeros(1, SEQUENCE_LEN);
    path_segment = zeros(SEQUENCE_LEN, 3);
    path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
    for k = 1:SEQUENCE_LEN
        pos_x = path_segment(k, 1);
        pos_y = path_segment(k, 2);
        if pos_x == 0 && pos_y == 0 
            live_sequence(k) = 0;
        else
            live_sequence(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                                       pos_x, pos_y, 'linear', 0);
        end
    end
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD;
    
    % iii. 调用 *2D* 粒子滤波器 (Alg 1, line 27) 
    if t == 2 && ~exist('dtw', 'file')
        error('函数 "dtw" 未找到。此脚本需要 Signal Processing Toolbox。');
    end
    % [P_t^out, X_hat_t] = Particle_Filter_DTW_Step(...) 
    [particles_out, best_guess_state] = Particle_Filter_DTW_Step_2D(particles, live_sequence, ...
                                    pdr_history_for_function, geo_map, process_noise, DTW_NOISE_STD);
    
    % --- [H-APF] 4e. H-APF M-控制器: 自适应调整 ---
    % (对应 Alg 1, block [f]) [cite: 29]
    
    % 检查是否到了调整窗口 (Alg 1, line 30)
    % (t=2 是第1步, t=W+1 是第 W 步)
    if (t-1 >= APF.W) && (mod(t-1, APF.W) == 0) %
        % (Alg 1, line 31)
        M_new = Adapt_Particle_Count(current_M, S_buffer, APF.K, APF.W, ...
                                   APF.p_l, APF.p_h, APF.M_min, APF.M_max);
        
        % (Alg 1, line 32)
        particles_next = Adjust_Particle_Set(particles_out, M_new);
    else
        % (Alg 1, line 34, 35)
        M_new = current_M; % M_{t+1} = M_t
        particles_next = particles_out; % P_t = P_t^out
    end
    
    % 为下一次迭代更新状态
    particles = particles_next; % P_t
    current_M = M_new;          % M_{t+1}
    
    % --- 4f. 存储结果 (用于绘图) ---
    % (对应 Alg 1, block [g])
    true_path_history(t, :) = true_state(1:2); 
    estimated_path_history(t, :) = best_guess_state(1:2); % 存储 X_hat_t
    M_history(t) = current_M; % 存储 M_{t+1}
    
    waitbar(t/NUM_STEPS, h_waitbar, sprintf('H-APF (M=%d)...', current_M));
end
close(h_waitbar);
fprintf('仿真完成.\n');
%% 5. 绘图 (Plotting)
% [H-APF] 修改 figure 尺寸以容纳 3 个子图
figure('Position', [100, 100, 1600, 500]);
% --- 图 1: 2D 路径图 ---
subplot(1, 3, 1); % [H-APF] 修改为 1x3 布局
imagesc(geo_map.X_grid(1,:), geo_map.Y_grid(:,1), geo_map.Mag_map);
hold on;
axis xy; 
colormap('jet');
colorbar;
title('2D 路径跟踪 (H-APF+DTW)');
xlabel('X 位置');
ylabel('Y 位置');
plot(particles(:, 1), particles(:, 2), 'k.', 'MarkerSize', 2, 'DisplayName', '最终粒子云');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', '真实路径');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'H-APF 估计');
plot(true_path_history(1, 1), true_path_history(1, 2), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', '真实起点');
plot(estimated_path_history(1, 1), estimated_path_history(1, 2), 'gs', 'MarkerSize', 12, 'DisplayName', '估计起点');
plot(true_path_history(end, 1), true_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', '真实终点');
plot(estimated_path_history(end, 1), estimated_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', '估计终点');
legend('show', 'Location', 'best');
axis equal; 
axis([0 MAP_X_LEN 0 MAP_Y_LEN]);
% --- 图 2: 误差图 ---
subplot(1, 3, 2); % [H-APF] 修改为 1x3 布局
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('定位误差 (欧氏距离)');
xlabel('时间步 (Time Step)');
ylabel('误差 (Error in meters/units)');
grid on;
ylim([0, Inf]);
% --- [H-APF] 图 3: 自适应粒子数 M_t ---
subplot(1, 3, 3);
plot(1:NUM_STEPS, M_history, 'b-', 'LineWidth', 2);
title('自适应粒子数量 (M_t)');
xlabel('时间步 (Time Step)');
ylabel('粒子数量 (M)');
grid on;
line([1, NUM_STEPS], [APF.M_min, APF.M_min], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'M_min');
line([1, NUM_STEPS], [APF.M_max, APF.M_max], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'M_max');
legend('M_t', '界限', 'Location', 'best');
ylim([APF.M_min*0.8, APF.M_max*1.2]);
%% 6. [H-APF] H-APF 辅助函数
% (以下是 PDF 中 Algorithms 2, 3, 4 以及状态转移模型的实现)
% -----------------------------------------------------------------
%% [H-APF] 辅助函数 1: 状态转移 (Propagate Particles)
% 对应 PDF 中的状态转移模型 (Eq 10) 
% 注意: PDF Eq 10 中 theta 的更新公式似乎有误 (显示为 'cos(...)')。
%       我们假设其意图是标准的角度累加: theta_t = theta_{t-1} + delta_theta_t + v_theta_t
function P_out = propagate_particles(P_in, u_t, process_noise)
    M = size(P_in, 1);
    
    % 从 P_in 中解包状态 [cite: 6]
    X_t_minus_1 = P_in(:, 1);
    Y_t_minus_1 = P_in(:, 2);
    Theta_t_minus_1 = P_in(:, 3);
    
    % 解包控制输入 u_t [cite: 7]
    step_t = u_t(1);
    d_theta_t = u_t(2);
    
    % 生成过程噪声 v_t [cite: 11]
    v_step = randn(M, 1) * process_noise.step_std;   % v_step,t
    v_theta = randn(M, 1) * process_noise.theta_std; % v_theta,t
    
    % 应用状态转移方程 (Eq 10) 
    % (使用 PDF 中的 sin/cos 约定)
    theta_new_noisy = Theta_t_minus_1 + d_theta_t + v_theta;
    step_new_noisy = step_t + v_step;
    
    X_t = X_t_minus_1 + sin(theta_new_noisy) .* step_new_noisy; % 
    Y_t = Y_t_minus_1 + cos(theta_new_noisy) .* step_new_noisy; % 
    
    % (修正后的 theta 更新)
    Theta_t = mod(theta_new_noisy + pi, 2*pi) - pi; % 保持在 [-pi, pi]
    
    P_out = [X_t, Y_t, Theta_t];
end
% -----------------------------------------------------------------
%% [H-APF] 辅助函数 2: Algorithm 2 (Generate AK Statistic) 
function a_t = Generate_AK_Statistic(P_t_bar, y_t_real, geo_map, sensor_noise_var, K)
    % 输入: P_t_bar (预测的粒子) 
    %       y_t_real (真实瞬时观测) 
    %       geo_map (h 函数) [cite: 15, 24]
    %       sensor_noise_var (sigma_e^2) 
    %       K (统计量参数) 
    
    M_t = size(P_t_bar, 1); % 
    Y_fict = zeros(K, 1);   % 虚拟观测集 
    sensor_noise_std = sqrt(sensor_noise_var);
    
    for k = 1:K % 
        % a. 均匀选择粒子 
        j = randi(M_t); % 
        X_j = P_t_bar(j, :); % 
        
        % b. 从 p(y|X_j) 采样 
        % mu_j = h(X_j) 
        mu_j = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                      X_j(1), X_j(2), 'linear', 0); % [cite: 15, 24]
        
        % y_tilde ~ N(mu_j, sigma_e^2) 
        y_tilde_k = mu_j + randn() * sensor_noise_std; 
        Y_fict(k) = y_tilde_k; % 
    end
    
    % c. 计算秩统计量 
    a_t = sum(Y_fict < y_t_real); % 
end
% -----------------------------------------------------------------
%% [H-APF] 辅助函数 3: Algorithm 3 (Adapt_Particle_Count) 
function M_new = Adapt_Particle_Count(M_curr, S, K, W, p_l, p_h, M_min, M_max)
    % (此函数需要 Statistics and Machine Learning Toolbox 才能使用 chi2cdf)
    if ~exist('chi2cdf', 'file')
        warning('H-APF:MissingToolbox', '函数 "chi2cdf" 未找到 (需要 Statistics Toolbox)。跳过粒子数自适应。');
        M_new = M_curr;
        return;
    end
    
    % a. 计算观测频数 O_j 
    O = zeros(K + 1, 1); % 
    for j = 0:K % 
        O(j + 1) = sum(S == j); % 
    end
    
    % b. 计算期望频数 E_j 
    E = W / (K + 1); % 
    
    % c. 计算 chi^2 统计量 
    chi2_stat = sum((O - E).^2 / E); % 
    
    % d. 计算 p-value 
    p_value = 1 - chi2cdf(chi2_stat, K); % 自由度为 K 
    
    % e. 决策 
    if p_value <= p_l % 收敛性差 
        M_new = min(M_curr * 2, M_max); % 
    elseif p_value >= p_h % 收敛性好, 可能过剩 
        M_new = max(M_curr / 2, M_min); % 
    else % 可接受范围 
        M_new = M_curr; % 
    end
    
    M_new = round(M_new); % 
end
% -----------------------------------------------------------------
%% [H-APF] 辅助函数 4: Algorithm 4 (Adjust_Particle_Set) 
function P_out = Adjust_Particle_Set(P_in, M_new)
    M_curr = size(P_in, 1); % 
    
    if M_new == M_curr % 
        P_out = P_in; % 
        return;
        
    elseif M_new > M_curr % 增加粒子 
        N_add = M_new - M_curr; % 
        
        % 带替换随机采样 
        indices = randi(M_curr, N_add, 1); % 
        P_add = P_in(indices, :);
        
        % 添加抖动 (Jitter) 以避免粒子完全重合 
        % (抖动幅度可调, 这里设一个较小的值)
        jitter_pos = randn(N_add, 2) * 0.1; 
        jitter_ang = randn(N_add, 1) * 0.01;
        P_add(:, 1:2) = P_add(:, 1:2) + jitter_pos;
        P_add(:, 3) = P_add(:, 3) + jitter_ang;
        
        P_out = [P_in; P_add]; % 
        
    else % M_new < M_curr, 减少粒子 
        
        % 不带替换随机采样 
        indices = randperm(M_curr, M_new); % 
        P_out = P_in(indices, :); % 
    end
end