clc;
clear;
close all;
%% 1. 仿真参数 (Parameters)
M_init = 2000;            % 初始粒子数量 (M_init)
NUM_STEPS = 100;          % 仿真步数 (T)
SEQUENCE_LEN = 50;        %[cite_start]% 序列长度 (L) [cite: 19]
% 地图尺寸
MAP_X_LEN = 50;
MAP_Y_LEN = 50;
% 噪声参数
SENSOR_NOISE_STD = 0.5;  % [cite_start]% 传感器噪声 [cite: 15]
% [修改 1: 大幅提高DTW噪声标准差]
DTW_NOISE_STD = 100;     
% [修改 2: 调大粒子过程噪声]
process_noise.step_std = 0.8;         
process_noise.theta_std = deg2rad(15); 
 pos_std = 10.0; % 5个单位的标准差
 ang_std = 0.3; % 0.5 弧度的标准差
% [APF-DTW] --------------------------------------------------
% [APF-DTW] 基于 DTW 距离的自适应参数
fprintf('加载 [DTW-based] APF 自适应参数...\n');
APF.M_min = 500;          % 最小粒子数 (M_min)
APF.M_max = 1000000;        % 最大粒子数 (M_max)

% [新增] 基于 DTW 距离的阈值
% (!! 注意: 这两个阈值是示例值, 您需要根据您的 'dist' 输出进行调试 !!)
APF.DTW_THRESH_HIGH = 25.0; % 如果 dist 高于此值, 增加 M
APF.DTW_THRESH_LOW  = 10.0;  % 如果 dist 低于此值, 减少 M
% [APF-DTW] --------------------------------------------------
%% 2. 生成2D地磁地图 (Map Generation)
fprintf('生成非对称平滑地磁地图...\n');
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
Mag = imgaussfilt(Mag_raw, 3.0); 
geo_map.X_grid = X;
geo_map.Y_grid = Y;
geo_map.Mag_map = Mag;
fprintf('地图加载完成. 尺寸: %d x %d\n', MAP_X_LEN, MAP_Y_LEN);
%% 3. 初始化 (Initialization)
fprintf('初始化仿真...\n');
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; 
INIT_POS_STD = 5.0; 
INIT_ANG_STD = 0.5; 
particles = zeros(M_init, 3);
particles(:, 1) = true_state(1) + randn(M_init, 1) * INIT_POS_STD;
particles(:, 2) = true_state(2) + randn(M_init, 1) * INIT_POS_STD;
particles(:, 3) = true_state(3) + randn(M_init, 1) * INIT_ANG_STD;
% 历史记录 (用于绘图)
full_true_path_history = zeros(NUM_STEPS, 3); 
full_pdr_step_history = zeros(NUM_STEPS, 2);  
true_path_history = zeros(NUM_STEPS, 2);      
estimated_path_history = zeros(NUM_STEPS, 2); 
% 存储初始状态 (t=1)
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);
estimated_path_history(1, :) = true_path_history(1, :);
% [APF-DTW] --------------------------------------------------
% [APF-DTW] 初始化自适应 M 所需的变量
current_M = M_init;                     % M_t
M_history = zeros(NUM_STEPS, 1);        % 存储 M_t 历史
M_history(1) = current_M;
% [APF-DTW] --------------------------------------------------
%% 4. 运行仿真 (Run Simulation)
fprintf('运行 %d 步仿真 (DTW-APF)...\n', NUM_STEPS);
h_waitbar = waitbar(0, '运行 DTW-APF 2D 粒子滤波...');
for t = 2:NUM_STEPS
   
    % --- 4a. 模拟真实运动 (PDR) ---
    [true_state, pdr_step] = get_next_step_random(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step;   
    full_true_path_history(t, :) = true_state; 
    
    % --- 4b. 准备函数输入 ---
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    pdr_history_for_function = zeros(SEQUENCE_LEN, 2); % U_t
    actual_len = end_idx - start_idx + 1;
    pdr_history_for_function(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
    
    live_sequence = zeros(1, SEQUENCE_LEN); % Y_t
    path_segment = zeros(SEQUENCE_LEN, 3);
    path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
    for k = 1:SEQUENCE_LEN
        pos_x = path_segment(k, 1);
        pos_y = path_segment(k, 2);
        if pos_x == 0 && pos_y == 0 
            live_sequence(k) = 0;
        else
            live_sequence(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                                       pos_x, pos_y, 'linear', 0); %[cite_start]% h(X_t) [cite: 15]
        end
    end
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD; %[cite_start]% + e_t [cite: 14]
    
    % --- 4c. 调用 *2D* 粒子滤波器 ---
    if t == 2 && ~exist('dtw', 'file')
        error('函数 "dtw" 未找到。此脚本需要 Signal Processing Toolbox。');
    end
    
    % [修改] 调用您修改后的函数，捕获 dist (最小DTW距离)
    [particles_out, best_guess_state, dist] = Particle_Filter_DTW_Step_2D(particles, live_sequence, ...
                                    pdr_history_for_function, geo_map, process_noise, DTW_NOISE_STD);
    
    % --- [APF-DTW] 4d. M-控制器: 基于 DTW 的自适应调整 ---
    M_current_step = size(particles_out, 1);
    N_reset = round(M_current_step * 0.05); % 5% 的粒子
    
    % 随机选择 5% 的粒子索引
    indices_to_reset = randperm(M_current_step, N_reset);
    
    % 将这些粒子重置到当前最佳估计的"附近" (或者全图随机)
    % (这里我们用 "在最佳估计附近" 作为示例)
    % 
    % 
    particles_out(indices_to_reset, 1) = best_guess_state(1) + randn(N_reset, 1) * pos_std;
    particles_out(indices_to_reset, 2) = best_guess_state(2) + randn(N_reset, 1) * pos_std;
    particles_out(indices_to_reset, 3) = best_guess_state(3) + randn(N_reset, 1) * ang_std;

    % particles(indices_to_reset, 1) = rand(N_reset, 1) * MAP_X_LEN;
    % particles(indices_to_reset, 2) = rand(N_reset, 1) * MAP_Y_LEN;
    % particles(indices_to_reset, 3) = rand(N_reset, 1) * 2 * pi - pi;
    % --- [新增结束] ---
    % [新增] 调用新的自适应函数 (见 %% 6.)
    M_new = Adapt_Particle_Count_DTW(current_M, dist, ...
                APF.DTW_THRESH_LOW, APF.DTW_THRESH_HIGH, ...
                APF.M_min, APF.M_max);
    
    % [新增] 仅在 M 确实改变时才调用 Algorithm 4
    if M_new ~= current_M
        particles_next = Adjust_Particle_Set(particles_out, M_new);
    else
        particles_next = particles_out; %
    end
    
    % [新增] 为下一次迭代更新状态
    particles = particles_next; 
    current_M = size(particles, 1); % M_{t+1}
                                    
    % --- 4e. 存储结果 (用于绘图) ---
    true_path_history(t, :) = true_state(1:2); 
    estimated_path_history(t, :) = best_guess_state(1:2); % 存储 X_hat_t
    M_history(t) = current_M; % 存储 M_{t+1}
    
    waitbar(t/NUM_STEPS, h_waitbar, sprintf('DTW-APF (M=%d, dist=%.1f)', current_M, dist));
end
close(h_waitbar);
fprintf('仿真完成.\n');
%% 5. 绘图 (Plotting)
figure('Position', [100, 100, 1600, 500]);
% --- 图 1: 2D 路径图 ---
subplot(1, 3, 1); 
imagesc(geo_map.X_grid(1,:), geo_map.Y_grid(:,1), geo_map.Mag_map);
hold on;
axis xy; 
colormap('jet');
colorbar;
title('2D 路径跟踪 (DTW-APF) - 随机路径'); 
xlabel('X 位置');
ylabel('Y 位置');
plot(particles(:, 1), particles(:, 2), 'k.', 'MarkerSize', 2, 'DisplayName', '最终粒子云');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', '真实路径');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'DTW-APF 估计'); 
plot(true_path_history(1, 1), true_path_history(1, 2), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', '真实起点');
plot(estimated_path_history(1, 1), estimated_path_history(1, 2), 'gs', 'MarkerSize', 12, 'DisplayName', '估计起点');
plot(true_path_history(end, 1), true_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', '真实终点');
plot(estimated_path_history(end, 1), estimated_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', '估计终点');
legend('show', 'Location', 'best');
axis equal; 
axis([0 MAP_X_LEN 0 MAP_Y_LEN]);
% --- 图 2: 误差图 ---
subplot(1, 3, 2); 
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('定位误差 (欧氏距离)');
xlabel('时间步 (Time Step)');
ylabel('误差 (Error in meters/units)');
grid on;
ylim([0, Inf]);
% --- [APF-DTW] 图 3: 自适应粒子数 M_t ---
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
%% 6. [APF-DTW] 辅助函数
% -----------------------------------------------------------------
%% [APF-DTW] 辅助函数 1: [新] 基于 DTW 的 M 调整
% (此函数替代了 PDF 中的 Algorithm 3)
function M_new = Adapt_Particle_Count_DTW(M_curr, dist, Thresh_Low, Thresh_High, M_min, M_max)
    %
    
    % 这是您建议的新逻辑:
    if dist > Thresh_High 
        % 误差太大 (匹配很差)，增加粒子
        M_new = min(M_curr * 2, M_max); %
    elseif dist < Thresh_Low 
        % 误差很小 (匹配很好)，减少粒子
        M_new = max(M_curr / 2, M_min); %
    else
        % 误差在可接受范围
        M_new = M_curr; %
    end
    
    M_new = round(M_new); %
end
% -----------------------------------------------------------------
%% [APF-DTW] 辅助函数 2: Algorithm 4 (Adjust_Particle_Set)
% (此函数来自 PDF, 用于增加或减少粒子)
function P_out = Adjust_Particle_Set(P_in, M_new)
    %
    M_curr = size(P_in, 1); %
    
    if M_new == M_curr %
        P_out = P_in; %
        return;
        
    elseif M_new > M_curr % 增加粒子
        N_add = M_new - M_curr; %
        
        % 带替换随机采样
        indices = randi(M_curr, N_add, 1); 
        P_add = P_in(indices, :);
        
        % 添加抖动 (Jitter)
        jitter_pos = randn(N_add, 2) * 0.1; 
        jitter_ang = randn(N_add, 1) * 0.01;
        P_add(:, 1:2) = P_add(:, 1:2) + jitter_pos;
        P_add(:, 3) = P_add(:, 3) + jitter_ang;
        
        P_out = [P_in; P_add]; %
        
    else % M_new < M_curr, 减少粒子
        
        % 不带替换随机采样
        indices = randperm(M_curr, M_new); 
        P_out = P_in(indices, :); %
    end
end