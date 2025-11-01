%==========================================================================
% SCRIPT: MAIN_2D_PFDTW_GPU (Hybrid GPU/CPU-Parallel Robustness Test)
%
% DESCRIPTION:
%   执行鲁棒性分析实验 (GPU/CPU 加速版)。
%   本脚本循环运行一个简化的 MAIN 脚本, 每次使用不同的 DTW_NOISE_STD 值,
%   并记录 t=10 时的估计位置, 以测试算法对该参数的敏感性。
%
% 依赖:
%   - Signal Processing Toolbox (for 'dtw' function)
%   - Parallel Computing Toolbox (for 'parfor' and 'parpool')
%   - A CUDA-enabled NVIDIA GPU (for 'gpuArray')
%==========================================================================
clc;
clear;
close all;
% [GPU/PARALLEL] 启动 CPU 并行池
if isempty(gcp('nocreate'))
    parpool; % <--- 使用 'Processes' 模式
end
fprintf('CPU Parallel Pool (Processes) started.\n');
%% 1. 实验设置 (ROBUSTNESS TEST)
%==========================================================================
% 定义要测试的 DTW_NOISE_STD (sigma) 范围
std_values_to_test = 20:10:90; % [20, 30, ..., 90]

% 固定的仿真步数
NUM_STEPS_FOR_TEST = 10; % <--- [修改] 运行到 t=10

% 结果存储 (在 CPU 上)
results_dist_from_origin = zeros(length(std_values_to_test), 1);
true_dist_from_origin = 0; % 用于存储 t=10 时的真实距离

fprintf('开始执行鲁棒性分析 (共 %d 次运行)...\n', length(std_values_to_test));

%% 2. MAP GENERATION
%==========================================================================
fprintf('生成共享地磁图 (并发送到 GPU)...\n');
MAP_X_LEN = 50;
MAP_Y_LEN = 50;
[X_cpu, Y_cpu] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag_raw_cpu = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
Mag_cpu = imgaussfilt(Mag_raw_cpu, 3.0); 

% [GPU] 将地图数据发送到 GPU
geo_map.X_grid = gpuArray(X_cpu);
geo_map.Y_grid = gpuArray(Y_cpu);
geo_map.Mag_map = gpuArray(Mag_cpu);
fprintf('Map generation complete (on GPU).\n');

%% 3. 循环执行仿真
%==========================================================================
h_waitbar_main = waitbar(0, '正在测试鲁棒性 (GPU/CPU)...');

for k_run = 1:length(std_values_to_test)
    
    waitbar(k_run / length(std_values_to_test), h_waitbar_main, ...
        sprintf('运行 %d/%d (sigma = %.1f)', k_run, length(std_values_to_test), std_values_to_test(k_run)));
    
    % --- 3a. 重置环境 (关键!) ---
    
    % [关键!] 重置随机种子, 确保真实路径/噪声/初始化在每次运行时都相同
    rng(123); 
    
    % --- 3b. 加载参数 ---
    
    % [被测参数]
    DTW_NOISE_STD_val = std_values_to_test(k_run); % <--- 这是我们循环的变量
    DTW_NOISE_STD_gpu = gpuArray(DTW_NOISE_STD_val); % 发送到 GPU
    
    M_init = 500; % <--- [修改] 从 100 增加到 500, 与您之前的实验一致
    NUM_STEPS = NUM_STEPS_FOR_TEST; 
    SEQUENCE_LEN = 50;
    SENSOR_NOISE_STD = 0.5;
    process_noise.step_std = gpuArray(0.5);
    process_noise.theta_std = gpuArray(deg2rad(10));
    pos_std_gpu = gpuArray(5.0);
    ang_std_gpu = gpuArray(0.2);
    SENSOR_NOISE_STD_gpu = gpuArray(SENSOR_NOISE_STD);
    
    APF.M_min = 500;
    APF.M_max = 10000000;
    APF.DTW_THRESH_HIGH = 100;
    APF.DTW_THRESH_LOW  = 10.0;
    
    % --- 3c. 运行初始化 (%% 3) ---
    true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)];
    INIT_POS_STD = 5.0;
    INIT_ANG_STD = 0.5;
    
    % [GPU] 直接在 GPU 上创建粒子
    particles = zeros(M_init, 3, 'gpuArray');
    particles(:, 1) = true_state(1) + randn(M_init, 1, 'gpuArray') * INIT_POS_STD;
    particles(:, 2) = true_state(2) + randn(M_init, 1, 'gpuArray') * INIT_POS_STD;
    particles(:, 3) = true_state(3) + randn(M_init, 1, 'gpuArray') * INIT_ANG_STD;
    
    % (历史记录保留在 CPU 上)
    full_true_path_history = zeros(NUM_STEPS, 3);
    full_pdr_step_history = zeros(NUM_STEPS, 2);
    full_true_path_history(1, :) = true_state;
    current_M = M_init;
    
    % --- 3d. 运行仿真 (%% 4) ---
    for t = 2:NUM_STEPS
        [true_state, pdr_step] = get_next_step_random(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
        full_pdr_step_history(t, :) = pdr_step;
        full_true_path_history(t, :) = true_state;
        
        start_idx = max(1, t - SEQUENCE_LEN + 1);
        end_idx = t;
        
        pdr_history_cpu = zeros(SEQUENCE_LEN, 2);
        actual_len = end_idx - start_idx + 1;
        pdr_history_cpu(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
        pdr_history_for_function = gpuArray(pdr_history_cpu);
        
        live_sequence = zeros(1, SEQUENCE_LEN, 'gpuArray');
        path_segment = zeros(SEQUENCE_LEN, 3);
        path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
        
        for k = 1:SEQUENCE_LEN
            pos_x = path_segment(k, 1);
            pos_y = path_segment(k, 2);
            if pos_x ~= 0 || pos_y ~= 0
                live_sequence(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, pos_x, pos_y, 'linear', 0);
            end
        end
        live_sequence = live_sequence + randn(1, SEQUENCE_LEN, 'gpuArray') * SENSOR_NOISE_STD_gpu;
        
        % [GPU] 调用混合函数
        [particles_out, best_guess_state, dist] = Particle_Filter_DTW_Step_2D(particles, live_sequence, ...
                                        pdr_history_for_function, geo_map, process_noise, DTW_NOISE_STD_gpu);
        
        M_current_step = size(particles_out, 1);
        N_reset = round(M_current_step * 0.05);
        indices_to_reset = randperm(M_current_step, N_reset);
        
        % [GPU] 重置在 GPU 上运行
        particles_out(indices_to_reset, 1) = best_guess_state(1) + randn(N_reset, 1, 'gpuArray') * pos_std_gpu;
        particles_out(indices_to_reset, 2) = best_guess_state(2) + randn(N_reset, 1, 'gpuArray') * pos_std_gpu;
        particles_out(indices_to_reset, 3) = best_guess_state(3) + randn(N_reset, 1, 'gpuArray') * ang_std_gpu;
        
        % [CPU] dist 是从 PF 函数返回的 CPU 标量, 直接使用
        M_new = Adapt_Particle_Count_DTW(current_M, dist, ...
                    APF.DTW_THRESH_LOW, APF.DTW_THRESH_HIGH, ...
                    APF.M_min, APF.M_max);
        
        if M_new ~= current_M
            particles_next = Adjust_Particle_Set(particles_out, M_new);
        else
            particles_next = particles_out;
        end
        
        particles = particles_next;
        current_M = size(particles, 1);
    end
    
    % --- 3e. 记录 t=10 (最后一步) 的结果 ---
    % [GPU -> CPU] 从 GPU 收集最终的估计值
    best_guess_cpu = gather(best_guess_state);
    x_est = best_guess_cpu(1);
    y_est = best_guess_cpu(2);
    results_dist_from_origin(k_run) = sqrt(x_est^2 + y_est^2);
    
    % 记录一次 t=10 时的真实距离 (在所有运行中都相同)
    if k_run == 1
        x_true = true_state(1);
        y_true = true_state(2);
        true_dist_from_origin = sqrt(x_true^2 + y_true^2);
    end
    
    
    clearvars -except std_values_to_test results_dist_from_origin true_dist_from_origin k_run geo_map h_waitbar_main NUM_STEPS_FOR_TEST MAP_X_LEN MAP_Y_LEN
end
close(h_waitbar_main);
fprintf('鲁棒性分析完成。\n');

%% --- 4. 绘图 (PLOTTING) ---
%==========================================================================
figure('Position', [200, 200, 900, 500]);
plot(std_values_to_test, results_dist_from_origin, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% hold on;
% line([min(std_values_to_test), max(std_values_to_test)], ...
%      [true_dist_from_origin, true_dist_from_origin], ...
%      'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
 
title('算法对 \sigma_{dtw} 参数的鲁棒性分析 (t=10, GPU/CPU版)');
xlabel('DTW 噪声标准差 (\sigma_{dtw})');
ylabel('距离原点的欧氏距离 (sqrt(x_{est}² + y_{est}²))');
legend('估计距离 (Est. Distance)','Location', 'best');
grid on;
set(gca, 'FontSize', 12);

%% 
%==========================================================================
% --- LOCAL HELPER FUNCTIONS (与您提供的代码相同) ---
%==========================================================================

function [particles_out, best_guess, d_min] = Particle_Filter_DTW_Step_2D(particles_in, live_sequence, ...
                                            pdr_history, geo_map, process_noise, DTW_NOISE_STD)
% 描述: 执行一步混合加速 (GPU + CPU parfor) 的粒子滤波
%
% [GPU] 输入:
%   - particles_in (gpuArray)
%   - live_sequence (gpuArray)
%   - pdr_history (gpuArray)
%   - geo_map (struct of gpuArrays)
%   - process_noise (struct of gpuArrays)
%   - DTW_NOISE_STD (gpuArray)
% [GPU] 输出:
%   - particles_out (gpuArray)
%   - best_guess (gpuArray)
%   - d_min (scalar on CPU)
%
    
    M = size(particles_in, 1);
    L = size(pdr_history, 1);
    
    % --- 1. 传播 (Prediction/Propagation) ---
    % [GPU] 此部分完全在 GPU 上矢量化运行 (非常快)
    last_pdr_step = pdr_history(end, :);
    
    step_noise = randn(M, 1, 'gpuArray') * process_noise.step_std;
    theta_noise = randn(M, 1, 'gpuArray') * process_noise.theta_std;
    
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
    
    propagated_particles = [new_x, new_y, new_theta]; % (仍在 GPU 上)
    
    % --- 2. 加权 (Weighting) ---
    % [CPU-PARFOR] 这是瓶颈。 'dtw' 无法在 GPU 运行.
    % 我们将此循环移至 CPU 并行运行.
    
    % [GPU -> CPU] 从 GPU 收集 (gather) 循环所需的数据
    propagated_particles_cpu = gather(propagated_particles);
    live_sequence_cpu = gather(live_sequence);
    pdr_history_cpu = gather(pdr_history);
    geo_map_cpu.X_grid = gather(geo_map.X_grid);
    geo_map_cpu.Y_grid = gather(geo_map.Y_grid);
    geo_map_cpu.Mag_map = gather(geo_map.Mag_map);
    dtw_variance_cpu = gather(DTW_NOISE_STD^2);
    
    % (为 parfor 预分配 CPU 数组)
    weights_cpu = zeros(M, 1);
    all_distances_cpu = zeros(M, 1);
    
    % [CPU-PARFOR] 在所有 CPU 核心上并行运行 for 循环
    parfor m = 1:M
        % --- 2a. 重建粒子路径 (Path Reconstruction) ---
        particle_path = zeros(L, 2); % (在 CPU worker 上)
        
        current_pos = propagated_particles_cpu(m, 1:2);
        current_theta = propagated_particles_cpu(m, 3);
        particle_path(end, :) = current_pos;
        
        for k = L-1:-1:1
            pdr_step_len = pdr_history_cpu(k+1, 1);
            pdr_d_theta = pdr_history_cpu(k+1, 2);
            
            prev_theta = current_theta - pdr_d_theta;
            prev_x = current_pos(1) - sin(current_theta) * pdr_step_len;
            prev_y = current_pos(2) - cos(current_theta) * pdr_step_len;
            
            particle_path(k, :) = [prev_x, prev_y];
            current_pos = [prev_x, prev_y];
            current_theta = prev_theta;
        end
        
        % --- 2b. 生成地图序列 ---
        map_sequence = interp2(geo_map_cpu.X_grid, geo_map_cpu.Y_grid, geo_map_cpu.Mag_map, ...
                              particle_path(:, 1), particle_path(:, 2), 'linear', 0);
                          
        % --- 2c. 计算 DTW 距离并转为权重 ---
        distance = dtw(live_sequence_cpu(:), map_sequence(:));
        
        % 存储结果 (parfor 需要这样做)
        all_distances_cpu(m) = distance;
        weights_cpu(m) = exp(-distance^2 / (2 * dtw_variance_cpu));
    end
    
    % [CPU] 找到最小距离 (在主脚本中需要它)
    d_min = min(all_distances_cpu); 
    
    % [CPU -> GPU] 将结果权重发送回 GPU
    weights = gpuArray(weights_cpu);
    
    % --- 3. 重采样 (Resampling) ---
    % [GPU] 归一化在 GPU 上运行
    sum_weights = sum(weights);
    if sum_weights > 1e-15
        weights = weights / sum_weights;
    else
        weights = ones(M, 1, 'gpuArray') / M;
    end
    
    % [GPU] 在 GPU 上进行重采样 (使用快速的矢量化方法)
    cdf = cumsum(weights);
    r_0 = rand(1, 'gpuArray') / M;
    U = r_0 + (gpuArray.colon(0, M-1)') / M; % (在 GPU 上创建 U 向量)
    
    % 使用 'discretize' (比 for 循环快得多)
    [~, indices] = histc(U, [0; cdf]);
    indices = min(indices, M); % (确保索引在范围内)
    
    particles_out = propagated_particles(indices, :); % (在 GPU 上索引)
    
    % --- 4. 估计 (Estimation) ---
    % [GPU] 均值在 GPU 上计算
    best_guess = mean(particles_out, 1);
    best_guess(3) = atan2(mean(sin(particles_out(:, 3))), mean(cos(particles_out(:, 3))));
end
%% 2. 粒子数调整函数 (GPU Version)
function P_out = Adjust_Particle_Set(P_in, M_new)
% 描述: 在 GPU 上增加或减少粒子集
% [GPU] P_in 是一个 gpuArray
% [GPU] P_out 将是一个 gpuArray
    
    M_curr = size(P_in, 1);
    
    if M_new == M_curr
        P_out = P_in;
        return;
        
    elseif M_new > M_curr % 增加粒子
        N_add = M_new - M_curr;
        
        % [GPU] 在 GPU 上随机采样
        indices = randi(M_curr, N_add, 1, 'gpuArray'); 
        P_add = P_in(indices, :);
        
        % [GPU] 在 GPU 上添加抖动 (Jitter)
        jitter_pos = randn(N_add, 2, 'gpuArray') * 0.1; 
        jitter_ang = randn(N_add, 1, 'gpuArray') * 0.01;
        P_add(:, 1:2) = P_add(:, 1:2) + jitter_pos;
        P_add(:, 3) = P_add(:, 3) + jitter_ang;
        
        P_out = [P_in; P_add]; % [GPU] 在 GPU 上拼接
        
    else % M_new < M_curr, 减少粒子
        
        % [GPU] 在 GPU 上随机采样
        indices = randperm(M_curr, M_new, 'gpuArray'); 
        P_out = P_in(indices, :);
    end
end
%% 3. 自适应计数函数 (CPU Version)
function M_new = Adapt_Particle_Count_DTW(M_curr, dist, Thresh_Low, Thresh_High, M_min, M_max)
    % (此函数在 CPU 上运行，因为它只处理标量)
    if dist > Thresh_High 
        % 误差太大 (匹配很差)，增加粒子
        M_new = min(M_curr * 2, M_max);
    elseif dist < Thresh_Low 
        % 误差很小 (匹配很好)，减少粒子
        M_new = max(M_curr / 2, M_min);
    else
        % 误差在可接受范围
        M_new = M_curr;
    end
    
    M_new = round(M_new);
end
%% 4. [STUB] 运动模型函数
function [true_state, pdr_step] = get_next_step_random(prev_state, MAP_X_LEN, MAP_Y_LEN)
    % (这是一个占位符 STUB 函数, 请使用您自己的版本)
    % 简单的“反弹”逻辑
    step_len = 0.8 + randn()*0.1;
    d_theta = deg2rad(2.0) + randn()*deg2rad(5);
    
    state = prev_state;
    state(3) = state(3) + d_theta; % 更新角度
    state(1) = state(1) + sin(state(3)) * step_len;
    state(2) = state(2) + cos(state(3)) * step_len;
    
    % 边界检查
    if state(1) < 1 || state(1) > MAP_X_LEN || state(2) < 1 || state(2) > MAP_Y_LEN
        state = prev_state; % 撞墙, 重置
        state(3) = prev_state(3) + deg2rad(90); % 转弯
    end
    
    true_state = state;
    pdr_step = [step_len, d_theta];
end
%% 5. [STUB] 地图生成器函数
function Mag_raw = Geometric_Map_Generator(seed, map_size)
    % (这是一个占位符 STUB 函数, 请使用您自己的版本)
    % 生成简单的梯度 + 噪声地图
    rng(seed);
    [X, Y] = meshgrid(1:map_size(2), 1:map_size(1));
    Mag_raw = (X / map_size(2)) + (Y / map_size(1)) + randn(map_size) * 0.5;
    Mag_raw = (Mag_raw - min(Mag_raw(:))) / (max(Mag_raw(:)) - min(Mag_raw(:)));
end