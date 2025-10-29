function MAIN_2D_PFDTW_GPU
%==========================================================================
% SCRIPT: MAIN_2D_PFDTW (Hybrid GPU/CPU-Parallel Single-File Version)
%
% DESCRIPTION:
%   Main simulation script for a 2D Adaptive Particle Filter (APF)
%   aided by PDR and geomagnetic matching using DTW.
%   This version is accelerated using a hybrid approach:
%     - GPU: Particle propagation, resampling, map storage.
%     - CPU (Parallel): DTW calculation (the main bottleneck).
%
% TO RUN:
%   1. Save this entire file as "MAIN_2D_PFDTW_GPU.m"
%   2. Type "MAIN_2D_PFDTW_GPU" in the MATLAB Command Window.
%
% DEPENDENCIES:
%   - Signal Processing Toolbox (for 'dtw' function)
%   - Parallel Computing Toolbox (for 'parfor' and 'parpool')
%   - A CUDA-enabled NVIDIA GPU (for 'gpuArray')
%==========================================================================
clc;
clear;
close all;

% [GPU/PARALLEL] 启动 CPU 并行池
if isempty(gcp('nocreate'))
    parpool; % <--- 改回：从 'threads' 改回默认的 'Processes'
end
fprintf('CPU Parallel Pool (Processes) started.\n');

%% 1. SIMULATION PARAMETERS
%==========================================================================
% (您的所有参数... 保持不变)
M_init = 100;            % Initial particle count
NUM_STEPS = 100;          % Number of simulation steps
SEQUENCE_LEN = 50;        % Sequence length for DTW
MAP_X_LEN = 50;           % Map width (X)
MAP_Y_LEN = 50;           % Map height (Y)
SENSOR_NOISE_STD = 0.5;   % Std dev of sensor noise (added to live sequence)
DTW_NOISE_STD = 60;      % Std dev for weight calculation (converts DTW dist to weight)
process_noise.step_std = 0.5;          % Process noise: std dev for step length
process_noise.theta_std = deg2rad(10); % Process noise: std dev for heading
pos_std = 5.0;           % Std dev for re-injected particle position
ang_std = 0.2;            % Std dev for re-injected particle angle
APF.M_min = 500;          % Minimum allowed particle count
APF.M_max = 10000000;      % Maximum allowed particle count
APF.DTW_THRESH_HIGH = 100; % DTW distance above which M is increased
APF.DTW_THRESH_LOW  = 10.0; % DTW distance below which M is decreased
%% 2. MAP GENERATION
%==========================================================================
fprintf('Generating geomagnetic map (and sending to GPU)...\n');
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
Mag = imgaussfilt(Mag_raw, 3.0); 

% Store map in a struct
% [GPU] 将地图数据发送到 GPU
geo_map.X_grid = gpuArray(X);
geo_map.Y_grid = gpuArray(Y);
geo_map.Mag_map = gpuArray(Mag);
fprintf('Map generation complete (on GPU).\n');
%% 3. INITIALIZATION
%==========================================================================
fprintf('Initializing simulation (on GPU)...\n');
% --- Initial True State ---
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; 
% --- Initial Particle Set ---
INIT_POS_STD = 5.0; 
INIT_ANG_STD = 0.5; 
% [GPU] 直接在 GPU 上创建粒子
particles = zeros(M_init, 3, 'gpuArray');
particles(:, 1) = true_state(1) + randn(M_init, 1, 'gpuArray') * INIT_POS_STD;
particles(:, 2) = true_state(2) + randn(M_init, 1, 'gpuArray') * INIT_POS_STD;
particles(:, 3) = true_state(3) + randn(M_init, 1, 'gpuArray') * INIT_ANG_STD;

% --- History Logs (for plotting) ---
% (历史记录保留在 CPU 上)
full_true_path_history = zeros(NUM_STEPS, 3); 
full_pdr_step_history = zeros(NUM_STEPS, 2);  
true_path_history = zeros(NUM_STEPS, 2);      
estimated_path_history = zeros(NUM_STEPS, 2); 
% --- Store Initial State (t=1) ---
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);
estimated_path_history(1, :) = true_path_history(1, :);
% --- Initialize Adaptive M Variables ---
current_M = M_init;                     
M_history = zeros(NUM_STEPS, 1);        
M_history(1) = current_M;

% [GPU] 将噪声参数也发送到 GPU
process_noise.step_std = gpuArray(process_noise.step_std);
process_noise.theta_std = gpuArray(process_noise.theta_std);
pos_std_gpu = gpuArray(pos_std);
ang_std_gpu = gpuArray(ang_std);
DTW_NOISE_STD_gpu = gpuArray(DTW_NOISE_STD);
SENSOR_NOISE_STD_gpu = gpuArray(SENSOR_NOISE_STD);

%% 4. RUN SIMULATION
%==========================================================================
fprintf('Running %d simulation steps (Hybrid DTW-APF)...\n', NUM_STEPS);
h_waitbar = waitbar(0, 'Running DTW-APF (GPU + CPU Parallel)...');
for t = 2:NUM_STEPS
   
    % --- 4a. Simulate True Motion (PDR) ---
    [true_state, pdr_step] = get_next_step_random(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step;
    full_true_path_history(t, :) = true_state;
    
    % --- 4b. Prepare Function Inputs ---
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    
    % Input 1: PDR History (U_t)
    pdr_history_cpu = zeros(SEQUENCE_LEN, 2);
    actual_len = end_idx - start_idx + 1;
    pdr_history_cpu(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
    pdr_history_for_function = gpuArray(pdr_history_cpu); % [GPU] 发送

    % Input 2: "Live" Sensor Sequence (Y_t)
    live_sequence = zeros(1, SEQUENCE_LEN, 'gpuArray'); % [GPU] 直接创建
    path_segment = zeros(SEQUENCE_LEN, 3); % (CPU, for loop)
    path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
    
    for k = 1:SEQUENCE_LEN
        pos_x = path_segment(k, 1);
        pos_y = path_segment(k, 2);
        
        if pos_x ~= 0 || pos_y ~= 0
            % [GPU] interp2 在 GPU 上运行
            live_sequence(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                                       pos_x, pos_y, 'linear', 0);
        end
    end
    % [GPU] 在 GPU 上添加噪声
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN, 'gpuArray') * SENSOR_NOISE_STD_gpu;
    
    % --- 4c. Call the 2D Particle Filter Step ---
    if t == 2 && ~exist('dtw', 'file')
        error('Function "dtw" not found. This script requires the Signal Processing Toolbox.');
    end
    
    % [GPU] 调用修改后的混合函数 (它在下面定义为本地函数)
    [particles_out, best_guess_state, dist] = Particle_Filter_DTW_Step_2D(particles, live_sequence, ...
                                    pdr_history_for_function, geo_map, process_noise, DTW_NOISE_STD_gpu);
    
    
    % --- 4d. APF Control: Re-injection and Adaptation ---
    
    % Particle Re-injection (prevents particle depletion)
    M_current_step = size(particles_out, 1);
    N_reset = round(M_current_step * 0.05); 
    indices_to_reset = randperm(M_current_step, N_reset);
    
    % [GPU] 重置在 GPU 上运行
    particles_out(indices_to_reset, 1) = best_guess_state(1) + randn(N_reset, 1, 'gpuArray') * pos_std_gpu;
    particles_out(indices_to_reset, 2) = best_guess_state(2) + randn(N_reset, 1, 'gpuArray') * pos_std_gpu;
    particles_out(indices_to_reset, 3) = best_guess_state(3) + randn(N_reset, 1, 'gpuArray') * ang_std_gpu;

    % M-Adaptation (adjust particle count)
    % [GPU] dist 是从 PF 函数返回的 CPU 标量, 直接使用
    M_new = Adapt_Particle_Count_DTW(current_M, dist, ...
                APF.DTW_THRESH_LOW, APF.DTW_THRESH_HIGH, ...
                APF.M_min, APF.M_max);
    
    % Apply the new particle count
    if M_new ~= current_M
        % [GPU] 调用修改后的 GPU 版本函数 (它在下面定义为本地函数)
        particles_next = Adjust_Particle_Set(particles_out, M_new);
    else
        particles_next = particles_out; % No change
    end
    
    % Update state for the next iteration
    particles = particles_next;
    current_M = size(particles, 1);
                                    
    % --- 4e. Store Results ---
    % [GPU] 从 GPU "gather" (收集) 回 CPU 以便存储
    true_path_history(t, :) = true_state(1:2);
    estimated_path_history(t, :) = gather(best_guess_state(1:2));
    M_history(t) = current_M;
    
    % --- 4f. Update Waitbar ---
    waitbar(t/NUM_STEPS, h_waitbar, sprintf('DTW-APF (M=%d, dist=%.1f) [GPU/CPU]', current_M, dist));
end
close(h_waitbar);
fprintf('Simulation complete.\n');
%% 5. PLOTTING
%==========================================================================
figure('Position', [100, 100, 1600, 500]);
% --- Plot 1: 2D Path ---
subplot(1, 3, 1); 
% [GPU] 从 GPU 收集数据以便绘图
imagesc(gather(geo_map.X_grid(1,:)), gather(geo_map.Y_grid(:,1)), gather(geo_map.Mag_map));
hold on;
axis xy; 
colormap('jet');
colorbar;
title('2D Path Tracking  - Random Path'); 
xlabel('X Position');
ylabel('Y Position');

% [GPU] 收集最终的粒子云
particles_cpu = gather(particles);
plot(particles_cpu(:, 1), particles_cpu(:, 2), 'k.', 'MarkerSize', 2, 'DisplayName', 'Final Particle Cloud');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', 'True Path');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'DTW-APF Estimate'); 
plot(true_path_history(1, 1), true_path_history(1, 2), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'True Start');
plot(estimated_path_history(1, 1), estimated_path_history(1, 2), 'gs', 'MarkerSize', 12, 'DisplayName', 'Est. Start');
plot(true_path_history(end, 1), true_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'True End');
plot(estimated_path_history(end, 1), estimated_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'Est. End');
legend('show', 'Location', 'best');
axis equal; 
axis([0 MAP_X_LEN 0 MAP_Y_LEN]);
% --- Plot 2: Error ---
subplot(1, 3, 2); 
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('Localization Error (Euclidean Distance)');
xlabel('Time Step');
ylabel('Error (in meters/units)');
grid on;
ylim([0, Inf]);
% --- Plot 3: Adaptive Particle Count (M_t) ---
subplot(1, 3, 3);
plot(1:NUM_STEPS, M_history, 'b-', 'LineWidth', 2);
title('Adaptive Particle Count (M_t)');
xlabel('Time Step');
ylabel('Particle Count (M)');
grid on;
line([1, NUM_STEPS], [APF.M_min, APF.M_min], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'M_min');
line([1, NUM_STEPS], [APF.M_max, APF.M_max], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'M_max');
legend('show', 'Location', 'best');
ylim([APF.M_min*0.8, APF.M_max*1.2]);

end % <--- [!!] 结束主函数 [!!]

%% 
%==========================================================================
% --- LOCAL HELPER FUNCTIONS ---
%==========================================================================

%% 1. 粒子滤波核心函数 (Hybrid GPU/CPU-Parallel)
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