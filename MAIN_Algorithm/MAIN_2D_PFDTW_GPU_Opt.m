%==========================================================================
% SCRIPT: MAIN_2D_PFDTW_GPU_V11_2 (DDTW-KLD-Denoise-GPF-HPC Arch)
%
%
% 依赖: 
%   - Parallel Computing Toolbox (parpool)
%   - Signal Processing Toolbox (dtw)
%   - Wavelet Toolbox (wavedec, wmaxlev, mad, wthresh, waverec)
%   - A CUDA-enabled NVIDIA GPU (gpuArray)
%==========================================================================
clc;
clear;
close all;
% [GPU/PARALLEL] 启动 CPU 并行池
if isempty(gcp('nocreate'))
    parpool;
end
fprintf('CPU Parallel Pool (Processes) started.\n');

% --- 0. GPU 和精度初始化 ---
fprintf('正在检查GPU...\n');
try
    gpu = gpuDevice(1);
    fprintf('在 %s (GPU) 上运行。\n', gpu.Name);
    disp('使用 "single" 精度 (float) 运行.');
    use_gpu = true;
    fp_precision = 'single';
catch ME
    warning('未找到GPU。将回退到CPU运行 (double 精度)。');
    use_gpu = false;
    fp_precision = 'double';
    gpuArray = @(x) x; 
    gather = @(x) x;
    randn = @(varargin) randn(varargin{:}, fp_precision);
    zeros = @(varargin) zeros(varargin{:}, fp_precision);
end

%% 1. SIMULATION PARAMETERS (V11: DDTW + Low Noise)
%==========================================================================
M_init = 500;            
NUM_STEPS = 100;         
SEQUENCE_LEN = 50;       
MAP_X_LEN = 50;          
MAP_Y_LEN = 50;           
SENSOR_NOISE_STD = 0.5;   


% 低过程噪声.
process_noise.step_std = 0.05;          
process_noise.theta_std = deg2rad(1.0); 

pos_std = 5.0;            % 局部注入 (Locked-on) 时的噪声
ang_std = 0.2;            % 局部注入 (Locked-on) 时的噪声

% 权重尺度 (用于柯西权重 1 / (1 + (dist/scale)^2))
DDTW_WEIGHT_SCALE = 20.0; % (V10.3 调谐值)

% --- KLD Sampling Parameters (基于 DDTW 距离) ---
KLD.M_min = 500;          
KLD.M_max = 100000;       
KLD.bin_size_xy = 0.5;    
KLD.bin_size_theta = deg2rad(10); 
KLD.epsilon = 0.05;       
KLD.delta = 0.01;         
%  "跟丢" 触发器 (基于 DDTW *距离*)
KLD.RECOVERY_DIST_THRESHOLD = 50.0; % (V10.3 调谐值)

%% 2. MAP GENERATION 
%==========================================================================
fprintf('生成地磁地图并发送到GPU (%s)...\n', fp_precision);
[X_cpu, Y_cpu] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag_raw_cpu = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
Mag_cpu = imgaussfilt(Mag_raw_cpu, 3.0); 

geo_map.X_grid = gpuArray(cast(X_cpu, fp_precision));
geo_map.Y_grid = gpuArray(cast(Y_cpu, fp_precision));
geo_map.Mag_map = gpuArray(cast(Mag_cpu, fp_precision));
fprintf('地图已在GPU上 (%s 精度)。\n', fp_precision);

%% 3. INITIALIZATION 
%==========================================================================
fprintf('在GPU上初始化粒子 (%s)...\n', fp_precision);
true_state_cpu = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; 
INIT_POS_STD = 5.0; 
INIT_ANG_STD = 0.5; 

particles = zeros(M_init, 3, fp_precision, 'gpuArray'); 
particles(:, 1) = cast(true_state_cpu(1), fp_precision) + randn(M_init, 1, fp_precision, 'gpuArray') * INIT_POS_STD;
particles(:, 2) = cast(true_state_cpu(2), fp_precision) + randn(M_init, 1, fp_precision, 'gpuArray') * INIT_POS_STD;
particles(:, 3) = cast(true_state_cpu(3), fp_precision) + randn(M_init, 1, fp_precision, 'gpuArray') * INIT_ANG_STD;

% --- History Logs ---
full_true_path_history = zeros(NUM_STEPS, 3); 
full_pdr_step_history = zeros(NUM_STEPS, 2);  
true_path_history = zeros(NUM_STEPS, 2);      
estimated_path_history = zeros(NUM_STEPS, 2); 
full_true_path_history(1, :) = true_state_cpu;
true_path_history(1, :) = true_state_cpu(1:2);
estimated_path_history(1, :) = true_path_history(1, :);
current_M = M_init;                     
M_history = zeros(NUM_STEPS, 1);        
M_history(1) = current_M;

% --- 将噪声参数发送到 GPU ---
process_noise_gpu.step_std = gpuArray(cast(process_noise.step_std, fp_precision));
process_noise_gpu.theta_std = gpuArray(cast(process_noise.theta_std, fp_precision));
pos_std_gpu = gpuArray(cast(pos_std, fp_precision));
ang_std_gpu = gpuArray(cast(ang_std, fp_precision));
SENSOR_NOISE_STD_gpu = gpuArray(cast(SENSOR_NOISE_STD, fp_precision));
DDTW_WEIGHT_SCALE_gpu = gpuArray(cast(DDTW_WEIGHT_SCALE, fp_precision)); % [V11] 新增

%% 4. RUN SIMULATION
%==========================================================================
fprintf('运行 %d 模拟步 (DDTW-KLD-GPF V11.2)...\n', NUM_STEPS);
h_waitbar = waitbar(0, 'Running DDTW-KLD-GPF (V11.2)...');
for t = 2:NUM_STEPS
   
    % --- 4a. Simulate True Motion (CPU) ---
    [true_state_cpu, pdr_step_cpu] = get_next_step_random(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step_cpu;
    full_true_path_history(t, :) = true_state_cpu;
    
    % --- 4b. Prepare Function Inputs  ---
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    
    pdr_history_cpu = zeros(SEQUENCE_LEN, 2);
    actual_len = end_idx - start_idx + 1;
    pdr_history_cpu(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
    pdr_history_for_function = gpuArray(cast(pdr_history_cpu, fp_precision));
    
    live_sequence_gpu = zeros(1, SEQUENCE_LEN, fp_precision, 'gpuArray');
    path_segment_cpu = zeros(SEQUENCE_LEN, 3);
    path_segment_cpu(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
    
    for k = 1:SEQUENCE_LEN
        pos_x = path_segment_cpu(k, 1);
        pos_y = path_segment_cpu(k, 2);
        
        if pos_x ~= 0 || pos_y ~= 0
            live_sequence_gpu(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                                           pos_x, pos_y, 'linear', 0);
        end
    end
    live_sequence_gpu = live_sequence_gpu + randn(1, SEQUENCE_LEN, fp_precision, 'gpuArray') * SENSOR_NOISE_STD_gpu;
    
    % --- 小波降噪 ---
    live_sequence_noisy_cpu = gather(live_sequence_gpu);
    live_sequence_denoised_cpu = Denoise_Visushrink_CPU(live_sequence_noisy_cpu, 'db4');
    live_sequence_for_function = gpuArray(live_sequence_denoised_cpu);
    
    % --- 4c. Call the 2D Particle Filter Step ---
    if t == 2 && ~exist('dtw', 'file') % [V11] 检查 dtw
        error('Function "dtw" (for DDTW) not found. This script requires the Signal Processing Toolbox.');
    end
    if t == 2 && ~exist('wavedec', 'file')
        warning('Function "wavedec" not found. Wavelet Denoising will be skipped. Requires Wavelet Toolbox.');
        live_sequence_for_function = live_sequence_gpu;
    end
    
    % 调用 (DDTW) 核心函数
    [particles_out, best_guess_state_gpu, d_min] = Particle_Filter_Step_V11_DDTW(particles, live_sequence_for_function, ...
                                    pdr_history_for_function, geo_map, process_noise_gpu, DDTW_WEIGHT_SCALE_gpu);
    
    
    % --- 4d. APF Control (基于 DDTW 距离) ---
    M_current_step = size(particles_out, 1);
    best_guess_state_cpu = gather(best_guess_state_gpu); 
    
    if d_min > KLD.RECOVERY_DIST_THRESHOLD % (d_min > 50.0)
        % --- 状态 1: 我们跟丢了 (Diverged) ---
        
        N_reset = round(M_current_step * 0.80); 
        indices_to_reset = randperm(M_current_step, N_reset);
        
        map_x_max_cpu = gather(geo_map.X_grid(1,end));
        map_y_max_cpu = gather(geo_map.Y_grid(end,1));
        
        random_x = 1 + rand(N_reset, 1, fp_precision, 'gpuArray') * (map_x_max_cpu - 1);
        random_y = 1 + rand(N_reset, 1, fp_precision, 'gpuArray') * (map_y_max_cpu - 1);
        random_theta = rand(N_reset, 1, fp_precision, 'gpuArray') * 2 * pi;
    
        particles_out(indices_to_reset, 1) = random_x;
        particles_out(indices_to_reset, 2) = random_y;
        particles_out(indices_to_reset, 3) = random_theta;
        
    else
        % --- 状态 2: 我们锁定了 (Locked-on) ---
        N_reset = round(M_current_step * 0.05); 
        indices_to_reset = randperm(M_current_step, N_reset);
        
        particles_out(indices_to_reset, 1) = best_guess_state_cpu(1) + randn(N_reset, 1, fp_precision, 'gpuArray') * pos_std_gpu;
        particles_out(indices_to_reset, 2) = best_guess_state_cpu(2) + randn(N_reset, 1, fp_precision, 'gpuArray') * pos_std_gpu;
        particles_out(indices_to_reset, 3) = best_guess_state_cpu(3) + randn(N_reset, 1, fp_precision, 'gpuArray') * ang_std_gpu;
    end
    
    
    % --- KLD 自适应 ( 基于 DDTW 距离) ---
    particles_for_kld = gather(particles_out);
    M_new = Adapt_Particle_Count_KLD_V11(particles_for_kld, d_min, KLD);
    
    if M_new ~= current_M
        particles_next = Adjust_Particle_Set_GPU(particles_out, M_new);
    else
        particles_next = particles_out; 
    end
    
    particles = particles_next;
    current_M = size(particles, 1);
                                    
    % --- 4e. Store Results ---
    true_path_history(t, :) = true_state_cpu(1:2);
    estimated_path_history(t, :) = best_guess_state_cpu(1:2);
    M_history(t) = current_M;
    
    % --- 4f. Update Waitbar ---
    waitbar(t/NUM_STEPS, h_waitbar, sprintf('DDTW-KLD-GPF (M=%d, dist=%.2f) [V11.2]', current_M, d_min));
end
close(h_waitbar);
fprintf('Simulation complete.\n');
%% 5. PLOTTING 
%==========================================================================
figure('Position', [100, 100, 1600, 500]);
% --- Plot 1: 2D Path ---
subplot(1, 3, 1); 
imagesc(X_cpu(1, :), Y_cpu(:, 1), Mag_cpu); 
hold on;
axis xy; 
colormap('jet');
colorbar;
title('2D Path Tracking - (DDTW-KLD-GPF V11.2)'); 
xlabel('X Position');
ylabel('Y Position');
particles_cpu = gather(particles); 
plot(particles_cpu(:, 1), particles_cpu(:, 2), 'k.', 'MarkerSize', 2, 'DisplayName', 'Final Particle Cloud');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', 'True Path');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'DDTW-APF Estimate'); 
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
title('Adaptive Particle Count (KLD)');
xlabel('Time Step');
ylabel('Particle Count (M)');
grid on;
line([1, NUM_STEPS], [KLD.M_min, KLD.M_min], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'M_min');
line([1, NUM_STEPS], [KLD.M_max, KLD.M_max], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'M_max');
legend('show', 'Location', 'best');
ylim([KLD.M_min*0.8, KLD.M_max*1.2]);
 % <--- [!!] 结束主函数 [!!]
%% 
%==========================================================================
% --- LOCAL HELPER FUNCTIONS  ---
%==========================================================================
%% 1. 粒子滤波核心函数 (DDTW + 低噪声 + HPC 架构)
function [particles_out, best_guess, d_min] = Particle_Filter_Step_V11_DDTW(particles_in, live_sequence, ...
                                            pdr_history, geo_map, process_noise, ddtw_weight_scale)
% 描述: (DDTW 架构)
%   1. [HPC] 路径重建在 *GPU* 上矢量化执行.
%   2. [parfor] CPU 并行循环 *只* 运行 DDTW 权重计算.
%   3. [Algo] 使用 DDTW (dtw(diff(...))) 算法.
%
    
    M = size(particles_in, 1);
    L = size(pdr_history, 1);
    fp_precision = underlyingType(particles_in); 
    
    % --- 1. 传播 (Prediction/Propagation) [ON GPU] ---
    % process_noise.theta_std 现在是 1.0 deg (非常小)
    last_pdr_step = pdr_history(end, :); 
    
    step_noise = randn(M, 1, fp_precision, 'gpuArray') * process_noise.step_std;
    theta_noise = randn(M, 1, fp_precision, 'gpuArray') * process_noise.theta_std;
    
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
    
    % --- 2. [HPC 优化] 加权 (Weighting) ---
    
    % --- 2a. 矢量化路径重建 [ON GPU] ---
    all_particle_paths_x_gpu = zeros(M, L, fp_precision, 'gpuArray');
    all_particle_paths_y_gpu = zeros(M, L, fp_precision, 'gpuArray');
    all_particle_thetas_gpu = zeros(M, L, fp_precision, 'gpuArray');
    
    all_particle_paths_x_gpu(:, L) = propagated_particles(:, 1);
    all_particle_paths_y_gpu(:, L) = propagated_particles(:, 2);
    all_particle_thetas_gpu(:, L) = propagated_particles(:, 3);
    
    pdr_history_cpu = gather(pdr_history); 
    
    for k = L-1:-1:1
        pdr_step_len = pdr_history_cpu(k+1, 1);
        pdr_d_theta = pdr_history_cpu(k+1, 2);
        
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
    
    % --- 2b. 矢量化地图采样 [ON GPU] ---
    all_map_sequences_gpu = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                                  all_particle_paths_x_gpu, all_particle_paths_y_gpu, 'linear', 0);
                              
    % --- 2c. [GPU -> CPU] 收集数据 ---
    all_map_sequences_cpu = gather(all_map_sequences_gpu);
    live_sequence_cpu = gather(live_sequence); % 已降噪
    ddtw_scale_cpu = gather(ddtw_weight_scale); % 权重尺度
    
    % --- 2d. 并行 DDTW 计算 [ON CPU] ---
    weights_cpu = zeros(M, 1, fp_precision);
    all_distances_cpu = zeros(M, 1, fp_precision);
    
    parfor m = 1:M
        map_sequence_raw = all_map_sequences_cpu(m, :);
                          
        % --- [算法: DDTW] ---
        distance = cast(Inf, fp_precision);
        weight = cast(0, fp_precision);
        
        % 确保序列足够长以计算 diff
        if length(map_sequence_raw) >= 2 && length(live_sequence_cpu) >= 2
            map_sequence_deriv = diff(map_sequence_raw);
            live_sequence_deriv = diff(live_sequence_cpu(:)'); % 确保是行向量
            
            if ~isempty(map_sequence_deriv) && ~isempty(live_sequence_deriv)
                % 1. 计算 DDTW 距离
                distance = dtw(map_sequence_deriv, live_sequence_deriv);
                
                % 2. 使用鲁棒的柯西权重 
                % 1 / (1 + (dist/scale)^2)
                weight = 1.0 / (1.0 + (distance / ddtw_scale_cpu)^2);
            end
        end
        
        all_distances_cpu(m) = distance;
        weights_cpu(m) = weight;
    end
    
    % d_min 是 KLD 和 "跟丢" 逻辑所需要的
    d_min = min(all_distances_cpu); 
    weights = gpuArray(weights_cpu);
    
    % --- 3. 重采样 (Resampling) [ON GPU] ---
    sum_weights = sum(weights);
    if sum_weights > 1e-15
        weights = weights / sum_weights;
    else
        weights = ones(M, 1, fp_precision, 'gpuArray') / M;
    end
    
    cdf = cumsum(weights);
    r_0 = rand(1, fp_precision, 'gpuArray') / M;
    U = r_0 + (gpuArray.colon(0, M-1)') / M;
    
    [~, indices] = histc(U, [0; cdf]);
    indices = min(indices, M); 
    indices(indices == 0) = 1; 
    
    particles_out = propagated_particles(indices, :); 
    
    % --- 4. 估计 (Estimation) [ON GPU] ---
    best_guess = mean(particles_out, 1);
    best_guess(3) = atan2(mean(sin(particles_out(:, 3))), mean(cos(particles_out(:, 3))));
end
%% 2. 粒子数调整函数 (GPU Version)
function P_out = Adjust_Particle_Set_GPU(P_in, M_new_double)

    
    fp_precision = underlyingType(P_in); 
    
    % M_curr 和 M_new 必须是 *整数* 才能用于 randperm 和索引.
    M_curr_double = size(P_in, 1);
    M_curr = int64(M_curr_double); % 强制转换为 int64
    M_new = int64(M_new_double);  % 强制转换为 int64
    
    if M_new == M_curr
        P_out = P_in;
        return;
        
    elseif M_new > M_curr % 增加粒子
        N_add = M_new - M_curr; % int64 运算
        
        % randi 接受 int64 作为输入 (但 M_curr_double 也可以)
        indices = randi(M_curr_double, N_add, 1, 'gpuArray'); 
        P_add = P_in(indices, :);
        
        jitter_pos = randn(N_add, 2, fp_precision, 'gpuArray') * 0.1; 
        jitter_ang = randn(N_add, 1, fp_precision, 'gpuArray') * 0.01;
        P_add(:, 1:2) = P_add(:, 1:2) + jitter_pos;
        P_add(:, 3) = P_add(:, 3) + jitter_ang;
        
        P_out = [P_in; P_add]; 
        
    else % M_new < M_curr, 减少粒子
        
        p_cpu = randperm(M_curr_double); % 1. 在 CPU 上安全运行
        indices_cpu = p_cpu(1:M_new);    % 2. 在 CPU 上索引 (M_new 是 int64)
        indices_gpu = gpuArray(indices_cpu); % 3. 仅发送索引回 GPU
        
        P_out = P_in(indices_gpu, :);    % 4. 在 GPU 上索引
    end
end
%% 3. KLD 自适应计数函数
function M_new = Adapt_Particle_Count_KLD_V11(particles_cpu, d_min, KLD)
% 描述: (基于 *DDTW 距离*)
    
    % 状态 1: 检查是否完全跟丢 (d_min > 50.0)
    if d_min > KLD.RECOVERY_DIST_THRESHOLD
        M_new = KLD.M_max;
        M_new = round(M_new);
        return;
    end
    
    % 状态 2: 未跟丢, 使用 KLD 统计计算
    
    % 1. 将状态空间离散化 (分箱)
    x_bins = ceil(particles_cpu(:, 1) / KLD.bin_size_xy);
    y_bins = ceil(particles_cpu(:, 2) / KLD.bin_size_xy);
    theta_normalized = particles_cpu(:, 3) + pi; % 映射到 [0, 2*pi]
    theta_bins = ceil(theta_normalized / KLD.bin_size_theta);
    
    % 2. 统计非空分箱的数量 'k'
    all_bins = [x_bins, y_bins, theta_bins];
    unique_bins = unique(all_bins, 'rows', 'stable');
    k = size(unique_bins, 1); 
    
    if k <= 1
        M_new = KLD.M_min;
    else
        % 3. 使用 KLD 的 chi-squared 近似公式计算 M_new
        z_delta = norminv(1 - KLD.delta); 
        term1 = (k - 1) / (2 * KLD.epsilon);
        term2 = 1 - (2 / (9 * (k - 1)));
        term3 = sqrt(2 / (9 * (k - 1))) * z_delta;
        M_new = term1 * (term2 + term3)^3;
    end
    
    % 4. 钳位 (Clamp) 到 [M_min, M_max] 范围
    M_new = max(KLD.M_min, min(KLD.M_max, M_new));
    M_new = round(M_new);
end

%% 4. 运动模型函数
function [true_state, pdr_step] = get_next_step_random(prev_state, MAP_X_LEN, MAP_Y_LEN)
    step_len = 0.8 + randn()*0.1;
    d_theta = deg2rad(2.0) + randn()*deg2rad(5);
    
    state = prev_state;
    state(3) = state(3) + d_theta; 
    state(1) = state(1) + sin(state(3)) * step_len;
    state(2) = state(2) + cos(state(3)) * step_len;
    
    if state(1) < 1 || state(1) > MAP_X_LEN || state(2) < 1 || state(2) > MAP_Y_LEN
        state = prev_state; 
        state(3) = prev_state(3) + deg2rad(90); 
    end
    
    true_state = state;
    pdr_step = [step_len, d_theta];
end
%% 5. 地图生成器函数
function Mag_raw = Geometric_Map_Generator(seed, map_size)
    rng(seed);
    [X, Y] = meshgrid(1:map_size(2), 1:map_size(1));
    Mag_raw = (X / map_size(2)) + (Y / map_size(1)) + randn(map_size) * 0.5;
    Mag_raw = (Mag_raw - min(Mag_raw(:))) / (max(Mag_raw(:)) - min(Mag_raw(:)));
end

%% 6. 小波降噪函数
function denoised_data = Denoise_Visushrink_CPU(data, wavelet)
%DENOISE_VISUSHRINK (CPU Version, 兼容 'single' 和 'double') 
%   在 CPU 上执行小波降噪 (需要 Wavelet Toolbox)
%
    n = length(data);
    fp_precision = class(data); 
    
    % --- 步骤 1: 小波分解 ---
    wlevel = wmaxlev(n, wavelet);
    if wlevel < 1
        denoised_data = data;
        return;
    end
    [c, l] = wavedec(data, wlevel, wavelet); 

    % --- 步骤 2: 计算"通用门槛" ---
    cD1_len = l(end-1);
    if cD1_len == 0
        denoised_data = data;
        return;
    end
    cD1_start_index = sum(l(1:end-2)) + 1;
    cD1_end_index = cD1_start_index + cD1_len - 1;
    
    cD1 = c(cD1_start_index : cD1_end_index);
    
    % 2b. 计算门槛
    sigma_est = mad(cD1, 1) / cast(0.6745, fp_precision); 
    threshold = sigma_est * sqrt(2 * log(cast(n, fp_precision)));

    % --- 步骤 3: 应用软阈值 ---
    cA = c(1:l(1)); 
    detail_coeffs = c(l(1)+1 : end); 
    
    denoised_details = wthresh(detail_coeffs, 's', threshold);
    
    denoised_c = [cA, denoised_details];

    % --- 步骤 4: 小波重构 ---
    denoised_data_temp = waverec(denoised_c, l, wavelet);
    
    if length(denoised_data_temp) > n
        denoised_data = denoised_data_temp(1:n);
    else
        % 确保返回长度与输入一致
        denoised_data = [denoised_data_temp, zeros(1, n-length(denoised_data_temp), fp_precision)];
        denoised_data = denoised_data(1:n);
    end
end
