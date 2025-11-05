%==========================================================================
% SCRIPT: MAIN_2D_PFDTW_KLD_GPU (已修复所有类型)
%
% 描述: 
%   此脚本已修改为在GPU上管理粒子和地图数据。
%   整个流程已统一为 'single' 精度以确保GPU类型兼容性。
%
% 依赖:
%   - Parallel Computing Toolbox (用于 parfor)
%   - Wavelet Toolbox (用于 Denoise_Visushrink)
%   - [你必须自己提供] Geometric_Map_Generator.m
%   - [你必须自己提供] Get_Next_Step_2D.m
%   - [你必须自己提供] Adapt_Particle_Count_DTW_KLD.m
%==========================================================================
clc;
clear;
close all;

% --- 0. GPU 和并行初始化 ---
% parpool; % <-- 启动并行池, parfor 需要它

fprintf('正在检查GPU...\n');
try
    gpu = gpuDevice(1); % 尝试获取第一个GPU
    fprintf('在 %s (GPU) 上运行。\n', gpu.Name);
catch ME
    warning('未找到GPU。将回退到CPU运行。');
    gpuArray = @(x) x; 
    gather = @(x) x;
    % 确保 randn 和 zeros 在CPU模式下也使用 'single'
    randn = @(varargin) randn(varargin{:}, 'single');
    zeros = @(varargin) zeros(varargin{:}, 'single');
end

%% 1. SIMULATION PARAMETERS (无变化)
%==========================================================================
M_init = 500;            
NUM_STEPS = 100;         
SEQUENCE_LEN = 100;       
MAP_X_LEN = 50;          
MAP_Y_LEN = 50;           
SENSOR_NOISE_STD = 0.5;   
DTW_NOISE_STD = 30;      
process_noise.step_std = 0.05;          
process_noise.theta_std = deg2rad(1); 
pos_std = 5.0;           
ang_std = 0.2;            
KLD.n_min = 500;          
KLD.n_max = 100000;       
KLD.bin_size_xy = 0.5;    
KLD.epsilon = 0.05;       
KLD.delta = 0.01;         
KLD.DTW_THRESH_HIGH = 63.9; 
sigma = 0.05;
eps = 0.5;

APFRATE = 0.08;

%% 2. MAP GENERATION (发送到GPU, 使用 'single')
%==========================================================================
fprintf('生成地磁地图并发送到GPU...\n');
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
Mag = imgaussfilt(Mag_raw, 3.0); 

% *** [GPU 类型修复] ***
% 将地图数据强制为 'single' 类型发送到GPU
geo_map.X_grid = gpuArray(single(X));
geo_map.Y_grid = gpuArray(single(Y));
geo_map.Mag_map = gpuArray(single(Mag));
fprintf('地图已在GPU上 (single 精度)。\n');

%% 3. INITIALIZATION (在GPU上创建 'single' 粒子)
%==========================================================================
fprintf('在GPU上初始化粒子...\n');
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; 
INIT_POS_STD = 5.0; 
INIT_ANG_STD = 0.5; 

% *** [GPU 类型修复] ***
% 直接在GPU上创建 'single' 类型的粒子
particles = zeros(M_init, 3, 'single', 'gpuArray'); % 在GPU上创建
particles(:, 1) = single(true_state(1)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_POS_STD;
particles(:, 2) = single(true_state(2)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_POS_STD;
particles(:, 3) = single(true_state(3)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_ANG_STD;

% --- History Logs (日志仍在CPU上) ---
% (CPU日志使用默认的 double 精度是OK的)
full_true_path_history = zeros(NUM_STEPS, 3); 
full_pdr_step_history = zeros(NUM_STEPS, 2);  
true_path_history = zeros(NUM_STEPS, 2);      
estimated_path_history = zeros(NUM_STEPS, 2); 
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);
estimated_path_history(1, :) = true_path_history(1, :);
current_M = M_init;                     
M_history = zeros(NUM_STEPS, 1);        
M_history(1) = current_M;

%% 4. RUN SIMULATION
%==========================================================================
fprintf('运行 %d 模拟步 (GPU-PF-DTW)...\n', NUM_STEPS);
h_waitbar = waitbar(0, 'Running GPU-PF-DTW 2D Particle Filter...');
for t = 2:NUM_STEPS
   
    % --- 4a. 模拟真实运动 (仍在CPU) ---
    [true_state, pdr_step] = Get_Next_Step_2D(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step;
    full_true_path_history(t, :) = true_state;
    
    % --- 4b. 准备函数输入 (在CPU上准备, 然后发送到GPU) ---
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    
    pdr_history_for_function = zeros(SEQUENCE_LEN, 2);
    actual_len = end_idx - start_idx + 1;
    pdr_history_for_function(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
    
    live_sequence = zeros(1, SEQUENCE_LEN);
    path_segment = zeros(SEQUENCE_LEN, 3);
    path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
    
    % 在CPU上采样 (X, Y, Mag 都是 double)
    for k = 1:SEQUENCE_LEN
        pos_x = path_segment(k, 1);
        pos_y = path_segment(k, 2);
        if pos_x ~= 0 || pos_y ~= 0
            live_sequence(k) = interp2(X, Y, Mag, pos_x, pos_y, 'linear', 0);
        end
    end
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD;
    
    % *** [修复：小波降噪在CPU上执行] ***
    % 1. 在CPU上降噪
    live_sequence_denoised = Denoise_Visushrink_CPU(single(live_sequence), 'db4');
    
    % 2. 然后发送到GPU
    live_sequence_denoised_gpu = gpuArray(live_sequence_denoised);
    
    % --- 4c. 调用 2D 粒子滤波步骤 ---
    [particles_out_gpu, best_guess_gpu, dist] = Particle_Filter_DTW_Step_2D_GPU( ...
                                    particles, ...
                                    live_sequence_denoised_gpu, ...
                                    pdr_history_for_function, ... 
                                    geo_map, ... 
                                    process_noise, DTW_NOISE_STD);
    
    
    % --- 4d. APF 控制: 重注入 (在GPU上) ---
    M_current_step = size(particles_out_gpu, 1);
    N_reset = round(M_current_step * APFRATE);

    % === [已修复] ===
    % 在CPU上安全地运行 randperm, 然后将索引发送到GPU
    p_cpu = randperm(M_current_step); 
    indices_cpu = p_cpu(1:N_reset);  
    indices_to_reset = gpuArray(indices_cpu); 
    % ===============
    
    best_guess_state = gather(best_guess_gpu);
    
    % *** [GPU 类型修复] ***
    % 确保添加的噪声也是 'single'
    particles_out_gpu(indices_to_reset, 1) = best_guess_state(1) + randn(N_reset, 1, 'single', 'gpuArray') * pos_std;
    particles_out_gpu(indices_to_reset, 2) = best_guess_state(2) + randn(N_reset, 1, 'single', 'gpuArray') * pos_std;
    particles_out_gpu(indices_to_reset, 3) = best_guess_state(3) + randn(N_reset, 1, 'single', 'gpuArray') * ang_std;
    
    % --- APF 控制: M-Adaptation (在CPU上) ---
    particles_for_kld = gather(particles_out_gpu);
    dist_cpu = gather(dist); 
    
    M_new = Adapt_Particle_Count_DTW_KLD(particles_for_kld, dist_cpu, KLD.bin_size_xy, ...
                                MAP_X_LEN, MAP_Y_LEN, ...
                                KLD.epsilon, KLD.delta, ...
                                KLD.n_min, KLD.n_max, ...
                                KLD.DTW_THRESH_HIGH);
    
    if M_new ~= current_M
        particles_next_gpu = Adjust_Particle_Set_GPU(particles_out_gpu, M_new);
    else
        particles_next_gpu = particles_out_gpu; 
    end
    
    particles = particles_next_gpu;
    current_M = size(particles, 1);
                                    
    % --- 4e. 存储结果 (从GPU gather) ---
    true_path_history(t, :) = true_state(1:2);
    estimated_path_history(t, :) = best_guess_state(1:2); 
    M_history(t) = current_M; 
    
    % --- 4f. Update Waitbar ---
    waitbar(t/NUM_STEPS, h_waitbar, sprintf('GPU-PF-DTW (M=%d, dist=%.1f)', current_M, dist_cpu));
end
close(h_waitbar);
fprintf('模拟完成。\n');

%% 5. PLOTTING (在绘图前 gather 最终粒子)
%==========================================================================
fprintf('正在从GPU收集数据用于绘图...\n');
particles_for_plot = gather(particles);
map_for_plot = gather(geo_map.Mag_map);
X_for_plot = gather(geo_map.X_grid);
Y_for_plot = gather(geo_map.Y_grid);

fprintf('正在绘图...\n');
figure('Position', [100, 100, 1600, 500]);
% --- Plot 1: 2D Path ---
subplot(1, 3, 1); 
imagesc(X_for_plot(1,:), Y_for_plot(:,1), map_for_plot);
hold on;
axis xy; colormap('jet'); colorbar;
title('2D Path Tracking (GPU-PF-DTW)'); 
xlabel('X Position'); ylabel('Y Position');
plot(particles_for_plot(:, 1), particles_for_plot(:, 2), 'k.', 'MarkerSize', 2, 'DisplayName', 'Final Particle Cloud');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', 'True Path');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'DTW-APF Estimate'); 
plot(true_path_history(1, 1), true_path_history(1, 2), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'True Start');
plot(estimated_path_history(1, 1), estimated_path_history(1, 2), 'gs', 'MarkerSize', 12, 'DisplayName', 'Est. Start');
plot(true_path_history(end, 1), true_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'True End');
plot(estimated_path_history(end, 1), estimated_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', 'Est. End');
legend('show', 'Location', 'best');
axis equal; axis([0 MAP_X_LEN 0 MAP_Y_LEN]);

% --- Plot 2: Error (已在CPU) ---
subplot(1, 3, 2); 
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('Localization Error');
xlabel('Time Step'); ylabel('Error'); grid on; ylim([0, Inf]);

% --- Plot 3: Adaptive Particle Count (M_t) (已在CPU) ---
subplot(1, 3, 3);
plot(1:NUM_STEPS, M_history, 'b-', 'LineWidth', 2);
title('Adaptive Particle Count (M_t)');
xlabel('Time Step'); ylabel('Particle Count (M)'); grid on;
line([1, NUM_STEPS], [KLD.n_min, KLD.n_min], 'Color', 'red', 'LineStyle', '--');
line([1, NUM_STEPS], [KLD.n_max, KLD.n_max], 'Color', 'red', 'LineStyle', '--');
ylim([KLD.n_min*0.8, KLD.n_max*1.2]);


%% ========================================================================
%  == 子函数定义区 ==
%  ========================================================================

function [particles_out_gpu, best_guess_gpu, d_min_gpu] = Particle_Filter_DTW_Step_2D_GPU(particles_in_gpu, live_sequence_gpu, ...
                                            pdr_history_cpu, geo_map_gpu, process_noise, DTW_NOISE_STD)
%
% Function: Particle_Filter_DTW_Step_2D_GPU
% Summary: [GPU/PARFOR VERSION] (已修复 'single' 类型)
%
    
    M = size(particles_in_gpu, 1);
    L = size(pdr_history_cpu, 1);
    
    % --- 1. Propagation (Prediction) [ON GPU] ---
    % 确保所有输入和噪声都是 'single'
    last_pdr_step_cpu = single(pdr_history_cpu(end, :)); 
    
    % *** [GPU 类型修复] ***
    step_noise_gpu = randn(M, 1, 'single', 'gpuArray') * single(process_noise.step_std);
    theta_noise_gpu = randn(M, 1, 'single', 'gpuArray') * single(process_noise.theta_std);
    
    steps_gpu = last_pdr_step_cpu(1) + step_noise_gpu;
    d_thetas_gpu = last_pdr_step_cpu(2) + theta_noise_gpu;
    
    old_theta_gpu = particles_in_gpu(:, 3);
    
    new_theta_gpu = old_theta_gpu + d_thetas_gpu;
    new_x_gpu = particles_in_gpu(:, 1) + sin(new_theta_gpu) .* steps_gpu;
    new_y_gpu = particles_in_gpu(:, 2) + cos(new_theta_gpu) .* steps_gpu;
    
    % Boundary check on GPU
    map_max_x = geo_map_gpu.X_grid(1,end);
    map_max_y = geo_map_gpu.Y_grid(end,1);
    new_x_gpu = max(1, min(map_max_x, new_x_gpu));
    new_y_gpu = max(1, min(map_max_y, new_y_gpu));
    
    propagated_particles_gpu = [new_x_gpu, new_y_gpu, new_theta_gpu];
    
    % --- 2. Weighting ---
    
    % --- 2a. Path Reconstruction (Vectorized on GPU) ---
    % *** [GPU 类型修复] ***
    all_particle_paths_x_gpu = zeros(M, L, 'single', 'gpuArray');
    all_particle_paths_y_gpu = zeros(M, L, 'single', 'gpuArray');
    all_particle_thetas_gpu = zeros(M, L, 'single', 'gpuArray');
    
    all_particle_paths_x_gpu(:, L) = propagated_particles_gpu(:, 1);
    all_particle_paths_y_gpu(:, L) = propagated_particles_gpu(:, 2);
    all_particle_thetas_gpu(:, L) = propagated_particles_gpu(:, 3);
    
    for k = L-1:-1:1
        pdr_step_len = single(pdr_history_cpu(k+1, 1));
        pdr_d_theta = single(pdr_history_cpu(k+1, 2));
        
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
    
    % --- 2b. Generate Map Sequences (Vectorized on GPU) ---
    all_map_sequences_gpu = interp2(geo_map_gpu.X_grid, geo_map_gpu.Y_grid, geo_map_gpu.Mag_map, ...
                                  all_particle_paths_x_gpu, all_particle_paths_y_gpu, 'linear', 0);
                              
    % --- 2c. Calculate DTW & Weights (Parallel on CPU) ---
    all_map_sequences_cpu = gather(all_map_sequences_gpu);
    live_sequence_cpu = gather(live_sequence_gpu);
    live_sequence_cpu = live_sequence_cpu(:)'; 
    
    % *** [GPU 类型修复] ***
    % 确保CPU数组也是 'single'，以减少转换开销
    weights_cpu = zeros(M, 1, 'single');
    all_distances_cpu = zeros(M, 1, 'single');
    dtw_variance = single(DTW_NOISE_STD^2);
    
    parfor m = 1:M
        map_seq = all_map_sequences_cpu(m, :);
        
        % DTW 在 'single' 上工作得很好
        distance = dtw(live_sequence_cpu, map_seq);
        
        all_distances_cpu(m) = distance;
        weights_cpu(m) = 1./(1+(distance./dtw_variance).^2);
    end
    
    % --- 2d. Normalize Weights (ON GPU) ---
    weights_gpu = gpuArray(weights_cpu);
    all_distances_gpu = gpuArray(all_distances_cpu);
    
    sum_weights = sum(weights_gpu);
    if sum_weights > 1e-15
        weights_gpu = weights_gpu / sum_weights;
    else
        % *** [GPU 类型修复] ***
        weights_gpu = ones(M, 1, 'single', 'gpuArray') / M;
    end
    
    % --- 3. Resampling (ON GPU) ---
    cdf_gpu = cumsum(weights_gpu);
    % *** [GPU 类型修复] ***
    r_0_gpu = rand(1, 'single', 'gpuArray') / M;
    U_gpu = r_0_gpu + (gpuArray(single(transpose(1:M))) - 1) / M;
    
    [~, indices_gpu] = histc(U_gpu, [0; cdf_gpu]);
    indices_gpu = indices_gpu + 1; 
    indices_gpu = min(indices_gpu, M); 
    
    particles_out_gpu = propagated_particles_gpu(indices_gpu, :);
    
    % --- 4. Estimation (ON GPU) ---
    mean_pos_gpu = mean(particles_out_gpu(:, 1:2), 1);
    mean_sin_gpu = mean(sin(particles_out_gpu(:, 3)));
    mean_cos_gpu = mean(cos(particles_out_gpu(:, 3)));
    mean_theta_gpu = atan2(mean_sin_gpu, mean_cos_gpu);
    
    best_guess_gpu = [mean_pos_gpu, mean_theta_gpu];
    d_min_gpu = min(all_distances_gpu);
end

% -------------------------------------------------------------------------

function particles_next_gpu = Adjust_Particle_Set_GPU(particles_curr_gpu, M_new_cpu)
%
% Function: Adjust_Particle_Set_GPU
% Summary: [GPU VERSION] (已修复 'single' 类型)
%
    
    M_new = round(M_new_cpu);
    M_curr = size(particles_curr_gpu, 1);
    
    if M_new == M_curr
        particles_next_gpu = particles_curr_gpu;
        return;
        
    elseif M_new > M_curr
        % --- 3. Increase particles (on GPU) ---
        N_add = M_new - M_curr;
        indices_gpu = randi(M_curr, N_add, 1, 'gpuArray'); 
        
        particles_to_add_gpu = particles_curr_gpu(indices_gpu, :);
        
        % (可选: 添加 'single' 噪声)
        % noise_pos = (rand(N_add, 2, 'single', 'gpuArray') - 0.5) * 0.1; 
        % noise_ang = (rand(N_add, 1, 'single', 'gpuArray') - 0.5) * deg2rad(1); 
        % ...
        
        particles_next_gpu = [particles_curr_gpu; particles_to_add_gpu];
        
    else % M_new < M_curr
        % --- 4. 减少粒子 (on GPU) ---
        
        % === [已修复] ===
        % 在CPU上安全地运行 randperm, 然后将索引发送到GPU
        p_cpu = randperm(M_curr); 
        indices_cpu = p_cpu(1:M_new); 
        indices_to_keep_gpu = gpuArray(indices_cpu); 
        % ===============
        
        particles_next_gpu = particles_curr_gpu(indices_to_keep_gpu, :);
    end
end

% -------------------------------------------------------------------------

function denoised_data = Denoise_Visushrink_CPU(data, wavelet)
%DENOISE_VISUSHRINK (CPU Version) 
%   修复版本：在CPU上执行小波降噪
%   输入 'data' 必须是 'single' 类型
%
    n = length(data);
    
    % --- 步骤 1: 小波分解 ---
    wlevel = wmaxlev(n, wavelet);
    [c, l] = wavedec(data, wlevel, wavelet); 

    % --- 步骤 2: 计算"通用门槛" ---
    
    % 2a. 估计噪音水平 sigma
    cD1_len = l(end-1);
    cD1_start_index = sum(l(1:end-2)) + 1;
    cD1_end_index = cD1_start_index + cD1_len - 1;
    
    cD1 = c(cD1_start_index : cD1_end_index);
    
    % 2b. 计算门槛
    sigma_est = mad(cD1, 1) / single(0.6745); 
    threshold = sigma_est * single(sqrt(2 * log(single(n))));

    % --- 步骤 3: 应用软阈值 ---
    cA = c(1:l(1));
    detail_coeffs = c(l(1)+1 : end);
    
    denoised_details = wthresh(detail_coeffs, 's', threshold);
    
    denoised_c = [cA, denoised_details];

    % --- 步骤 4: 小波重构 ---
    denoised_data_temp = waverec(denoised_c, l, wavelet);
    denoised_data = denoised_data_temp(1:n);
end