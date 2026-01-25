%==========================================================================
% SCRIPT: MAIN_2D_PFDTW_KLD_GPU (V5.1 - 最终索引修复版)
%
% 描述: 
%   此脚本为最终的、已修复所有问题的版本。
%   - 循环遍历所有 "tt" 文件
%   - 使用 V5 逻辑加载 "完美" PDR 步骤和 "真实" 磁场测量值
%   - 修复了 "循环真值" 问题 (移除了 interp2)
%   - 修复了 for 循环的 "end" 位置
%   - [V5.1] 修复了 t=2 时 data_start_idx=0 的索引错误
%
% 依赖:
%   - Parallel Computing Toolbox
%   - Wavelet Toolbox
%   - Get_Next_Step_2D_UJI.m (V5 版本 - 您必须有这个文件)
%   - Adapt_Particle_Count_DTW_KLD.m (您必须有这个文件)
%   - my_hybrid_map.mat (您必须有这个文件)
%==========================================================================
clc;
clear;
close all;

% --- 0. GPU 和并行初始化 ---
fprintf('正在检查GPU...\n');
try
    gpu = gpuDevice(1); 
    fprintf('在 %s (GPU) 上运行。\n', gpu.Name);
catch ME
    warning('未找到GPU。将回退到CPU运行。');
    gpuArray = @(x) x; 
    gather = @(x) x;
    randn = @(varargin) randn(varargin{:}, 'single');
    zeros = @(varargin) zeros(varargin{:}, 'single');
end

%% 1. SIMULATION PARAMETERS
%==========================================================================
M_init = 500;            
NUM_STEPS_Fallback = 500; % [重命名] 仅在 TESTING=0 时使用
SEQUENCE_LEN = 100;       
MAP_X_LEN = 50;          
MAP_Y_LEN = 50;           
SENSOR_NOISE_STD = 0.5;   % 这是添加到 *真实* 噪声上的额外噪声
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
KLD.DTW_THRESH_HIGH = 1010.9; 
sigma = 0.05;
eps = 0.5;
APFRATE = 0.1;
TESTING = 1; % 1 is using the UJI while 0 is the total fake one.

% --- [V5 修复: 使用 string (双引号) 定义路径] ---
root = ".\tests";
filenames = ["tt01.txt", "tt02.txt", "tt03.txt","tt04.txt", "tt05.txt", "tt06.txt", "tt07.txt", "tt08.txt", "tt09.txt", "tt10.txt", "tt11.txt"]; 

%% 2. MAP GENERATION (V-Final: "拉伸"真实地图以适应模拟器)
%==========================================================================
fprintf('正在从 my_uji_map.mat 加载真实地图...\n');
MAP_FILENAME = 'my_hybrid_map.mat';
if ~exist(MAP_FILENAME, 'file')
    error('未找到地图文件: %s. \n请先运行 Build_UJI_Map.m 脚本来创建它。', MAP_FILENAME);
end
load(MAP_FILENAME); % 加载 'geo_map_cpu'
fprintf('...已加载真实地图 (原始尺寸: %d x %d)。\n', ...
        size(geo_map_cpu.Mag_map, 1), size(geo_map_cpu.Mag_map, 2));

% [!! 警告 !!]：如前所述, imresize 会扭曲地图。
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag = imresize(geo_map_cpu.Mag_map, [MAP_Y_LEN, MAP_X_LEN]);
Mag = imgaussfilt(Mag, 1.0); 
fprintf('...真实地图已被 "拉伸" 到 %d x %d 以匹配模拟器。\n', ...
        MAP_Y_LEN, MAP_X_LEN);

% 发送 *拉伸后* 的地图到 GPU
fprintf('发送 "拉伸后" 的地图到 GPU...\n');
geo_map.X_grid = gpuArray(single(X));
geo_map.Y_grid = gpuArray(single(Y));
geo_map.Mag_map = gpuArray(single(Mag));
fprintf('地图已在GPU上 (single 精度)。\n');

%==========================================================================
% [!! 关键结构修复: 循环从这里开始 !!]
%==========================================================================
for i = 1:length(filenames)
    
    % [V5 修复: 正确拼接 string]
    filename = root + "\" + filenames(i);
    
    fprintf('\n======================================================\n');
    fprintf('开始测试: %s\n', filename);
    
%% 3. INITIALIZATION (在GPU上创建 'single' 粒子)
%==========================================================================
fprintf('在GPU上初始化粒子...\n');
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; 
INIT_POS_STD = 5.0; 
INIT_ANG_STD = 0.5; 
% 直接在GPU上创建 'single' 类型的粒子
particles = zeros(M_init, 3, 'single', 'gpuArray'); % 在GPU上创建
particles(:, 1) = single(true_state(1)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_POS_STD;
particles(:, 2) = single(true_state(2)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_POS_STD;
particles(:, 3) = single(true_state(3)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_ANG_STD;

% --- History Logs (日志仍在CPU上) ---
if TESTING == 1
    % [V5 修复: 加载所有 UJI 数据]
    UJI_Data_CPU = Get_Next_Step_2D_UJI('init', filename);
    
    if UJI_Data_CPU.num_steps == 0
        warning('在 %s 中未检测到步伐。跳过此文件。', filename);
        continue; % 跳到下一个 for 循环 (下一个文件)
    end
    
    NUM_STEPS = UJI_Data_CPU.num_steps + 1; % (例如 73 + 1 = 74)
    
    % 存储 *真实* 的磁场读数历史 (用于 'live_sequence')
    full_real_mag_history_CPU = UJI_Data_CPU.mag_readings;
    
    % 存储 *完美* 的 PDR 步骤历史 (用于粒子传播 和 true_path 模拟)
    full_perfect_pdr_history_CPU = UJI_Data_CPU.pdr_steps;
    
else
    NUM_STEPS = NUM_STEPS_Fallback;
    full_real_mag_history_CPU = randn(NUM_STEPS, 3) * 50; 
    full_perfect_pdr_history_CPU = [ones(NUM_STEPS, 1)*0.7, zeros(NUM_STEPS, 1)];
end

% [V5 修复: 使用正确的 NUM_STEPS 初始化]
full_true_path_history = zeros(NUM_STEPS, 3); 
full_pdr_step_history = zeros(NUM_STEPS, 2);  % [V5.1] 此变量现在也由 "完美" PDR 填充
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
fprintf('运行 %d 模拟步 (GPU-PF-DTW)...\n', NUM_STEPS - 1);
h_waitbar = waitbar(0, 'Running GPU-PF-DTW 2D Particle Filter...');

for t = 2:NUM_STEPS
   
    % --- 4a. 模拟真实运动 (仍在CPU) ---
    % [!] 这部分现在 *只* 用于生成 "Ground Truth" 路径以供绘图
    % 它使用 V5 函数的 "执行模式"
    [true_state, pdr_step_simulated] = Get_Next_Step_2D_UJI(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step_simulated; % [V5.1] 改为存储模拟的PDR
    full_true_path_history(t, :) = true_state;
    
    
    % --- 4b. 准备函数输入 (V5.1 索引修复) ---
    
    % [V5.1 修复] 我们需要索引 (t-1) 的数据点
    % data_end_idx 是我们 PDR/Mag 列表中的最新索引
    data_end_idx = t - 1; 
    % data_start_idx 是序列的起始索引, 最小为 1
    data_start_idx = max(1, data_end_idx - SEQUENCE_LEN + 1);
    
    % --- 准备 PDR (控制) 输入 ---
    pdr_history_for_function = zeros(SEQUENCE_LEN, 2, 'single');
    
    % 提取 "完美" PDR 步骤的真实片段
    pdr_segment = full_perfect_pdr_history_CPU(data_start_idx:data_end_idx, :);
    actual_len = size(pdr_segment, 1);
    
    % 将其放入 pdr_history_for_function 的末尾
    pdr_history_for_function(end-actual_len+1:end, :) = single(pdr_segment);
    
    
    % --- 准备 Mag (测量) 输入 [核心修复] ---
    live_sequence = zeros(1, SEQUENCE_LEN);
    
    % [V5.1 修复] 使用 *完全相同* 的索引
    real_mag_data_segment = full_real_mag_history_CPU(data_start_idx:data_end_idx, :);
    real_mag_norms = sqrt(sum(real_mag_data_segment.^2, 2));
    
    % [V5.1 修复] 使用 pdr_segment 中的 actual_len
    live_sequence(end-actual_len+1:end) = real_mag_norms'; 
    
    % (可选) 添加额外的传感器"电子"噪声
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD; 
    
    
    % --- 降噪 ---
    live_sequence_denoised = Denoise_Visushrink_CPU(single(live_sequence), 'db4');
    live_sequence_denoised_gpu = gpuArray(live_sequence_denoised);
    
    % --- 4c. 调用 2D 粒子滤波步骤 ---
     % [particles_out_gpu, best_guess_gpu, dist] = Particle_Filter_DTW_Step_2D_GPU( ...
     %                                particles, ...
     %                                live_sequence_denoised_gpu, ...
     %                                pdr_history_for_function, ... 
     %                                geo_map, ... 
     %                                process_noise, DTW_NOISE_STD);
     [particles_out_gpu, best_guess_gpu, dist] = Particle_Filter_DTW_Step_2D_CPFTest( ...
                                    particles, ...
                                    live_sequence_denoised_gpu, ...
                                    pdr_history_for_function, ... 
                                    geo_map, ... 
                                    process_noise, DTW_NOISE_STD);

                                
    % --- 4d. APF 控制: 重注入 (在GPU上) ---
    M_current_step = size(particles_out_gpu, 1);
    N_reset = round(M_current_step * APFRATE);
    p_cpu = randperm(M_current_step); 
    indices_cpu = p_cpu(1:N_reset);  
    indices_to_reset = gpuArray(indices_cpu); 
    
    best_guess_state = gather(best_guess_gpu);
    
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

% [V5] 为每个文件创建一个新图像
figure('Name', char(filename), 'Position', [100, 100, 1600, 500]);

% --- Plot 1: 2D Path ---
subplot(1, 3, 1); 
imagesc(X_for_plot(1,:), Y_for_plot(:,1), map_for_plot);
hold on;
axis xy; colormap('jet'); colorbar;
title(['2D Path Tracking: ' char(filenames(i))]); 
xlabel('X Position'); ylabel('Y Position');
plot(particles_for_plot(:, 1), particles_for_plot(:, 2), 'k.', 'MarkerSize', 2, 'DisplayName', 'Final Particle Cloud');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', 'True Path (Simulated)');
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
mean_errors = mean(errors);
disp(["文件 " + filenames(i) + " 的平均误差为: ", mean_errors])

% --- Plot 3: Adaptive Particle Count (M_t) (已在CPU) ---
subplot(1, 3, 3);
plot(1:NUM_STEPS, M_history, 'b-', 'LineWidth', 2);
title('Adaptive Particle Count (M_t)');
xlabel('Time Step'); ylabel('Particle Count (M)'); grid on;
line([1, NUM_STEPS], [KLD.n_min, KLD.n_min], 'Color', 'red', 'LineStyle', '--');
line([1, NUM_STEPS], [KLD.n_max, KLD.n_max], 'Color', 'red', 'LineStyle', '--');
ylim([KLD.n_min*0.8, KLD.n_max*1.2]);

% --- [V5 修复: 保存图像] ---
file_number_str = sprintf('%02d', i);
save_filename = "DDTWCPF_tt" + file_number_str + ".png";
fprintf('正在保存图像: %s\n', save_filename);
exportgraphics(gcf, save_filename, 'Resolution', 150);
    
%==========================================================================
% [!! 关键结构修复: FOR 循环在这里结束 !!]
%==========================================================================
end 


%% ========================================================================
%  == 子函数定义区 (必须在主脚本代码和循环 *之后*) ==
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
    last_pdr_step_cpu = single(pdr_history_cpu(end, :)); 
    
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
    all_particle_paths_x_gpu = zeros(M, L, 'single', 'gpuArray');
    all_particle_paths_y_gpu = zeros(M, L, 'single', 'gpuArray');
    all_particle_thetas_gpu = zeros(M, L, 'single', 'gpuArray');
    
    all_particle_paths_x_gpu(:, L) = propagated_particles_gpu(:, 1);
    all_particle_paths_y_gpu(:, L) = propagated_particles_gpu(:, 2);
    all_particle_thetas_gpu(:, L) = propagated_particles_gpu(:, 3);
    
    % [V5.1] pdr_history_cpu 是 double, 但子函数需要 'single'
    pdr_history_single = single(pdr_history_cpu);

    for k = L-1:-1:1
        pdr_step_len = pdr_history_single(k+1, 1);
        pdr_d_theta = pdr_history_single(k+1, 2);
        
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
    
    weights_cpu = zeros(M, 1, 'single');
    all_distances_cpu = zeros(M, 1, 'single');
    dtw_variance = single(DTW_NOISE_STD^2);
    
    parfor m = 1:M
        map_seq = all_map_sequences_cpu(m, :);
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
        weights_gpu = ones(M, 1, 'single', 'gpuArray') / M;
    end
    
    % --- 3. Resampling (ON GPU) ---
    cdf_gpu = cumsum(weights_gpu);
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
        particles_next_gpu = [particles_curr_gpu; particles_to_add_gpu];
        
    else % M_new < M_curr
        % --- 4. 减少粒子 (on GPU) ---
        p_cpu = randperm(M_curr); 
        indices_cpu = p_cpu(1:M_new); 
        indices_to_keep_gpu = gpuArray(indices_cpu); 
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
    if wlevel < 1 % 确保至少有 1 级
        denoised_data = data;
        return;
    end
    [c, l] = wavedec(data, wlevel, wavelet); 
    
    % --- 步骤 2: 计算"通用门槛" ---
    cD1_len = l(end-1);
    cD1_start_index = sum(l(1:end-2)) + 1;
    cD1_end_index = cD1_start_index + cD1_len - 1;
    
    % 确保 cD1 至少有 2 个元素以供 mad 使用
    if cD1_len >= 2 && cD1_end_index <= length(c)
        cD1 = c(cD1_start_index : cD1_end_index);
        sigma_est = mad(cD1, 1) / single(0.6745); 
        threshold = sigma_est * single(sqrt(2 * log(single(n))));
    else
        threshold = single(0); % 如果 cD1 太短, 不进行阈值处理
    end

    % --- 步骤 3: 应用软阈值 ---
    cA = c(1:l(1));
    detail_coeffs = c(l(1)+1 : end);
    
    denoised_details = wthresh(detail_coeffs, 's', threshold);
    
    denoised_c = [cA, denoised_details];
    
    % --- 步骤 4: 小波重构 ---
    denoised_data_temp = waverec(denoised_c, l, wavelet);
    % 确保输出长度与输入一致
    if length(denoised_data_temp) > n
        denoised_data = denoised_data_temp(1:n);
    else
        denoised_data = denoised_data_temp;
    end
end