%==========================================================================
% SCRIPT: MAIN_2D_PFDTW_KLD_GPU (A/B 对比版本)
%
% 描述: 
%   此脚本已修改为在GPU上管理粒子和地图数据。
%   它将并行运行两个粒子滤波函数 (PF-DTW-GPU 和 PF-DDTW-CPFTest)
%   以便在相同的地图和路径下进行直接比较。
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
DTW_NOISE_STD = 30; % 这个值用于两个滤波器
process_noise.step_std = 0.05;          
process_noise.theta_std = deg2rad(1); 
pos_std = 5.0;           
ang_std = 0.2;            
KLD.n_min = 500;          
KLD.n_max = 100000;       
KLD.bin_size_xy = 0.5;    
KLD.epsilon = 0.05;       
KLD.delta = 0.01;         
KLD.DTW_THRESH_HIGH = 83.9; % KLD 阈值 (DDTW的阈值可能需要调整)
sigma = 0.05;
eps = 0.5;
APFRATE = 0.08;
% 2. MAP GENERATION (发送到GPU, 使用 'single')
%==========================================================================
% fprintf('生成地磁地图并发送到GPU...\n');
% [X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
% Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
% Mag = imgaussfilt(Mag_raw, 3.0); 
% % *** [GPU 类型修复] ***
% geo_map.X_grid = gpuArray(single(X));
% geo_map.Y_grid = gpuArray(single(Y));
% geo_map.Mag_map = gpuArray(single(Mag));
% fprintf('地图已在GPU上 (single 精度)。\n');
% %% 2. MAP GENERATION (V-Final: "拉伸"真实地图以适应模拟器)
% %==========================================================================
fprintf('正在从 my_uji_map.mat 加载真实地图...\n');
MAP_FILENAME = 'my_hybrid_map.mat';

% 1. 检查地图文件是否存在
if ~exist(MAP_FILENAME, 'file')
    error('未找到地图文件: %s. \n请先运行 Build_UJI_Map.m 脚本来创建它。', MAP_FILENAME);
end

% 2. [关键] 加载 .mat 文件
%    这会在工作区中创建 'geo_map_cpu' 变量
load(MAP_FILENAME); 

fprintf('...已加载真实地图 (原始尺寸: %d x %d)。\n', ...
        size(geo_map_cpu.Mag_map, 1), size(geo_map_cpu.Mag_map, 2));

% 3. [!! 核心修复 !!] "拉伸"真实地图以匹配模拟器
%    我们使用 imresize 将加载的真实地图 "拉伸" 到
%    你在 %% 1. 中定义的 50x50 模拟尺寸。

%    从 %% 1. 获取模拟网格尺寸
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);

%    拉伸地图
Mag = imresize(geo_map_cpu.Mag_map, [MAP_Y_LEN, MAP_X_LEN]);
Mag = imgaussfilt(Mag, 1.0); % 轻微平滑一下拉伸后的图像

fprintf('...真实地图已被 "拉伸" 到 %d x %d 以匹配模拟器。\n', ...
        MAP_Y_LEN, MAP_X_LEN);

% 4. [兼容性] 创建你的 %% 4. (模拟器) 需要的变量
%    现在 X, Y, Mag 都有了, 第 161 行的 interp2 不会报错
X_for_interp = X; % (保留 X, Y, Mag 以便你将来调试)
Y_for_interp = Y;
Mag_for_interp = Mag;

% 5. [兼容性] 创建你的 *子函数* 需要的 'geo_map' 结构体
%    (你的子函数 Particle_Filter_... 也需要这个)
fprintf('发送 "拉伸后" 的地图到 GPU...\n');
geo_map.X_grid = gpuArray(single(X));
geo_map.Y_grid = gpuArray(single(Y));
geo_map.Mag_map = gpuArray(single(Mag));
fprintf('地图已在GPU上 (single 精度)。\n');

%% 3. INITIALIZATION (为两个滤波器创建 'single' 粒子)
%==========================================================================
fprintf('在GPU上为两个滤波器初始化粒子...\n');
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; 
INIT_POS_STD = 5.0; 
INIT_ANG_STD = 0.5; 

% --- 滤波器 1 (GPU-DTW) 的粒子和历史 ---
particles_gpu = zeros(M_init, 3, 'single', 'gpuArray'); 
particles_gpu(:, 1) = single(true_state(1)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_POS_STD;
particles_gpu(:, 2) = single(true_state(2)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_POS_STD;
particles_gpu(:, 3) = single(true_state(3)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_ANG_STD;
current_M_gpu = M_init;
M_history_gpu = zeros(NUM_STEPS, 1);
M_history_gpu(1) = current_M_gpu;
estimated_path_history_gpu = zeros(NUM_STEPS, 2); 
estimated_path_history_gpu(1, :) = true_state(1:2);

% --- 滤波器 2 (CPF-DDTW) 的粒子和历史 ---
% 使用完全相同的初始状态
particles_cpf = particles_gpu; 
current_M_cpf = M_init;
M_history_cpf = zeros(NUM_STEPS, 1);
M_history_cpf(1) = current_M_cpf;
estimated_path_history_cpf = zeros(NUM_STEPS, 2); 
estimated_path_history_cpf(1, :) = true_state(1:2);


% --- 共享的真实路径历史 (CPU) ---
full_true_path_history = zeros(NUM_STEPS, 3); 
full_pdr_step_history = zeros(NUM_STEPS, 2);  
true_path_history = zeros(NUM_STEPS, 2);      
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);

%% 4. RUN SIMULATION (并行对比)
%==========================================================================
fprintf('运行 %d 模拟步 (A/B 对比)...\n', NUM_STEPS);
h_waitbar = waitbar(0, '运行 A/B 对比...');
for t = 2:NUM_STEPS
   
    % --- 4a. 模拟真实运动 (共享) ---
    [true_state, pdr_step] = Get_Next_Step_2D(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step;
    full_true_path_history(t, :) = true_state;
    true_path_history(t, :) = true_state(1:2); % 存储真实路径
    
    % --- 4b. 准备函数输入 (共享) ---
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    
    pdr_history_for_function = zeros(SEQUENCE_LEN, 2);
    actual_len = end_idx - start_idx + 1;
    pdr_history_for_function(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
    
    live_sequence = zeros(1, SEQUENCE_LEN);
    path_segment = zeros(SEQUENCE_LEN, 3);
    path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
    
    for k = 1:SEQUENCE_LEN
        pos_x = path_segment(k, 1);
        pos_y = path_segment(k, 2);
        if pos_x ~= 0 || pos_y ~= 0
            live_sequence(k) = interp2(X, Y, Mag, pos_x, pos_y, 'linear', 0);
        end
    end
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD;
    live_sequence_denoised = Denoise_Visushrink_CPU(single(live_sequence), 'db4');
    live_sequence_denoised_gpu = gpuArray(live_sequence_denoised);
    
    % --- 4c. 【滤波器 1: GPU-DTW】 ---
    [particles_out_gpu, best_guess_gpu, dist_gpu] = Particle_Filter_DTW_Step_2D_GPU( ...
                                    particles_gpu, ...
                                    live_sequence_denoised_gpu, ...
                                    pdr_history_for_function, ... 
                                    geo_map, ... 
                                    process_noise, DTW_NOISE_STD);
    
    % --- 4d. 【滤波器 2: CPF-DDTW】 ---
    % *** 解开注释 ***
     [particles_out_cpf, best_guess_cpf, dist_cpf] = Particle_Filter_DTW_Step_2D_CPFTest( ...
                                    particles_cpf, ...
                                    live_sequence_denoised_gpu, ...
                                    pdr_history_for_function, ... 
                                    geo_map, ... 
                                    process_noise, DTW_NOISE_STD);
    
    % --- 4e. 【滤波器 1: APF 和 KLD】 ---
    N_reset_gpu = round(current_M_gpu * APFRATE);
    p_cpu_gpu = randperm(current_M_gpu); 
    indices_to_reset_gpu = gpuArray(p_cpu_gpu(1:N_reset_gpu)); 
    best_guess_state_gpu = gather(best_guess_gpu);
    
    particles_out_gpu(indices_to_reset_gpu, 1) = best_guess_state_gpu(1) + randn(N_reset_gpu, 1, 'single', 'gpuArray') * pos_std;
    particles_out_gpu(indices_to_reset_gpu, 2) = best_guess_state_gpu(2) + randn(N_reset_gpu, 1, 'single', 'gpuArray') * pos_std;
    particles_out_gpu(indices_to_reset_gpu, 3) = best_guess_state_gpu(3) + randn(N_reset_gpu, 1, 'single', 'gpuArray') * ang_std;
    
    dist_cpu_gpu = gather(dist_gpu); 
    M_new_gpu = Adapt_Particle_Count_DTW_KLD(gather(particles_out_gpu), dist_cpu_gpu, KLD.bin_size_xy, ...
                                MAP_X_LEN, MAP_Y_LEN, KLD.epsilon, KLD.delta, ...
                                KLD.n_min, KLD.n_max, KLD.DTW_THRESH_HIGH);
    
    if M_new_gpu ~= current_M_gpu
        particles_gpu = Adjust_Particle_Set_GPU(particles_out_gpu, M_new_gpu);
    else
        particles_gpu = particles_out_gpu; 
    end
    current_M_gpu = size(particles_gpu, 1);
    
    % --- 4f. 【滤波器 2: APF 和 KLD】 ---
    N_reset_cpf = round(current_M_cpf * APFRATE);
    p_cpu_cpf = randperm(current_M_cpf); 
    indices_to_reset_cpf = gpuArray(p_cpu_cpf(1:N_reset_cpf)); 
    best_guess_state_cpf = gather(best_guess_cpf);
    
    particles_out_cpf(indices_to_reset_cpf, 1) = best_guess_state_cpf(1) + randn(N_reset_cpf, 1, 'single', 'gpuArray') * pos_std;
    particles_out_cpf(indices_to_reset_cpf, 2) = best_guess_state_cpf(2) + randn(N_reset_cpf, 1, 'single', 'gpuArray') * pos_std;
    particles_out_cpf(indices_to_reset_cpf, 3) = best_guess_state_cpf(3) + randn(N_reset_cpf, 1, 'single', 'gpuArray') * ang_std;
    
    dist_cpu_cpf = gather(dist_cpf); 
    % 注意: DDTW的距离阈值 (KLD.DTW_THRESH_HIGH) 可能需要针对DDTW进行调整
    M_new_cpf = Adapt_Particle_Count_DTW_KLD(gather(particles_out_cpf), dist_cpu_cpf, KLD.bin_size_xy, ...
                                MAP_X_LEN, MAP_Y_LEN, KLD.epsilon, KLD.delta, ...
                                KLD.n_min, KLD.n_max, KLD.DTW_THRESH_HIGH);
    
    if M_new_cpf ~= current_M_cpf
        particles_cpf = Adjust_Particle_Set_GPU(particles_out_cpf, M_new_cpf);
    else
        particles_cpf = particles_out_cpf; 
    end
    current_M_cpf = size(particles_cpf, 1);

    % --- 4g. 存储结果 (从GPU gather) ---
    estimated_path_history_gpu(t, :) = best_guess_state_gpu(1:2); 
    M_history_gpu(t) = current_M_gpu; 
    
    estimated_path_history_cpf(t, :) = best_guess_state_cpf(1:2); 
    M_history_cpf(t) = current_M_cpf; 
    
    % --- 4h. Update Waitbar ---
    waitbar(t/NUM_STEPS, h_waitbar, sprintf('对比 (M_{gpu}=%d, M_{cpf}=%d)', current_M_gpu, current_M_cpf));
end
close(h_waitbar);
fprintf('模拟完成。\n');
%% 5. PLOTTING (对比图)
%==========================================================================
fprintf('正在从GPU收集数据用于绘图...\n');
% 只需要收集一个滤波器的最终粒子云 (或不收集)
% particles_for_plot = gather(particles_cpf); % 以CPF的为例
map_for_plot = gather(geo_map.Mag_map);
X_for_plot = gather(geo_map.X_grid);
Y_for_plot = gather(geo_map.Y_grid);
fprintf('正在绘图...\n');

figure('Position', [100, 100, 1600, 600]); % 调整图窗大小

% --- Plot 1: 2D Path (对比) ---
subplot(1, 3, 1); 
imagesc(X_for_plot(1,:), Y_for_plot(:,1), map_for_plot);
hold on;
axis xy; colormap('jet'); colorbar;
title('2D 路径追踪对比'); 
xlabel('X Position'); ylabel('Y Position');
% plot(particles_for_plot(:, 1), particles_for_plot(:, 2), 'k.', 'MarkerSize', 1, 'DisplayName', 'Final Cloud (CPF)');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', 'True Path (真实路径)');
plot(estimated_path_history_gpu(:, 1), estimated_path_history_gpu(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'Est: GPU (标准DTW)'); 
plot(estimated_path_history_cpf(:, 1), estimated_path_history_cpf(:, 2), 'g-.*', 'LineWidth', 1.5, 'DisplayName', 'Est: CPF (DDTW)'); 
legend('show', 'Location', 'best');
axis equal; axis([0 MAP_X_LEN 0 MAP_Y_LEN]);

% --- Plot 2: Error (对比) ---
subplot(1, 3, 2); 
errors_gpu = sqrt(sum((true_path_history - estimated_path_history_gpu).^2, 2));
errors_cpf = sqrt(sum((true_path_history - estimated_path_history_cpf).^2, 2));
plot(1:NUM_STEPS, errors_gpu, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Error: GPU (标准DTW)');
hold on;
plot(1:NUM_STEPS, errors_cpf, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Error: CPF (DDTW)');
title('定位误差对比');
xlabel('Time Step'); ylabel('Error'); grid on; ylim([0, Inf]);
legend('show', 'Location', 'best');

% --- Plot 3: Adaptive Particle Count (M_t) (对比) ---
subplot(1, 3, 3);
plot(1:NUM_STEPS, M_history_gpu, 'r--', 'LineWidth', 2, 'DisplayName', 'M_t (GPU)');
hold on;
plot(1:NUM_STEPS, M_history_cpf, 'g-', 'LineWidth', 2, 'DisplayName', 'M_t (CPF)');
title('自适应粒子数量 (M_t)');
xlabel('Time Step'); ylabel('Particle Count (M)'); grid on;
line([1, NUM_STEPS], [KLD.n_min, KLD.n_min], 'Color', 'black', 'LineStyle', '--');
line([1, NUM_STEPS], [KLD.n_max, KLD.n_max], 'Color', 'black', 'LineStyle', '--');
legend('show', 'Location', 'best');
% ylim([KLD.n_min*0.8, KLD.n_max*1.2]); % 调整Y轴限制

%% 6. 输出统计信息
%==========================================================================
% 计算从第100步到结束的统计数据 (忽略初始收敛阶段)
start_step_for_stats = 100;
errors_gpu_stable = errors_gpu(start_step_for_stats:end);
errors_cpf_stable = errors_cpf(start_step_for_stats:end);

fprintf('\n--- 滤波器性能对比 (从第 %d 步开始) ---\n', start_step_for_stats);

fprintf('方法 1: GPU (标准 DTW)\n');
fprintf('  - 平均误差 (Mean):   %.4f\n', mean(errors_gpu_stable));
fprintf('  - 中位数误差 (Median): %.4f\n', median(errors_gpu_stable));
fprintf('  - 误差标准差 (Std):  %.4f\n', std(errors_gpu_stable));
fprintf('  - 误差方差 (Var):    %.4f\n', var(errors_gpu_stable));
fprintf('  - 均方根误差 (RMSE): %.4f\n', sqrt(mean(errors_gpu_stable.^2)));
fprintf('  - 最大误差 (Max):    %.4f\n', max(errors_gpu_stable));

fprintf('\n方法 2: CPF (DDTW)\n');
fprintf('  - 平均误差 (Mean):   %.4f\n', mean(errors_cpf_stable));
fprintf('  - 中位数误差 (Median): %.4f\n', median(errors_cpf_stable));
fprintf('  - 误差标准差 (Std):  %.4f\n', std(errors_cpf_stable));
fprintf('  - 误差方差 (Var):    %.4f\n', var(errors_cpf_stable));
fprintf('  - 均方根误差 (RMSE): %.4f\n', sqrt(mean(errors_cpf_stable.^2)));
fprintf('  - 最大误差 (Max):    %.4f\n', max(errors_cpf_stable));

fprintf('-----------------------------------------------\n');


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