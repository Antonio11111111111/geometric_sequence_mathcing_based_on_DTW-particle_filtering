%==========================================================================
% SCRIPT: MAIN_2D_PFDTW_KLD_GPU (V10.1 - 终极去零动态窗口修复版)
%==========================================================================
clc;
clear;
close all;

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
SEQUENCE_LEN = 20;       
MAP_X_LEN = 50;  
MAP_Y_LEN = 50;  
SENSOR_NOISE_STD = 0.5;   
DTW_NOISE_STD = 30;      

% --- 回归理性的黄金参数 ---
process_noise.step_std = 0.1;          
process_noise.theta_std = deg2rad(5);  
pos_std = 1.0;           
ang_std = 0.2;           
APFRATE = 0.05;          
INIT_POS_STD = 2.0; 
INIT_ANG_STD = 0.5; % 这里把漏掉的初始航向角标准差补回来了
% --------------------------

KLD.n_min = 500;          
KLD.n_max = 100000;       
KLD.bin_size_xy = 0.5;    
KLD.epsilon = 0.05;       
KLD.delta = 0.01;         
KLD.DTW_THRESH_HIGH = 10086; 
sigma = 0.05;
eps = 0.5;
TESTING = 1; 
NUM_STEPS_Fallback = 500; 

root = 'D:\dachuang\geometric_sequence_mathcing_based_on_DTW-particle_filtering\MAIN_Algorithm\tests'; 
filenames = ["tt01.txt", "tt02.txt", "tt03.txt","tt04.txt", "tt05.txt", "tt06.txt", "tt07.txt", "tt08.txt", "tt09.txt", "tt10.txt", "tt11.txt"]; 

%% 2. MAP GENERATION 
%==========================================================================
fprintf('正在从 my_uji_map.mat 加载真实地图...\n');
MAP_FILENAME = 'my_hybrid_map.mat';
if ~exist(MAP_FILENAME, 'file')
    error('未找到地图文件: %s. \n请先运行 Build_UJI_Map.m 脚本来创建它。', MAP_FILENAME);
end
load(MAP_FILENAME); 
fprintf('...已加载真实地图 (原始尺寸: %d x %d)。\n', ...
        size(geo_map_cpu.Mag_map, 1), size(geo_map_cpu.Mag_map, 2));

[MAP_Y_LEN_REAL, MAP_X_LEN_REAL] = size(geo_map_cpu.Mag_map); 
MAP_X_LEN = MAP_X_LEN_REAL; 
MAP_Y_LEN = MAP_Y_LEN_REAL; 
fprintf('...模拟器尺寸已重置为地图原始尺寸: %d x %d。\n', ...
        MAP_Y_LEN, MAP_X_LEN);

[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN); 
Mag = geo_map_cpu.Mag_map; 

fprintf('发送 "原始" 地图到 GPU...\n');
geo_map.X_grid = gpuArray(single(X));
geo_map.Y_grid = gpuArray(single(Y));
geo_map.Mag_map = gpuArray(single(Mag));
fprintf('地图已在GPU上 (single 精度)。\n');

%==========================================================================
for i = 1:length(filenames)
    
    filename = char(fullfile(root, filenames(i)));
    
    fprintf('\n======================================================\n');
    fprintf('开始测试: %s\n', filename);
    
%% 3. INITIALIZATION 
%==========================================================================
fprintf('在GPU上初始化粒子...\n');
true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)]; 

particles = zeros(M_init, 3, 'single', 'gpuArray'); 
particles(:, 1) = single(true_state(1)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_POS_STD;
particles(:, 2) = single(true_state(2)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_POS_STD;
particles(:, 3) = single(true_state(3)) + randn(M_init, 1, 'single', 'gpuArray') * INIT_ANG_STD;

if TESTING == 1
    UJI_Data_CPU = Get_Next_Step_2D_UJI('init', filename);
    if UJI_Data_CPU.num_steps == 0
        warning('在 %s 中未检测到步伐。跳过此文件。', filename);
        continue; 
    end
    NUM_STEPS = UJI_Data_CPU.num_steps + 1; 
    full_real_mag_history_CPU = UJI_Data_CPU.mag_readings;
    full_perfect_pdr_history_CPU = UJI_Data_CPU.pdr_steps;
else
    NUM_STEPS = NUM_STEPS_Fallback;
    full_real_mag_history_CPU = randn(NUM_STEPS, 3) * 50; 
    full_perfect_pdr_history_CPU = [ones(NUM_STEPS, 1)*0.7, zeros(NUM_STEPS, 1)];
end

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
fprintf('运行 %d 模拟步 (GPU-PF-DTW)...\n', NUM_STEPS - 1);
h_waitbar = waitbar(0, 'Running GPU-PF-DTW 2D Particle Filter...');

for t = 2:NUM_STEPS
   
    [true_state, pdr_step_simulated] = Get_Next_Step_2D_UJI(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step_simulated; 
    full_true_path_history(t, :) = true_state;
    
    data_end_idx = t - 1; 
    data_start_idx = max(1, data_end_idx - SEQUENCE_LEN + 1);
    
    % --- 核心修复：完全动态长度，彻底杜绝前置零填充 ---
    pdr_segment = full_perfect_pdr_history_CPU(data_start_idx:data_end_idx, :);
    pdr_history_for_function = single(pdr_segment);
    
    real_mag_data_segment = full_real_mag_history_CPU(data_start_idx:data_end_idx, :);
    real_mag_norms = sqrt(sum(real_mag_data_segment.^2, 2));
    live_sequence = real_mag_norms'; 
    live_sequence = live_sequence + randn(1, length(live_sequence)) * SENSOR_NOISE_STD; 
    % ----------------------------------------------------

    live_sequence_denoised = Denoise_Visushrink_CPU(single(live_sequence), 'db4');
    live_sequence_denoised_gpu = gpuArray(live_sequence_denoised);
    
    [particles_out_gpu, best_guess_gpu, dist] = Particle_Filter_DTW_Step_2D_CPFTest( ...
                                particles, ...
                                live_sequence_denoised_gpu, ...
                                pdr_history_for_function, ... 
                                geo_map, ... 
                                process_noise, DTW_NOISE_STD);
     
    M_current_step = size(particles_out_gpu, 1);
    N_reset = round(M_current_step * APFRATE);
    
    best_guess_state = gather(best_guess_gpu);
    
    if N_reset > 0
        p_cpu = randperm(M_current_step); 
        indices_cpu = p_cpu(1:N_reset);  
        indices_to_reset = gpuArray(indices_cpu); 
        particles_out_gpu(indices_to_reset, 1) = best_guess_state(1) + randn(N_reset, 1, 'single', 'gpuArray') * pos_std;
        particles_out_gpu(indices_to_reset, 2) = best_guess_state(2) + randn(N_reset, 1, 'single', 'gpuArray') * pos_std;
        particles_out_gpu(indices_to_reset, 3) = best_guess_state(3) + randn(N_reset, 1, 'single', 'gpuArray') * ang_std;
    end
    
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
                                    
    true_path_history(t, :) = true_state(1:2);
    estimated_path_history(t, :) = best_guess_state(1:2); 
    M_history(t) = current_M; 
    
    waitbar(t/NUM_STEPS, h_waitbar, sprintf('GPU-PF-DTW (M=%d, dist=%.1f)', current_M, dist_cpu));
end
close(h_waitbar);
fprintf('模拟完成。\n');

%% 5. PLOTTING
%==========================================================================
fprintf('正在从GPU收集数据用于绘图...\n');
particles_for_plot = gather(particles);
map_for_plot = gather(geo_map.Mag_map);
X_for_plot = gather(geo_map.X_grid);
Y_for_plot = gather(geo_map.Y_grid);

fprintf('正在绘图...\n');
figure('Name', char(filename), 'Position', [100, 100, 1600, 500]);

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

subplot(1, 3, 2); 
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('Localization Error');
xlabel('Time Step'); ylabel('Error'); grid on; ylim([0, Inf]);
mean_errors = mean(errors);
disp(["文件 " + filenames(i) + " 的平均误差为: ", mean_errors])

subplot(1, 3, 3);
plot(1:NUM_STEPS, M_history, 'b-', 'LineWidth', 2);
title('Adaptive Particle Count (M_t)');
xlabel('Time Step'); ylabel('Particle Count (M)'); grid on;
line([1, NUM_STEPS], [KLD.n_min, KLD.n_min], 'Color', 'red', 'LineStyle', '--');
line([1, NUM_STEPS], [KLD.n_max, KLD.n_max], 'Color', 'red', 'LineStyle', '--');
ylim([KLD.n_min*0.8, KLD.n_max*1.2]);

file_number_str = sprintf('%02d', i);
save_filename = "DDTWCPF_tt" + file_number_str + ".png";
fprintf('正在保存图像: %s\n', save_filename);
exportgraphics(gcf, save_filename, 'Resolution', 150);
    
end 

%% ========================================================================
%  == 子函数定义区 ==
%  ========================================================================
function [particles_out_gpu, best_guess_gpu, d_min_gpu] = Particle_Filter_DTW_Step_2D_CPFTest(particles_in_gpu, live_sequence_gpu, ...
                                            pdr_history_cpu, geo_map_gpu, process_noise, DTW_NOISE_STD)
    
    M = size(particles_in_gpu, 1);
    L = size(pdr_history_cpu, 1);
    
    last_pdr_step_cpu = single(pdr_history_cpu(end, :)); 
    
    step_noise_gpu = randn(M, 1, 'single', 'gpuArray') * single(process_noise.step_std);
    theta_noise_gpu = randn(M, 1, 'single', 'gpuArray') * single(process_noise.theta_std);
    
    steps_gpu = last_pdr_step_cpu(1) + step_noise_gpu;
    d_thetas_gpu = last_pdr_step_cpu(2) + theta_noise_gpu;
    
    old_theta_gpu = particles_in_gpu(:, 3);
    
    new_theta_gpu = old_theta_gpu + d_thetas_gpu;
    new_x_gpu = particles_in_gpu(:, 1) + sin(new_theta_gpu) .* steps_gpu;
    new_y_gpu = particles_in_gpu(:, 2) + cos(new_theta_gpu) .* steps_gpu;
    
    map_max_x = geo_map_gpu.X_grid(1,end);
    map_max_y = geo_map_gpu.Y_grid(end,1);
    new_x_gpu = max(1, min(map_max_x, new_x_gpu));
    new_y_gpu = max(1, min(map_max_y, new_y_gpu));
    
    propagated_particles_gpu = [new_x_gpu, new_y_gpu, new_theta_gpu];
    
    all_particle_paths_x_gpu = zeros(M, L, 'single', 'gpuArray');
    all_particle_paths_y_gpu = zeros(M, L, 'single', 'gpuArray');
    all_particle_thetas_gpu = zeros(M, L, 'single', 'gpuArray');
    
    all_particle_paths_x_gpu(:, L) = propagated_particles_gpu(:, 1);
    all_particle_paths_y_gpu(:, L) = propagated_particles_gpu(:, 2);
    all_particle_thetas_gpu(:, L) = propagated_particles_gpu(:, 3);
    
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
    
    all_map_sequences_gpu = interp2(geo_map_gpu.X_grid, geo_map_gpu.Y_grid, geo_map_gpu.Mag_map, ...
                                  all_particle_paths_x_gpu, all_particle_paths_y_gpu, 'linear', 0);
                              
    all_map_sequences_cpu = gather(all_map_sequences_gpu);
    live_sequence_cpu = gather(live_sequence_gpu);
    live_sequence_cpu = live_sequence_cpu(:)'; 
    
    all_distances_cpu = zeros(M, 1, 'single');
    
    parfor m = 1:M
        map_seq = all_map_sequences_cpu(m, :);
        distance = dtw(live_sequence_cpu, map_seq);
        all_distances_cpu(m) = distance;
    end
    
    min_dist_current_step = min(all_distances_cpu);
    disp(['当前帧最小DTW距离: ', num2str(min_dist_current_step)]); 
    
    relative_distances = all_distances_cpu - min_dist_current_step;
    
    alpha = 50; 
    weights_cpu = exp(-relative_distances / alpha);
    
    weights_cpu = double(weights_cpu) + 1e-12; 
    sum_weights = sum(weights_cpu);
    weights_cpu = weights_cpu / sum_weights;
    
    cdf_cpu = cumsum(weights_cpu);
    cdf_cpu(end) = 1.0; 
    
    r_0_cpu = rand(1, 'double') / M;
    U_cpu = r_0_cpu + (transpose(0:M-1)) / M;
    
    [~, indices_cpu] = histc(U_cpu, [0; cdf_cpu]);
    indices_cpu = indices_cpu + 1; 
    indices_cpu = min(indices_cpu, M); 
    
    indices_gpu = gpuArray(single(indices_cpu));
    particles_out_gpu = propagated_particles_gpu(indices_gpu, :);
    
    mean_pos_gpu = mean(particles_out_gpu(:, 1:2), 1);
    mean_sin_gpu = mean(sin(particles_out_gpu(:, 3)));
    mean_cos_gpu = mean(cos(particles_out_gpu(:, 3)));
    mean_theta_gpu = atan2(mean_sin_gpu, mean_cos_gpu);
    
    best_guess_gpu = [mean_pos_gpu, mean_theta_gpu];
    d_min_gpu = gpuArray(single(min_dist_current_step));
end

function particles_next_gpu = Adjust_Particle_Set_GPU(particles_curr_gpu, M_new_cpu)
    M_new = round(M_new_cpu);
    M_curr = size(particles_curr_gpu, 1);
    
    if M_new == M_curr
        particles_next_gpu = particles_curr_gpu;
        return;
    elseif M_new > M_curr
        N_add = M_new - M_curr;
        indices_gpu = randi(M_curr, N_add, 1, 'gpuArray'); 
        particles_to_add_gpu = particles_curr_gpu(indices_gpu, :);
        particles_next_gpu = [particles_curr_gpu; particles_to_add_gpu];
    else 
        p_cpu = randperm(M_curr); 
        indices_cpu = p_cpu(1:M_new); 
        indices_to_keep_gpu = gpuArray(indices_cpu); 
        particles_next_gpu = particles_curr_gpu(indices_to_keep_gpu, :);
    end
end

function denoised_data = Denoise_Visushrink_CPU(data, wavelet)
    n = length(data);
    wlevel = wmaxlev(n, wavelet);
    if wlevel < 1 
        denoised_data = data;
        return;
    end
    [c, l] = wavedec(data, wlevel, wavelet); 
    
    cD1_len = l(end-1);
    cD1_start_index = sum(l(1:end-2)) + 1;
    cD1_end_index = cD1_start_index + cD1_len - 1;
    
    if cD1_len >= 2 && cD1_end_index <= length(c)
        cD1 = c(cD1_start_index : cD1_end_index);
        sigma_est = mad(cD1, 1) / single(0.6745); 
        threshold = sigma_est * single(sqrt(2 * log(single(n))));
    else
        threshold = single(0); 
    end
    cA = c(1:l(1));
    detail_coeffs = c(l(1)+1 : end);
    denoised_details = wthresh(detail_coeffs, 's', threshold);
    denoised_c = [cA, denoised_details];
    
    denoised_data_temp = waverec(denoised_c, l, wavelet);
    if length(denoised_data_temp) > n
        denoised_data = denoised_data_temp(1:n);
    elseif length(denoised_data_temp) < n
        denoised_data = [denoised_data_temp, zeros(1, n-length(denoised_data_temp), 'single')];
    else
        denoised_data = denoised_data_temp;
    end
end