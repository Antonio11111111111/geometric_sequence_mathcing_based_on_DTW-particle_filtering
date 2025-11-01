%==========================================================================
% SCRIPT: RUN_ROBUSTNESS_TEST
%
% DESCRIPTION:
%   执行鲁棒性分析实验。
%   本脚本循环运行一个简化的 MAIN 脚本, 每次使用不同的 DTW_NOISE_STD 值,
%   并记录 t=50 时的估计位置, 以测试算法对该参数的敏感性。
%
% 依赖:
%   - MAIN_2D_PFDTW.m (中的所有依赖, 如 Particle_Filter_DTW_Step_2D 等)
%==========================================================================
clc; clear; close all;

%% --- 1. 实验设置 ---
% 定义要测试的 DTW_NOISE_STD (sigma) 范围
std_values_to_test = [20, 23, 25, 30, 35, 40, 50, 60, 70, 80, 90];

% 固定的仿真步数
NUM_STEPS_FOR_TEST = 10; 

% 结果存储
results_dist_from_origin = zeros(length(std_values_to_test), 1);
true_dist_from_origin = 0; % 用于存储 t=50 时的真实距离

fprintf('开始执行鲁棒性分析 (共 %d 次运行)...\n', length(std_values_to_test));

%% --- 2. 生成一次性地图 (保持地图不变) ---
fprintf('生成共享地磁图...\n');
MAP_X_LEN = 50;
MAP_Y_LEN = 50;
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);
Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 
Mag = imgaussfilt(Mag_raw, 3.0);
geo_map.X_grid = X;
geo_map.Y_grid = Y;
geo_map.Mag_map = Mag;

%% --- 3. 循环执行仿真 ---
h_waitbar_main = waitbar(0, '正在测试鲁棒性...');

for k_run = 1:length(std_values_to_test)
    
    waitbar(k_run / length(std_values_to_test), h_waitbar_main, ...
        sprintf('运行 %d/%d (sigma = %.1f)', k_run, length(std_values_to_test), std_values_to_test(k_run)));
    
    % --- 3a. 重置环境 (关键!) ---
    
    % [关键!] 重置随机种子, 确保真实路径/噪声/初始化在每次运行时都相同
    rng(123); 
    
    % --- 3b. 加载 MAIN 脚本中的参数 (节选) ---
    
    % [被测参数]
    DTW_NOISE_STD = std_values_to_test(k_run); % <--- 这是我们循环的变量
    
    M_init = 500;
    NUM_STEPS = NUM_STEPS_FOR_TEST; % <--- 使用我们固定的步数
    SEQUENCE_LEN = 50;
    SENSOR_NOISE_STD = 0.5;
    process_noise.step_std = 0.5;
    process_noise.theta_std = deg2rad(10);
    pos_std = 5.0;
    ang_std = 0.2;
    KLD.n_min = 500;
    KLD.n_max = 100000;
    KLD.bin_size_xy = 0.5;
    KLD.epsilon = 0.05;
    KLD.delta = 0.01;
    KLD.DTW_THRESH_HIGH = 23.9;
    
    % --- 3c. 运行 MAIN 脚本中的初始化 (%% 3) ---
    % [第 72 行] 此处需要 MAP_X_LEN 和 MAP_Y_LEN
    true_state = [MAP_X_LEN/2, MAP_Y_LEN/4, deg2rad(45)];
    INIT_POS_STD = 5.0;
    INIT_ANG_STD = 0.5;
    particles = zeros(M_init, 3);
    particles(:, 1) = true_state(1) + randn(M_init, 1) * INIT_POS_STD;
    particles(:, 2) = true_state(2) + randn(M_init, 1) * INIT_POS_STD;
    particles(:, 3) = true_state(3) + randn(M_init, 1) * INIT_ANG_STD;
    
    full_true_path_history = zeros(NUM_STEPS, 3);
    full_pdr_step_history = zeros(NUM_STEPS, 2);
    full_true_path_history(1, :) = true_state;
    current_M = M_init;
    
    % --- 3d. 运行 MAIN 脚本中的仿真 (%% 4) ---
    for t = 2:NUM_STEPS
        [true_state, pdr_step] = Get_Next_Step_2D(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
        full_pdr_step_history(t, :) = pdr_step;
        full_true_path_history(t, :) = true_state;
        
        start_idx = max(1, t - SEQUENCE_LEN + 1);
        end_idx = t;
        pdr_history_for_function = zeros(SEQUENCE_LEN, 2);
        actual_len = end_idx - start_idx + 1;
        pdr_history_for_function(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);
        
        live_sequence = zeros(1, SEQUENCE_LEN);
        path_segment = zeros(SEQUENCE_LEN, 3);
        path_segment(end-actual_len+1:end, :) = full_true_path_history(start_idx:end_idx, :);
        for k = 1:SEQUENCE_LEN
            pos_x = path_segment(k, 1); pos_y = path_segment(k, 2);
            if pos_x ~= 0 || pos_y ~= 0
                live_sequence(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, pos_x, pos_y, 'linear', 0);
            end
        end
        live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD;
        
        [particles_out, best_guess_state, dist] = Particle_Filter_DTW_Step_2D(particles, live_sequence, ...
                                        pdr_history_for_function, geo_map, process_noise, DTW_NOISE_STD);
        
        M_current_step = size(particles_out, 1);
        N_reset = round(M_current_step * 0.05);
        indices_to_reset = randperm(M_current_step, N_reset);
        particles_out(indices_to_reset, 1) = best_guess_state(1) + randn(N_reset, 1) * pos_std;
        particles_out(indices_to_reset, 2) = best_guess_state(2) + randn(N_reset, 1) * pos_std;
        particles_out(indices_to_reset, 3) = best_guess_state(3) + randn(N_reset, 1) * ang_std;
        
        M_new = Adapt_Particle_Count_DTW_KLD(particles_out, dist, KLD.bin_size_xy, ...
                                        MAP_X_LEN, MAP_Y_LEN, ...
                                        KLD.epsilon, KLD.delta, ...
                                        KLD.n_min, KLD.n_max, ...
                                        KLD.DTW_THRESH_HIGH);
        if M_new ~= current_M
            particles_next = Adjust_Particle_Set(particles_out, M_new);
        else
            particles_next = particles_out;
        end
        particles = particles_next;
        current_M = size(particles, 1);
    end
    
    % --- 3e. 记录 t=50 (最后一步) 的结果 ---
    x_est = best_guess_state(1);
    y_est = best_guess_state(2);
    results_dist_from_origin(k_run) = sqrt(x_est^2 + y_est^2);
    
    % 记录一次 t=50 时的真实距离 (在所有运行中都相同)
    if k_run == 1
        x_true = true_state(1);
        y_true = true_state(2);
        true_dist_from_origin = sqrt(x_true^2 + y_true^2);
    end
    
    % [最终修正!] 清理变量, 防止泄漏到下一次循环
    % 将 MAP_X_LEN 和 MAP_Y_LEN 添加到 "例外" 列表中
    clearvars -except std_values_to_test results_dist_from_origin true_dist_from_origin k_run geo_map h_waitbar_main NUM_STEPS_FOR_TEST MAP_X_LEN MAP_Y_LEN
end
close(h_waitbar_main);
fprintf('鲁棒性分析完成。\n');

%% --- 4. 绘图 ---
figure('Position', [200, 200, 900, 500]);
plot(std_values_to_test, results_dist_from_origin, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
% line([min(std_values_to_test), max(std_values_to_test)], ...
%      [true_dist_from_origin, true_dist_from_origin], ...
%      'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
 
title('算法对 \sigma_{dtw} 参数的鲁棒性分析 (t=50)');
xlabel('DTW 噪声标准差 (\sigma_{dtw})');
ylabel('距离原点的欧氏距离 (sqrt(x_{est}² + y_{est}²))');
legend('估计距离 (Est. Distance)', 'Location', 'best');
grid on;
set(gca, 'FontSize', 12);