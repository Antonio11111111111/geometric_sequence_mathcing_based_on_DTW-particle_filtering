
clc;
clear;
close all;

%% 1. 仿真参数 (Parameters)
N_PARTICLES = 2000;       % 粒子数量
NUM_STEPS = 100;          % 仿真步数
SEQUENCE_LEN = 25;        % 序列长度

% 地图尺寸
MAP_X_LEN = 50;
MAP_Y_LEN = 50;

% 噪声参数
SENSOR_NOISE_STD = 0.5;   % 传感器噪声

% [修改 1: 大幅提高DTW噪声标准差]
DTW_NOISE_STD = 25;     

% [修改 2: 调大粒子过程噪声]
process_noise.step_std = 0.3;         % (大于真实的 0.2)
process_noise.theta_std = deg2rad(12.5); % (大于真实的 10.0 度)


%% 2. 生成2D地磁地图 (Map Generation)
fprintf('生成非对称平滑地磁地图...\n');

% 1. 生成坐标网格 (X, Y grids)
[X, Y] = meshgrid(1:MAP_X_LEN, 1:MAP_Y_LEN);

% 2. 调用你的 Generator 生成 "白噪声"
Mag_raw = Geometric_Map_Generator(2, [MAP_X_LEN, MAP_Y_LEN]); 

% 3. [!!] 对白噪声地图进行高斯平滑 [!!]
% 3.0 的平滑半径会创造一个非常平滑、简单的地图, 适合调试
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

particles = zeros(N_PARTICLES, 3);
% 1. 在真实 X (true_state(1)) 附近初始化
particles(:, 1) = true_state(1) + randn(N_PARTICLES, 1) * INIT_POS_STD;
% 2. 在真实 Y (true_state(2)) 附近初始化
particles(:, 2) = true_state(2) + randn(N_PARTICLES, 1) * INIT_POS_STD;
% 3. 在真实 Theta (true_state(3)) 附近初始化
particles(:, 3) = true_state(3) + randn(N_PARTICLES, 1) * INIT_ANG_STD;

% 历史记录 (用于绘图)
full_true_path_history = zeros(NUM_STEPS, 3); % 存储 [x, y, theta]
full_pdr_step_history = zeros(NUM_STEPS, 2);  % 存储 [step_len, d_theta]
true_path_history = zeros(NUM_STEPS, 2);      % 存储 [x, y] 用于绘图
estimated_path_history = zeros(NUM_STEPS, 2); % 存储 [x, y] 用于绘图

% 存储初始状态 (t=1)
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);
estimated_path_history(1, :) = true_path_history(1, :);


%% 4. 运行仿真 (Run Simulation)
fprintf('运行 %d 步仿真...\n', NUM_STEPS);
h_waitbar = waitbar(0, '运行 2D 粒子滤波+DTW (已修复)...');

for t = 2:NUM_STEPS
   
    % --- 4a. 模拟真实运动 (PDR) ---
    [true_state, pdr_step] = get_next_step_random(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    
    % 存储完整的历史记录
    full_pdr_step_history(t, :) = pdr_step;   
    full_true_path_history(t, :) = true_state; 
    
    % --- 4b. 准备函数输入 ---
    
    % i. 获取PDR历史 (最后 SEQUENCE_LEN 步)
    start_idx = max(1, t - SEQUENCE_LEN + 1);
    end_idx = t;
    
    pdr_history_for_function = zeros(SEQUENCE_LEN, 2);
    actual_len = end_idx - start_idx + 1;
    pdr_history_for_function(end-actual_len+1:end, :) = full_pdr_step_history(start_idx:end_idx, :);

    % ii. 生成 "Live" 传感器序列
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
    
    % 添加传感器噪声
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD;

    % --- 4c. 调用 *2D* 粒子滤波器 ---
    [particles, best_guess] = Particle_Filter_DTW_Step_2D(particles, live_sequence, ...
                                    pdr_history_for_function, geo_map, process_noise, DTW_NOISE_STD);
    
    % --- 4d. 存储结果 (用于绘图) ---
    true_path_history(t, :) = true_state(1:2); 
    estimated_path_history(t, :) = best_guess; 
    
    waitbar(t/NUM_STEPS, h_waitbar);
end

close(h_waitbar);
fprintf('仿真完成.\n');

%% 5. 绘图 (Plotting)
figure('Position', [100, 100, 1200, 500]);
% --- 图 1: 2D 路径图 ---
subplot(1, 2, 1);
imagesc(geo_map.X_grid(1,:), geo_map.Y_grid(:,1), geo_map.Mag_map);
hold on;
axis xy; 
colormap('jet');
colorbar;
title('2D 路径跟踪 (PF+DTW) - 随机路径');
xlabel('X 位置');
ylabel('Y 位置');
plot(particles(:, 1), particles(:, 2), 'k.', 'MarkerSize', 2, 'DisplayName', '最终粒子云');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', '真实路径');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'PF+DTW 估计');

% --- [新代码] 标注开始和终止位置 ---
% 标注开始位置 (绿色方块)
plot(true_path_history(1, 1), true_path_history(1, 2), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', '真实起点');
plot(estimated_path_history(1, 1), estimated_path_history(1, 2), 'gs', 'MarkerSize', 12, 'DisplayName', '估计起点');
% 标注终止位置 (红色X)
plot(true_path_history(end, 1), true_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', '真实终点');
plot(estimated_path_history(end, 1), estimated_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', '估计终点');
% --- 结束 [新代码] ---

legend('show', 'Location', 'best');
axis equal; 
axis([0 MAP_X_LEN 0 MAP_Y_LEN]);
% --- 图 2: 误差图 ---
subplot(1, 2, 2);
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('定位误差 (欧氏距离)');
xlabel('时间步 (Time Step)');
ylabel('误差 (Error in meters/units)');
grid on;
ylim([0, Inf]); % (修复 绘图Bug)


% =========================================================================
% 本地函数定义 (LOCAL FUNCTIONS)
% =========================================================================

function geometric_map = Geometric_Map_Generator(D, map_size)
% (*** 这是你提供的函数 ***)
% (*** 注意: 在当前版本中, 它未被调用, 因为我们需要一个平滑的地图 ***)
switch D
    case 1
        if numel(map_size) < 1
            error('For D=1, size_vec must have at least one element specifying the length.');
        end
        map_length = map_size(1);
        geometric_map = rand(1, map_length);
        
    case 2
        if numel(map_size) < 2
            error('For D=2, size_vec must have at least two elements for length and width.');
        end
        map_rows = map_size(1);
        map_cols = map_size(2);
        geometric_map = rand(map_rows, map_cols);
        
    otherwise
        error('Invalid dimension specified. D must be 1 or 2.');
end
end

% =========================================================================

function [new_state, pdr_step] = get_next_step_random(old_state, MAP_X_LEN, MAP_Y_LEN)
% 模拟一步 *随机* 的 2D PDR 运动 (定义 "真实" 噪声)

    % 1. 定义 "真实世界" 的PDR噪声
    MEAN_STEP_LEN = 1.0;    % 平均步长
    STEP_LEN_STD = 0.2;     % 真实步长噪声
    
    MEAN_TURN_RAD = 0.0;    % 平均转向 (0 = 倾向于走直线)
    TURN_STD_RAD = deg2rad(10.0); % 真实转向噪声 (10.0 度)

    % 2. 生成随机的 PDR 步骤
    random_step_len = MEAN_STEP_LEN + randn() * STEP_LEN_STD;
    random_d_theta = MEAN_TURN_RAD + randn() * TURN_STD_RAD;
    
    if random_step_len < 0, random_step_len = 0; end
    
    pdr_step = [random_step_len, random_d_theta];
    
    % 4. 应用PDR运动模型
    old_x = old_state(1);
    old_y = old_state(2);
    old_theta = old_state(3);
    
    new_theta = old_theta + random_d_theta;
    new_x = old_x + sin(new_theta) * random_step_len;
    new_y = old_y + cos(new_theta) * random_step_len;
    
    % 5. 边界检查
    new_x = max(1, min(MAP_X_LEN, new_x));
    new_y = max(1, min(MAP_Y_LEN, new_y));
    
    new_state = [new_x, new_y, new_theta];
end

% =========================================================================
