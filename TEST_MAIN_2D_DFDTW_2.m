clc;
clear;
close all;
%% 1. 仿真参数 (Parameters)
% [修改 1: 粒子数参数变为自适应]
INITIAL_N_PARTICLES = 100;      % 初始粒子数量
M_MIN = 100;                     % 最小粒子数 (可以设更低)
M_MAX = 10000;                   % 最大粒子数 (可以设更高)
% [新参数: 论文 (Elvira et al.) 中的收敛性评估参数]
K_SAMPLES = 15;                   % (K) 虚拟观测样本数
WINDOW_SIZE = 50;                % (W) 统计窗口大小
P_LOW = 0.9;                     % (p_l) 低p-value阈值
P_HIGH = 1;                    % (p_h) 高p-value阈值
% [新参数: 强制多样化的抖动强度]
JITTER_POS_STD = 0.05;           % 位置抖动标准差 (单位: 米/单位)
JITTER_ANG_STD = deg2rad(0.5);   % 角度抖动标准差 (单位: 弧度)
% --- 仿真和地图 ---
NUM_STEPS = 100;                 % 仿真步数
SEQUENCE_LEN = 50;               % 序列长度
MAP_X_LEN = 50;
MAP_Y_LEN = 50;
% --- 噪声 ---
SENSOR_NOISE_STD = 0.5;          % 传感器噪声 (用于生成 live_sequence 和 虚拟观测)
DTW_NOISE_STD = 30;              % (这个参数现在只被 DTW 函数使用)
process_noise.step_std = 0.25;    % 粒子过程噪声 (大于真实的 0.2)
process_noise.theta_std = deg2rad(10.5); % (大于真实的 10.0 度)

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
% [修改 2: M_t 是动态的]
M_t = INITIAL_N_PARTICLES;  % M_t 是当前的粒子数
% [新: 用于评估的状态变量]
history_A_K = nan(WINDOW_SIZE, 1);
% 粒子: [Nx3] 矩阵, 每一行是 [x, y, theta]
INIT_POS_STD = 2.0;  
INIT_ANG_STD = 0.5;  
particles = zeros(M_t, 3);
particles(:, 1) = true_state(1) + randn(M_t, 1) * INIT_POS_STD;
particles(:, 2) = true_state(2) + randn(M_t, 1) * INIT_POS_STD;
particles(:, 3) = true_state(3) + randn(M_t, 1) * INIT_ANG_STD;
% 历史记录 (用于绘图)
full_true_path_history = zeros(NUM_STEPS, 3); 
full_pdr_step_history = zeros(NUM_STEPS, 2);  
true_path_history = zeros(NUM_STEPS, 2);      
estimated_path_history = zeros(NUM_STEPS, 2); 
full_M_history = zeros(NUM_STEPS, 1); 
% 存储初始状态 (t=1)
full_true_path_history(1, :) = true_state;
true_path_history(1, :) = true_state(1:2);
estimated_path_history(1, :) = true_path_history(1, :);
full_M_history(1) = M_t;

%% 4. 运行仿真 (Run Simulation)
fprintf('运行 %d 步仿真...\n', NUM_STEPS);
h_waitbar = waitbar(0, '运行 2D 粒子滤波+DTW (自适应 M, 强制多样化)...');

for t = 2:NUM_STEPS
   
    % --- 4a. 模拟真实运动 (PDR) ---
    [true_state, pdr_step] = get_next_step_random(full_true_path_history(t-1, :), MAP_X_LEN, MAP_Y_LEN);
    full_pdr_step_history(t, :) = pdr_step;   
    full_true_path_history(t, :) = true_state; 
    
    % --- 4b. 准备函数输入 ---
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
        
        if pos_x == 0 && pos_y == 0 
            live_sequence(k) = 0;
        else
            live_sequence(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                                       pos_x, pos_y, 'linear', 0);
        end
    end
    live_sequence = live_sequence + randn(1, SEQUENCE_LEN) * SENSOR_NOISE_STD;
    
    % --- [新: 4c. 自适应 M 评估 (基于论文) - 预测步骤] ---
    y_t_real = live_sequence(end);
    pdr_step_t = pdr_history_for_function(end, :);
    propagated_particles_for_eval = propagate_particles_manual(particles, pdr_step_t, process_noise, MAP_X_LEN, MAP_Y_LEN);
    
    % --- [新: 4d. 生成统计量 A_K,M,t] ---
    a_k_m_t = Generate_A_K_Statistic(propagated_particles_for_eval, y_t_real, geo_map, SENSOR_NOISE_STD, K_SAMPLES);
    
    % --- [新: 4e. 更新统计量历史窗口] ---
    history_A_K(mod(t-2, WINDOW_SIZE) + 1) = a_k_m_t;
    
    % --- [4f. 调用 *2D* 粒子滤波器 (主估计器)] ---
    % !! 注意: 现在调用的是修改后的、包含强制多样化的版本 !!
    [particles_out, best_guess, d_min] = Particle_Filter_DTW_Step_2D(particles, live_sequence, ...
                                    pdr_history_for_function, geo_map, process_noise, DTW_NOISE_STD, ...
                                    JITTER_POS_STD, JITTER_ANG_STD); % <-- 传入抖动参数
    
    % --- [4g. 存储结果 (用于绘图)] ---
    true_path_history(t, :) = true_state(1:2); 
    estimated_path_history(t, :) = best_guess(1:2);
    
    % --- [新: 4h. 执行自适应 M 调整 (基于论文)] ---
    if (t > WINDOW_SIZE) && (mod(t, WINDOW_SIZE) == 0) % 每 W 步评估一次
        % (1) 调用评估函数 (使用更大胆的调整逻辑)
        M_t = Adapt_Particle_Count(size(particles,1), history_A_K, K_SAMPLES, WINDOW_SIZE, P_LOW, P_HIGH, M_MIN, M_MAX);
        
        % (2) 调整粒子集大小 (为下一次迭代准备)
        particles = Adjust_Particle_Set(particles_out, M_t);
        
    else
        % 如果不调整, 就使用上一轮的输出
        particles = particles_out;
    end
    
    % --- [新: 4i. 记录M] ---
    full_M_history(t) = size(particles, 1); 
    
    waitbar(t/NUM_STEPS, h_waitbar);
end
close(h_waitbar);
fprintf('仿真完成.\n');
fprintf('平均粒子数: %.0f\n', mean(full_M_history(2:end)));

%% 5. 绘图 (Plotting)
figure('Position', [100, 100, 1800, 500]);
% --- 图 1: 2D 路径图 ---
subplot(1, 3, 1);
imagesc(geo_map.X_grid(1,:), geo_map.Y_grid(:,1), geo_map.Mag_map);
hold on; axis xy; colormap('jet'); colorbar;
title('2D 路径跟踪 (自适应 M + DTW + 强制多样化)');
xlabel('X 位置'); ylabel('Y 位置');
plot(particles(:, 1), particles(:, 2), 'k.', 'MarkerSize', 2, 'DisplayName', '最终粒子云');
plot(true_path_history(:, 1), true_path_history(:, 2), 'b-o', 'LineWidth', 2.5, 'DisplayName', '真实路径');
plot(estimated_path_history(:, 1), estimated_path_history(:, 2), 'r--*', 'LineWidth', 1.5, 'DisplayName', 'PF+DTW 估计');
plot(true_path_history(1, 1), true_path_history(1, 2), 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', '真实起点');
plot(estimated_path_history(1, 1), estimated_path_history(1, 2), 'gs', 'MarkerSize', 12, 'DisplayName', '估计起点');
plot(true_path_history(end, 1), true_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', '真实终点');
plot(estimated_path_history(end, 1), estimated_path_history(end, 2), 'rx', 'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', '估计终点');
legend('show', 'Location', 'best'); axis equal; axis([0 MAP_X_LEN 0 MAP_Y_LEN]);

% --- 图 2: 误差图 ---
subplot(1, 3, 2);
errors = sqrt(sum((true_path_history - estimated_path_history).^2, 2));
plot(1:NUM_STEPS, errors, 'k-', 'LineWidth', 1.5);
title('定位误差 (欧氏距离)'); xlabel('时间步 (Time Step)'); ylabel('误差 (Error in meters/units)');
grid on; ylim([0, Inf]); 

% --- 图 3: 粒子数 M 变化图 ---
subplot(1, 3, 3);
plot(1:NUM_STEPS, full_M_history, 'b-', 'LineWidth', 2);
hold on;
line([1, NUM_STEPS], [M_MIN, M_MIN], 'Color', 'r', 'LineStyle', '--');
line([1, NUM_STEPS], [M_MAX, M_MAX], 'Color', 'r', 'LineStyle', '--');
title('自适应粒子数 (M_t) 变化'); xlabel('时间步 (Time Step)'); ylabel('粒子数量 (M_t)');
legend('M_t', 'Min/Max Bounds', 'Location', 'best'); grid on; ylim([0, M_MAX * 1.1]);

% =========================================================================
% 本地函数定义 (LOCAL FUNCTIONS)
% =========================================================================
function geometric_map = Geometric_Map_Generator(D, map_size)
% [来自用户代码]
switch D, case 1, if numel(map_size)<1, error(''); end, map_length=map_size(1); geometric_map=rand(1,map_length); case 2, if numel(map_size)<2, error(''); end, map_rows=map_size(1); map_cols=map_size(2); geometric_map=rand(map_rows,map_cols); otherwise, error(''); end
end
% =========================================================================
function [new_state, pdr_step] = get_next_step_random(old_state, MAP_X_LEN, MAP_Y_LEN)
% [来自用户代码]
MEAN_STEP_LEN=1.0; STEP_LEN_STD=0.2; MEAN_TURN_RAD=0.0; TURN_STD_RAD=deg2rad(10.0); 
random_step_len=MEAN_STEP_LEN+randn()*STEP_LEN_STD; random_d_theta=MEAN_TURN_RAD+randn()*TURN_STD_RAD;
if random_step_len<0, random_step_len=0; end; pdr_step=[random_step_len, random_d_theta];
old_x=old_state(1); old_y=old_state(2); old_theta=old_state(3);
new_theta=old_theta+random_d_theta; new_x=old_x+sin(new_theta)*random_step_len; new_y=old_y+cos(new_theta)*random_step_len;
new_x=max(1,min(MAP_X_LEN,new_x)); new_y=max(1,min(MAP_Y_LEN,new_y)); new_state=[new_x,new_y,new_theta];
end
% =========================================================================
function [particles_out, best_guess, d_min] = Particle_Filter_DTW_Step_2D(particles_in, live_sequence, ...
                                            pdr_history, geo_map, process_noise, DTW_NOISE_STD, ...
                                            JITTER_POS_STD, JITTER_ANG_STD) % <-- 新增抖动参数
% [用户代码 + 修改]
% 摘要: 执行一步完整的“传播-加权-重采样-强制多样化”

    M = size(particles_in, 1);
    L = size(pdr_history, 1);
    
    % --- 1. Propagation (Prediction) ---
    last_pdr_step = pdr_history(end, :); 
    step_noise = randn(M, 1) * process_noise.step_std;
    theta_noise = randn(M, 1) * process_noise.theta_std;
    steps = last_pdr_step(1) + step_noise;
    d_thetas = last_pdr_step(2) + theta_noise;
    old_x = particles_in(:, 1); old_y = particles_in(:, 2); old_theta = particles_in(:, 3);
    new_theta = old_theta + d_thetas;
    new_x = old_x + sin(new_theta) .* steps; 
    new_y = old_y + cos(new_theta) .* steps; 
    new_x = max(1, min(geo_map.X_grid(1,end), new_x));
    new_y = max(1, min(geo_map.Y_grid(end,1), new_y));
    propagated_particles = [new_x, new_y, new_theta];
    
    % --- 2. Weighting ---
    weights = zeros(M, 1);
    dtw_variance = DTW_NOISE_STD^2;
    all_distances = zeros(M, 1); 
    
    for m = 1:M
        % --- 2a. Path Reconstruction ---
        particle_path = zeros(L, 2); 
        current_pos = propagated_particles(m, 1:2); current_theta = propagated_particles(m, 3);
        particle_path(end, :) = current_pos; 
        for k = L-1:-1:1
            pdr_step_len = pdr_history(k+1, 1); pdr_d_theta = pdr_history(k+1, 2);
            prev_theta = current_theta - pdr_d_theta;
            prev_x = current_pos(1) - sin(current_theta) * pdr_step_len;
            prev_y = current_pos(2) - cos(current_theta) * pdr_step_len;
            particle_path(k, :) = [prev_x, prev_y];
            current_pos = [prev_x, prev_y]; current_theta = prev_theta;
        end
        % --- 2b. Generate Map Sequence ---
        map_sequence = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                              particle_path(:, 1), particle_path(:, 2), 'linear', 0);
        % --- 2c. Calculate DTW Distance and Convert to Weight ---
        if exist('dtw', 'file')
             distance = dtw(live_sequence(:), map_sequence(:));
        else
             distance = dtw_cost(live_sequence(:), map_sequence(:));
             % if m == 1 && evalin('base','t') == 2, warning('Using local dtw_cost'); end % 避免重复警告
        end
        all_distances(m) = distance;
        weights(m) = exp(-distance^2 / (2 * dtw_variance));
    end
    % --- 2d. Normalize Weights ---
    sum_weights = sum(weights);
    if sum_weights > 1e-15, weights = weights / sum_weights; else weights = ones(M, 1) / M; end
    
    % --- 3. Resampling ---
    resampled_particles = zeros(M, 3); % <-- 先存到临时变量
    cdf = cumsum(weights); 
    r_0 = rand() / M; 
    idx_j = 1;
    for m = 1:M
        U = r_0 + (m-1)/M; 
        while U > cdf(idx_j), idx_j = idx_j + 1; end
        resampled_particles(m, :) = propagated_particles(idx_j, :);
    end

    % --- [新] 3b. Forced Diversification (Jittering) ---
    % 对所有重采样后的粒子添加少量噪声
    jitter_pos = randn(M, 2) * JITTER_POS_STD;
    jitter_ang = randn(M, 1) * JITTER_ANG_STD;
    
    particles_out = resampled_particles; % <-- 最终输出
    particles_out(:, 1:2) = particles_out(:, 1:2) + jitter_pos;
    particles_out(:, 3) = particles_out(:, 3) + jitter_ang;
    
    % 再次进行边界检查 (因为抖动可能导致出界)
    particles_out(:, 1) = max(1, min(geo_map.X_grid(1,end), particles_out(:, 1)));
    particles_out(:, 2) = max(1, min(geo_map.Y_grid(end,1), particles_out(:, 2)));

    % --- 4. Estimation ---
    best_guess(1:2) = mean(particles_out(:, 1:2), 1);
    best_guess(3) = atan2(mean(sin(particles_out(:, 3))), mean(cos(particles_out(:, 3))));
    d_min = min(all_distances);
end
% =========================================================================
% =========================================================================
% == 新增函数: 基于论文 (Elvira et al.) 的自适应逻辑 ==
% =========================================================================
function [propagated_particles] = propagate_particles_manual(particles_in, pdr_step, process_noise, MAP_X_LEN, MAP_Y_LEN)
% [新增函数]
M=size(particles_in,1); step_noise=randn(M,1)*process_noise.step_std; theta_noise=randn(M,1)*process_noise.theta_std;
steps=pdr_step(1)+step_noise; d_thetas=pdr_step(2)+theta_noise;
old_x=particles_in(:,1); old_y=particles_in(:,2); old_theta=particles_in(:,3);
new_theta=old_theta+d_thetas; new_x=old_x+sin(new_theta).*steps; new_y=old_y+cos(new_theta).*steps;
new_x=max(1,min(MAP_X_LEN,new_x)); new_y=max(1,min(MAP_Y_LEN,new_y)); propagated_particles=[new_x,new_y,new_theta];
end
% =========================================================================
function [a_k_m_t] = Generate_A_K_Statistic(propagated_particles, y_t_real, geo_map, SENSOR_NOISE_STD, K)
% [新增函数] (Algorithm 2)
M_t=size(propagated_particles,1); indices=randi(M_t,K,1);
predicted_means=interp2(geo_map.X_grid,geo_map.Y_grid,geo_map.Mag_map,propagated_particles(indices,1),propagated_particles(indices,2),'linear',0);
fictitious_obs=predicted_means+randn(K,1)*SENSOR_NOISE_STD; a_k_m_t=sum(fictitious_obs<y_t_real);
end
% =========================================================================
function [new_M] = Adapt_Particle_Count(current_M, history_S_t, K, W, p_low, p_high, M_min, M_max)
% [新增函数 + 修改] (Algorithm 3)
O_j=zeros(K+1,1); for j=0:K, O_j(j+1)=sum(history_S_t==j); end; E_j=W/(K+1);
chi2_stat=sum((O_j-E_j).^2/(E_j+1e-9)); p_value=1-chi2cdf(chi2_stat,K);
if p_value<=p_low
    % [修改] 更大胆的增加: * 3
    new_M=min(current_M*3, M_max); 
    fprintf(' (t=%d, p=%.2f <= %.2f) 增加 M -> %d\n', evalin('base','t'), p_value, p_low, new_M);
elseif p_value>=p_high
    % [修改] 更大胆的减少: / 3
    new_M=max(current_M/3, M_min); 
    fprintf(' (t=%d, p=%.2f >= %.2f) 减少 M -> %d\n', evalin('base','t'), p_value, p_high, new_M);
else, new_M=current_M; 
    % fprintf(' (t=%d, p=%.2f) 保持 M = %d\n', evalin('base','t'), p_value, new_M); 
end; new_M=round(new_M);
end
% =========================================================================
function [new_particle_set] = Adjust_Particle_Set(old_particle_set, new_M)
% [新增函数] (Algorithm 4)
current_M=size(old_particle_set,1); if new_M==current_M, new_particle_set=old_particle_set; return; end
if new_M>current_M, num_to_add=new_M-current_M; indices_to_add=randi(current_M,num_to_add,1); new_particles=old_particle_set(indices_to_add,:);
jitter_pos=randn(num_to_add,2)*0.1; jitter_ang=randn(num_to_add,1)*0.01; new_particles(:,1:2)=new_particles(:,1:2)+jitter_pos; new_particles(:,3)=new_particles(:,3)+jitter_ang;
new_particle_set=[old_particle_set; new_particles];
else, indices_to_keep=randperm(current_M,new_M); new_particle_set=old_particle_set(indices_to_keep,:); end
end
% =========================================================================
function cost = dtw_cost(s1, s2)
% [新增函数] (PF-DTW 的备用辅助函数)
n=length(s1); m=length(s2); D=zeros(n+1,m+1); D(1:end,1)=Inf; D(1,1:end)=Inf; D(1,1)=0; s1=s1(:); s2=s2(:);
for i=2:n+1, for j=2:m+1, cost=(s1(i-1)-s2(j-1)).^2; D(i,j)=cost+min([D(i-1,j),D(i,j-1),D(i-1,j-1)]); end; end; cost=sqrt(D(n+1,m+1));
end
% =========================================================================