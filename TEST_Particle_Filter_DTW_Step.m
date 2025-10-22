% this is the test of Particle_Filter_DTW_Step function, which is 1
% dimensional
clc;
clear;
close all;


N_PARTICLES = 1000;       
NUM_STEPS = 100;          
WINDOW_RADIUS = 15;        
SENSOR_NOISE_STD = 0.5;   
DTW_NOISE_STD = 1.0;      
PROCESS_NOISE_STD = 1;   

map_len = 50;
geo_map = Geometric_Map_Generator(1, map_len);
fprintf('地图加载完成. 长度: %d\n', map_len);


fprintf('初始化仿真...\n');

true_location = randi([WINDOW_RADIUS+1, map_len-WINDOW_RADIUS]); % 随机开始 (避开边缘)

particles = randi(map_len, N_PARTICLES, 1); 

true_path_history = zeros(NUM_STEPS, 1);
estimated_path_history = zeros(NUM_STEPS, 1);
true_path_history(1) = true_location;
estimated_path_history(1) = mean(particles);


fprintf('运行 %d 步仿真...\n', NUM_STEPS);
h_waitbar = waitbar(0, '运行 粒子滤波+DTW...');

for t = 2:NUM_STEPS
   
    true_location = Get_Next_Step(true_location, map_len, WINDOW_RADIUS);

    live_sequence_indices = (true_location - WINDOW_RADIUS) : (true_location + WINDOW_RADIUS);
    live_sequence = geo_map(live_sequence_indices) + randn(1, 2*WINDOW_RADIUS+1) * SENSOR_NOISE_STD;

    

    
    [particles, best_guess] = Particle_Filter_DTW_Step(particles, live_sequence, ...
                                    geo_map, WINDOW_RADIUS, PROCESS_NOISE_STD, DTW_NOISE_STD);
    
    % --- [存储结果] ---
    true_path_history(t) = true_location;
    estimated_path_history(t) = best_guess;
    
    waitbar(t/NUM_STEPS, h_waitbar);
end

close(h_waitbar);
fprintf('仿真完成.\n');

figure;
plot(1:NUM_STEPS, true_path_history, 'b-o', 'LineWidth', 2.5, 'DisplayName', '真实路径');
hold on;
plot(1:NUM_STEPS, estimated_path_history, 'r--*', 'LineWidth', 1.5, 'DisplayName', 'PF+DTW 估计路径');
title('1D 粒子滤波 (PF+DTW) 跟踪');
xlabel('时间步 (Time Step)');
ylabel('地图位置 (Map Index)');
legend('show', 'Location', 'best');
grid on;
ylim([0, map_len + 1]); 

