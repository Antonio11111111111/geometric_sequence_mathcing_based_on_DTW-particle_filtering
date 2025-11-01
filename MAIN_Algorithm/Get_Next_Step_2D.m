function [new_state, pdr_step] = Get_Next_Step_2D(old_state, MAP_X_LEN, MAP_Y_LEN)
% 模拟一步 *随机* 的 2D PDR 运动 (随机行走模型)
%
%   old_state:  [1x3] 向量 [x, y, theta]
%   MAP_X_LEN:  地图 X 边界
%   MAP_Y_LEN:  地图 Y 边界
%
%   new_state:  [1x3] 向量 [new_x, new_y, new_theta]
%   pdr_step:   [1x2] 向量 [random_step_len, random_d_theta]

    % 1. 定义随机行走的参数
    % 你可以调整这些值来改变 "随机" 的程度
    
    MEAN_STEP_LEN = 1.0;    % 平均步长 (单位/米)
    STEP_LEN_STD = 0.2;     % 步长标准差 (步长波动的幅度)
    
    MEAN_TURN_RAD = 0.0;    % 平均转向 (0 = 倾向于走直线)
    TURN_STD_RAD = deg2rad(10.0); % 转向标准差 (每次转向波动的幅度, e.g., +/- 10度)

    % 2. 生成随机的 PDR 步骤
    % randn() 生成符合高斯分布(正态分布)的随机数
    random_step_len = MEAN_STEP_LEN + randn() * STEP_LEN_STD;
    random_d_theta = MEAN_TURN_RAD + randn() * TURN_STD_RAD;
    
    % (可选) 确保步长不会是负数
    if random_step_len < 0
        random_step_len = 0;
    end
    
    % 3. 存储这个随机生成的PDR步骤 (用于粒子滤波器的历史记录)
    pdr_step = [random_step_len, random_d_theta];
    
    % 4. 应用PDR运动模型
    old_x = old_state(1);
    old_y = old_state(2);
    old_theta = old_state(3);
    
    new_theta = old_theta + random_d_theta;
    
    % [sin(theta_new); cos(theta_new)]
    new_x = old_x + sin(new_theta) * random_step_len;
    new_y = old_y + cos(new_theta) * random_step_len;
    
    % 5. 边界检查
    new_x = max(1, min(MAP_X_LEN, new_x));
    new_y = max(1, min(MAP_Y_LEN, new_y));
    
    new_state = [new_x, new_y, new_theta];
end