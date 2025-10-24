function [new_state, pdr_step] = get_next_step_random(old_state, MAP_X_LEN, MAP_Y_LEN)

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
