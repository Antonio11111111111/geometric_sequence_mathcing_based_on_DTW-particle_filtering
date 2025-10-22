function next_location = Get_Next_Step(last_location, map_len, margin)
 
        % 随机选择一个方向: -1 或 +1
        direction = randi(2) * 2 - 3; 

    
    next_location = last_location + direction;
    
    % 边界检查: 如果撞墙了，就朝反方向走
    if next_location <= margin || next_location >= (map_len - margin)
        next_location = last_location - direction;
    end
end