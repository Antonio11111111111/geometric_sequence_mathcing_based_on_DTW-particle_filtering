% % function [new_state, pdr_step] = Get_Next_Step_2D_UJI(old_state, MAP_X_LEN, MAP_Y_LEN)
% % % Get_Next_Step_2D_UJI (V3 - 基于地面真实坐标点的版本)
% % %
% % %   此版本结合了 UJI 文件的两个部分：
% % %   1. 传感器数据 (Part 1): 用于检测 *何时* 发生步伐 (findpeaks)。
% % %   2. 路径点数据 (Part 2): 用于计算每个走廊的 *完美角度*。
% % %
% % %   结果是一个 "完美" 的 PDR 步骤列表，由直线 (d_theta=0) 和
% % %   瞬时转弯 (d_theta=完美角度差) 组成。
% % %
% % 
% % % --- [ PDR 调整参数 ] ---
% % persistent PDR_PARAMS
% % if isempty(PDR_PARAMS)
% %     PDR_PARAMS.MIN_PEAK_HEIGHT = 0.8; % 最小峰值高度 (m/s^2)
% %     PDR_PARAMS.MIN_PEAK_DIST = 4;   % 最小峰值距离 (0.4s)
% %     PDR_PARAMS.FIXED_STEP_LEN = 0.7;% (米) 假设的固定步长
% % end
% % % --------------------------
% % 
% % persistent pdr_steps_list step_counter
% % 
% % %==========================================================================
% % % 1. 初始化模式 (Initialization Mode)
% % %==========================================================================
% % if ischar(old_state) && strcmp(old_state, 'init')
% %     filename = MAP_X_LEN; 
% %     fprintf('Initializing UJI PDR (V3 - Ground Truth)... Loading file: %s\n', filename);
% % 
% %     % 加载并处理 UJI 文件 (V3)
% %     [pdr_steps_list, num_detected] = load_and_process_uji_file_v3(filename, PDR_PARAMS);
% % 
% %     step_counter = 1; 
% %     new_state = num_detected; % 这是返回给主脚本的步数
% %     pdr_step = [];    
% % 
% %     fprintf('UJI PDR initialized. %d "perfect" steps synthesized.\n', new_state);
% % 
% %     if new_state == 0
% %         warning('!!! V3: 检测到 0 步. 模拟循环将不会运行。');
% %         warning('!!! 请尝试在 Get_Next_Step_2D_UJI.m 中调低 PDR_PARAMS。');
% %     end
% %     return;
% % end
% % 
% % %==========================================================================
% % % 2. 执行模式 (Execution Mode)
% % %==========================================================================
% % if isempty(pdr_steps_list)
% %     error('UJI PDR not initialized. Call Get_Next_Step_2D_UJI(''init'', ''filename.txt'') before the simulation loop.');
% % end
% % 
% % % 从预先计算好的 "完美" 列表中获取下一步
% % if step_counter > size(pdr_steps_list, 1)
% %     pdr_step = [0, 0]; 
% % else
% %     pdr_step = pdr_steps_list(step_counter, :);
% %     step_counter = step_counter + 1;
% % end
% % 
% % %==========================================================================
% % % 3. 应用PDR运动模型 (无变化)
% % %==========================================================================
% % step_len = pdr_step(1);
% % d_theta = pdr_step(2);
% % old_x = old_state(1);
% % old_y = old_state(2);
% % old_theta = old_state(3);
% % new_theta = old_theta + d_theta;
% % new_x = old_x + sin(new_theta) * step_len;
% % new_y = old_y + cos(new_theta) * step_len;
% % new_x = max(1, min(MAP_X_LEN, new_x));
% % new_y = max(1, min(MAP_Y_LEN, new_y));
% % new_state = [new_x, new_y, new_theta];
% % end
% % 
% % 
% % %==========================================================================
% % % 内部辅助函数 (Local Helper Function) - (V3)
% % %==========================================================================
% % function [pdr_steps_list, N_total_steps] = load_and_process_uji_file_v3(filename, PDR_PARAMS)
% % 
% %     fid = fopen(filename);
% %     if fid == -1, error('Cannot open file: %s', filename); end
% % 
% %     % --- 1. 读取 Part 1 (传感器数据) ---
% %     sensor_lines = {};
% %     while ~feof(fid)
% %         line = fgetl(fid);
% %         if isempty(line) || startsWith(strtrim(line), '<'), break; end
% %         sensor_lines{end+1} = line;
% %     end
% %     sensor_data = cell2mat(cellfun(@(x) sscanf(x, '%f')', sensor_lines, 'UniformOutput', false)');
% %     if isempty(sensor_data)
% %         error('No sensor data found in file: %s', filename); 
% %     end
% % 
% %     % --- 2. 检测步伐 (来自传感器) ---
% %     accel_x = sensor_data(:, 5);
% %     accel_y = sensor_data(:, 6);
% %     accel_z = sensor_data(:, 7);
% %     mag_accel = sqrt(accel_x.^2 + accel_y.^2 + accel_z.^2);
% % 
% %     [pks, step_indices] = findpeaks(mag_accel, 'MinPeakHeight', PDR_PARAMS.MIN_PEAK_HEIGHT, 'MinPeakDistance', PDR_PARAMS.MIN_PEAK_DIST);
% % 
% %     N_total_steps = length(step_indices);
% %     if N_total_steps == 0
% %         pdr_steps_list = zeros(0, 2); % 返回空矩阵
% %         fclose(fid);
% %         return;
% %     end
% % 
% %     % --- 3. 读取 Part 2 (路径点数据) ---
% %     % (此时, fid 停在 <m> 标记行)
% %     if feof(fid)
% %         warning('未找到路径点数据 (Part 2)，将退回到 V2（嘈杂的）PDR 逻辑。');
% %         fclose(fid);
% %         % [回退逻辑 - 从 V2 复制而来]
% %         pdr_steps_list = zeros(N_total_steps, 2);
% %         azimuth_rad = deg2rad(sensor_data(:, 8));
% %         last_step_azimuth = azimuth_rad(step_indices(1));
% %         pdr_steps_list(1, 1) = PDR_PARAMS.FIXED_STEP_LEN;
% %         pdr_steps_list(1, 2) = 0;
% %         for i = 2:N_total_steps
% %             current_step_azimuth = azimuth_rad(step_indices(i));
% %             d_theta = current_step_azimuth - last_step_azimuth;
% %             pdr_steps_list(i, 1) = PDR_PARAMS.FIXED_STEP_LEN;
% %             pdr_steps_list(i, 2) = atan2(sin(d_theta), cos(d_theta));
% %             last_step_azimuth = current_step_azimuth;
% %         end
% %         return;
% %     end
% % 
% %     num_segments_str = strtrim(line);
% %     m = sscanf(num_segments_str, '<%d>');
% % 
% %     segment_data = zeros(m, 6); % lat1, lon1, lat2, lon2, FS, LS
% %     for i = 1:m
% %         line = fgetl(fid);
% %         segment_data(i, :) = sscanf(line, '%f')';
% %     end
% %     fclose(fid);
% % 
% %     % --- 4. 计算每个走廊的"完美"角度 ---
% %     segment_headings = zeros(m, 1);
% %     for i = 1:m
% %         lat1 = segment_data(i, 1);
% %         lon1 = segment_data(i, 2);
% %         lat2 = segment_data(i, 3);
% %         lon2 = segment_data(i, 4);
% %         segment_headings(i) = calculate_bearing(lat1, lon1, lat2, lon2);
% %     end
% % 
% %     % --- 5. 合成 "完美" PDR 步骤列表 ---
% %     pdr_steps_list = zeros(N_total_steps, 2);
% %     step_list_counter = 1; % 跟踪我们填充到 pdr_steps_list 的位置
% % 
% %     % 我们假设第一个走廊的角度是我们的 "初始" 角度
% %     last_heading = segment_headings(1); 
% % 
% %     for i = 1:m % 遍历每个走廊 (Segment)
% %         current_heading = segment_headings(i);
% % 
% %         % UJI 索引从 0 开始, MATLAB 索引从 1 开始
% %         segment_start_index = segment_data(i, 5) + 1;
% %         segment_end_index = segment_data(i, 6) + 1;
% % 
% %         % 找到所有在此走廊中发生的步伐
% %         steps_in_this_segment_mask = (step_indices >= segment_start_index & step_indices <= segment_end_index);
% %         num_steps_in_segment = sum(steps_in_this_segment_mask);
% % 
% %         if num_steps_in_segment == 0
% %             continue; % 这个走廊很短，没有检测到步伐
% %         end
% % 
% %         % --- 这是关键逻辑 ---
% %         % 计算从上一个走廊到这个走廊的转向角
% %         d_theta = current_heading - last_heading;
% %         d_theta = atan2(sin(d_theta), cos(d_theta)); % 环绕 -pi 到 +pi
% % 
% %         % 1. 此走廊的 *第一步* 承担所有转向
% %         pdr_steps_list(step_list_counter, 1) = PDR_PARAMS.FIXED_STEP_LEN;
% %         pdr_steps_list(step_list_counter, 2) = d_theta;
% %         step_list_counter = step_list_counter + 1;
% % 
% %         % 2. 此走廊的 *所有后续步伐* 都是直线
% %         for j = 2:num_steps_in_segment
% %             if step_list_counter > N_total_steps, break; end
% %             pdr_steps_list(step_list_counter, 1) = PDR_PARAMS.FIXED_STEP_LEN;
% %             pdr_steps_list(step_list_counter, 2) = 0; % 完美直线
% %             step_list_counter = step_list_counter + 1;
% %         end
% % 
% %         % 更新 "上一个" 角度，为下一次循环做准备
% %         last_heading = current_heading;
% %     end
% % end
% % 
% % 
% % function bearing_rad = calculate_bearing(lat1, lon1, lat2, lon2)
% % % 计算从 (lat1, lon1) 到 (lat2, lon2) 的方位角 (bearing)
% % % 角度以弧度为单位, 0 = 北, +pi/2 = 东
% % 
% %     lat1r = deg2rad(lat1);
% %     lon1r = deg2rad(lon1);
% %     lat2r = deg2rad(lat2);
% %     lon2r = deg2rad(lon2);
% % 
% %     dLon = lon2r - lon1r;
% % 
% %     y = sin(dLon) * cos(lat2r);
% %     x = cos(lat1r) * sin(lat2r) - sin(lat1r) * cos(lat2r) * cos(dLon);
% % 
% %     bearing_rad = atan2(y, x);
% % end
% 
% 
% function [new_state, pdr_step] = Get_Next_Step_2D_UJI(old_state, MAP_X_LEN, MAP_Y_LEN)
% % Get_Next_Step_2D_UJI (V4 - 可变步长 + 地面真实角度)
% %
% %   此版本结合了 UJI 文件的两个部分：
% %   1. 传感器数据: 
% %      - (V3) 用于检测 *何时* 发生步伐 (findpeaks)。
% %      - (V4 新增) 用于估算每一步的 *可变步长* (Weinberg SLE)。
% %   2. 路径点数据: 
% %      - (V3) 用于计算每个走廊的 *完美角度*。
% %
% %   结果是一个 "完美" 的 PDR 步骤列表，由 (可变步长, 0) 
% %   和 (可变步长, 完美转角) 组成。
% %
% 
% % --- [ PDR 调整参数 ] ---
% persistent PDR_PARAMS
% if isempty(PDR_PARAMS)
%     PDR_PARAMS.MIN_PEAK_HEIGHT = 0.8; % 最小峰值高度 (m/s^2)
%     PDR_PARAMS.MIN_PEAK_DIST = 4;   % 最小峰值距离 (0.4s)
% 
%     % --- [V4 步长模型参数] ---
%     % PDR_PARAMS.FIXED_STEP_LEN = 0.7; % <-- V3 的做法 (已移除)
%     PDR_PARAMS.SLE_K = 0.45; % 步长估计 (SLE) 的 K 因子. 0.45 是一个合理的初始值
%     % --------------------------
% end
% % --------------------------
% 
% persistent pdr_steps_list step_counter
% 
% %==========================================================================
% % 1. 初始化模式 (Initialization Mode)
% %==========================================================================
% if ischar(old_state) && strcmp(old_state, 'init')
%     filename = MAP_X_LEN; 
%     fprintf('Initializing UJI PDR (V4 - Variable Step Length)... Loading file: %s\n', filename);
% 
%     % 加载并处理 UJI 文件 (V4)
%     [pdr_steps_list, num_detected] = load_and_process_uji_file_v4(filename, PDR_PARAMS);
% 
%     step_counter = 1; 
%     new_state = num_detected; % 这是返回给主脚本的步数
%     pdr_step = [];    
% 
%     fprintf('UJI PDR initialized. %d "perfect" variable steps synthesized.\n', new_state);
% 
%     if new_state == 0
%         warning('!!! V4: 检测到 0 步. 模拟循环将不会运行。');
%         warning('!!! 请尝试在 Get_Next_Step_2D_UJI.m 中调低 PDR_PARAMS。');
%     end
%     return;
% end
% 
% %==========================================================================
% % 2. 执行模式 (Execution Mode)
% %==========================================================================
% if isempty(pdr_steps_list)
%     error('UJI PDR not initialized. Call Get_Next_Step_2D_UJI(''init'', ''filename.txt'') before the simulation loop.');
% end
% 
% if step_counter > size(pdr_steps_list, 1)
%     pdr_step = [0, 0]; 
% else
%     pdr_step = pdr_steps_list(step_counter, :);
%     step_counter = step_counter + 1;
% end
% 
% %==========================================================================
% % 3. 应用PDR运动模型 (无变化)
% %==========================================================================
% step_len = pdr_step(1);
% d_theta = pdr_step(2);
% old_x = old_state(1);
% old_y = old_state(2);
% old_theta = old_state(3);
% new_theta = old_theta + d_theta;
% new_x = old_x + sin(new_theta) * step_len;
% new_y = old_y + cos(new_theta) * step_len;
% new_x = max(1, min(MAP_X_LEN, new_x));
% new_y = max(1, min(MAP_Y_LEN, new_y));
% new_state = [new_x, new_y, new_theta];
% end
% 
% 
% %==========================================================================
% % 内部辅助函数 (Local Helper Function) - (V4)
% %==========================================================================
% function [pdr_steps_list, N_total_steps] = load_and_process_uji_file_v4(filename, PDR_PARAMS)
% 
%     fid = fopen(filename);
%     if fid == -1, error('Cannot open file: %s', filename); end
% 
%     % --- 1. 读取 Part 1 (传感器数据) ---
%     sensor_lines = {};
%     while ~feof(fid)
%         line = fgetl(fid);
%         if isempty(line) || startsWith(strtrim(line), '<'), break; end
%         sensor_lines{end+1} = line;
%     end
%     sensor_data = cell2mat(cellfun(@(x) sscanf(x, '%f')', sensor_lines, 'UniformOutput', false)');
%     if isempty(sensor_data), error('No sensor data found in file: %s', filename); end
% 
%     % --- 2a. 检测步伐 (来自传感器) ---
%     accel_x = sensor_data(:, 5);
%     accel_y = sensor_data(:, 6);
%     accel_z = sensor_data(:, 7);
%     mag_accel = sqrt(accel_x.^2 + accel_y.^2 + accel_z.^2);
% 
%     % `pks` 是峰值 (acc_max), `step_indices` 是峰值的位置
%     [pks, step_indices] = findpeaks(mag_accel, 'MinPeakHeight', PDR_PARAMS.MIN_PEAK_HEIGHT, 'MinPeakDistance', PDR_PARAMS.MIN_PEAK_DIST);
% 
%     N_total_steps = length(step_indices);
%     if N_total_steps == 0
%         pdr_steps_list = zeros(0, 2); % 返回空矩阵
%         fclose(fid);
%         return;
%     end
% 
%     % --- 2b. [V4 新增] 估算可变步长 ---
%     step_lengths_list = zeros(N_total_steps, 1);
%     for i = 1:N_total_steps
%         acc_max = pks(i); % findpeaks 已经为我们找到了峰值
% 
%         % 寻找此步伐对应的"谷值" (acc_min)
%         if i == 1
%             search_start_idx = 1; % 对于第一步, 从头开始搜索
%         else
%             search_start_idx = step_indices(i-1); % 从上一个峰值开始搜索
%         end
%         search_end_idx = step_indices(i); % 搜索到当前峰值
% 
%         acc_min = min(mag_accel(search_start_idx : search_end_idx));
% 
%         % Weinberg 模型
%         step_len = PDR_PARAMS.SLE_K * (acc_max - acc_min)^(1/4);
% 
%         step_lengths_list(i) = step_len;
%     end
% 
%     % --- 3. 读取 Part 2 (路径点数据) ---
%     if feof(fid)
%         % [回退逻辑] 未找到路径点数据，使用 V2 (嘈杂的) 角度
%         warning('未找到路径点数据 (Part 2)，将退回到 V2（嘈杂的）PDR 角度。');
%         fclose(fid);
% 
%         pdr_steps_list = zeros(N_total_steps, 2);
%         azimuth_rad = deg2rad(sensor_data(:, 8));
%         last_step_azimuth = azimuth_rad(step_indices(1));
% 
%         for i = 1:N_total_steps
%             current_step_azimuth = azimuth_rad(step_indices(i));
%             d_theta = current_step_azimuth - last_step_azimuth;
% 
%             pdr_steps_list(i, 1) = step_lengths_list(i); % <-- [V4] 使用可变步长
%             pdr_steps_list(i, 2) = atan2(sin(d_theta), cos(d_theta));
% 
%             last_step_azimuth = current_step_azimuth;
%         end
%         return;
%     end
% 
%     % [V3/V4 正常逻辑] 继续读取路径点
%     num_segments_str = strtrim(line);
%     m = sscanf(num_segments_str, '<%d>');
%     segment_data = zeros(m, 6);
%     for i = 1:m
%         line = fgetl(fid);
%         segment_data(i, :) = sscanf(line, '%f')';
%     end
%     fclose(fid);
% 
%     % --- 4. 计算每个走廊的"完美"角度 ---
%     segment_headings = zeros(m, 1);
%     for i = 1:m
%         segment_headings(i) = calculate_bearing(segment_data(i, 1), segment_data(i, 2), segment_data(i, 3), segment_data(i, 4));
%     end
% 
%     % --- 5. 合成 "完美" PDR 步骤列表 ---
%     pdr_steps_list = zeros(N_total_steps, 2);
%     last_heading = segment_headings(1);
%     last_segment_idx = 1;
% 
%     for i = 1:N_total_steps % 遍历检测到的每一步
% 
%         current_step_sample_index = step_indices(i);
% 
%         % 找到这一步属于哪一个走廊 (segment)
%         current_segment_idx = find(current_step_sample_index >= (segment_data(:, 5)+1) & ...
%                                    current_step_sample_index <= (segment_data(:, 6)+1), 1, 'first');
% 
%         if isempty(current_segment_idx)
%             current_segment_idx = last_segment_idx; % 容错：假设还在上一个走廊
%         end
% 
%         current_heading = segment_headings(current_segment_idx);
%         current_step_length = step_lengths_list(i); % <-- [V4] 使用可变步长
% 
%         % 检查这是否是进入新走廊的第一步
%         if current_segment_idx ~= last_segment_idx
%             % 是的, 计算完美转弯
%             d_theta = current_heading - last_heading;
%             d_theta = atan2(sin(d_theta), cos(d_theta)); % 环绕 -pi 到 +pi
%         else
%             % 不是, 这是一个直线步伐
%             d_theta = 0;
%         end
% 
%         % 存储最终的 [步长, 转向角]
%         pdr_steps_list(i, 1) = current_step_length;
%         pdr_steps_list(i, 2) = d_theta;
% 
%         % 为下一次循环更新状态
%         last_heading = current_heading;
%         last_segment_idx = current_segment_idx;
%     end
% end
% 
% 
% function bearing_rad = calculate_bearing(lat1, lon1, lat2, lon2)
% % 计算从 (lat1, lon1) 到 (lat2, lon2) 的方位角 (bearing)
% % 角度以弧度为单位, 0 = 北, +pi/2 = 东
%     lat1r = deg2rad(lat1); lon1r = deg2rad(lon1);
%     lat2r = deg2rad(lat2); lon2r = deg2rad(lon2);
%     dLon = lon2r - lon1r;
%     y = sin(dLon) * cos(lat2r);
%     x = cos(lat1r) * sin(lat2r) - sin(lat1r) * cos(lat2r) * cos(dLon);
%     bearing_rad = atan2(y, x);
% end
function [new_state, pdr_step] = Get_Next_Step_2D_UJI(old_state, MAP_X_LEN, MAP_Y_LEN)
% Get_Next_Step_2D_UJI (V5 - 真实测量值)
%
%   此版本在 'init' 模式下，会返回一个包含 *所有* 必需数据的结构体：
%   1. UJI_Data.pdr_steps: (Nx2) "完美" PDR 步骤 [步长, 转向角]
%   2. UJI_Data.mag_readings: (Nx3) *真实* 磁场读数 [mx, my, mz]
%   3. UJI_Data.num_steps: (N) 检测到的总步数
%
%   执行模式 ('Execution Mode') 保持不变，仅用于生成 'true_state' 路径
%

% --- [ PDR 调整参数 ] ---
persistent PDR_PARAMS
if isempty(PDR_PARAMS)
    PDR_PARAMS.MIN_PEAK_HEIGHT = 0.8; 
    PDR_PARAMS.MIN_PEAK_DIST = 4;   
    PDR_PARAMS.SLE_K = 0.45; 
end
% --------------------------

% --- [ V5 持久化变量 ] ---
% pdr_steps_list 现在只用于 "true_state" 的模拟
persistent pdr_steps_list step_counter UJI_Data_Cache

%==========================================================================
% 1. 初始化模式 (Initialization Mode)
%==========================================================================
if ischar(old_state) && strcmp(old_state, 'init')
    filename_str = char(MAP_X_LEN); % 转换 string 为 char
    fprintf('Initializing UJI PDR (V5 - Real Measurements)... Loading file: %s\n', filename_str);
    
    % 加载并处理 UJI 文件 (V5)
    [pdr_steps_list_local, N_total_steps, mag_readings_at_steps] = load_and_process_uji_file_v5(filename_str, PDR_PARAMS);
    
    % --- 存储 V5 数据 ---
    % 存储用于执行模式 (模拟 true_path)
    pdr_steps_list = pdr_steps_list_local;
    step_counter = 1; 
    
    % 存储用于返回给主脚本 (用于 'live_sequence')
    UJI_Data_Cache.pdr_steps = pdr_steps_list_local;
    UJI_Data_Cache.mag_readings = mag_readings_at_steps;
    UJI_Data_Cache.num_steps = N_total_steps;
    
    % 返回这个包含所有数据的结构体
    new_state = UJI_Data_Cache; 
    pdr_step = [];    
    
    fprintf('UJI PDR initialized. %d variable steps and real mag readings synthesized.\n', N_total_steps);
    
    if N_total_steps == 0
        warning('!!! V5: 检测到 0 步. 模拟循环将不会运行。');
    end
    return;
end

%==========================================================================
% 2. 执行模式 (Execution Mode - 用于生成 "True Path")
%==========================================================================
% 此部分保持不变。它的 *唯一* 目的现在是为绘图生成 "true_path_history"
% 它从 'init' 期间创建的 pdr_steps_list 中提取数据。
%
if isempty(pdr_steps_list)
    error('UJI PDR not initialized. Call Get_Next_Step_2D_UJI(''init'', ''filename.txt'') before the simulation loop.');
end

if step_counter > size(pdr_steps_list, 1)
    pdr_step = [0, 0]; 
else
    pdr_step = pdr_steps_list(step_counter, :);
    step_counter = step_counter + 1;
end

%==========================================================================
% 3. 应用PDR运动模型 (无变化)
%==========================================================================
step_len = pdr_step(1);
d_theta = pdr_step(2);
old_x = old_state(1);
old_y = old_state(2);
old_theta = old_state(3);
new_theta = old_theta + d_theta;
new_x = old_x + sin(new_theta) * step_len;
new_y = old_y + cos(new_theta) * step_len;
new_x = max(1, min(MAP_X_LEN, new_x));
new_y = max(1, min(MAP_Y_LEN, new_y));
new_state = [new_x, new_y, new_theta];
end


%==========================================================================
% 内部辅助函数 (Local Helper Function) - (V5)
%==========================================================================
function [pdr_steps_list, N_total_steps, mag_readings_at_steps] = load_and_process_uji_file_v5(filename, PDR_PARAMS)
    
    fid = fopen(filename);
    if fid == -1, error('Cannot open file: %s', filename); end
    
    % --- 1. 读取 Part 1 (传感器数据) ---
    sensor_lines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(line) || startsWith(strtrim(line), '<'), break; end
        sensor_lines{end+1} = line;
    end
    sensor_data = cell2mat(cellfun(@(x) sscanf(x, '%f')', sensor_lines, 'UniformOutput', false)');
    if isempty(sensor_data), error('No sensor data found in file: %s', filename); end
    
    % --- 2a. 检测步伐 (来自传感器) ---
    accel_x = sensor_data(:, 5);
    accel_y = sensor_data(:, 6);
    accel_z = sensor_data(:, 7);
    mag_accel = sqrt(accel_x.^2 + accel_y.^2 + accel_z.^2);
    
    [pks, step_indices] = findpeaks(mag_accel, 'MinPeakHeight', PDR_PARAMS.MIN_PEAK_HEIGHT, 'MinPeakDistance', PDR_PARAMS.MIN_PEAK_DIST);
    
    N_total_steps = length(step_indices);
    if N_total_steps == 0
        pdr_steps_list = zeros(0, 2);
        mag_readings_at_steps = zeros(0, 3);
        fclose(fid);
        return;
    end

    % --- 2b. [V4] 估算可变步长 ---
    step_lengths_list = zeros(N_total_steps, 1);
    for i = 1:N_total_steps
        acc_max = pks(i); 
        if i == 1, search_start_idx = 1; else, search_start_idx = step_indices(i-1); end
        search_end_idx = step_indices(i);
        acc_min = min(mag_accel(search_start_idx : search_end_idx));
        step_len = PDR_PARAMS.SLE_K * (acc_max - acc_min)^(1/4);
        step_lengths_list(i) = step_len;
    end
    
    % --- 2c. [V5 新增] 提取真实的磁场读数 ---
    % sensor_data 列: 2=mx, 3=my, 4=mz
    mag_readings_at_steps = sensor_data(step_indices, [2, 3, 4]);

    % --- 3. 读取 Part 2 (路径点数据) ---
    if feof(fid)
        % [回退逻辑]
        warning('未找到路径点数据 (Part 2)，将退回到 V2（嘈杂的）PDR 角度。');
        fclose(fid);
        
        pdr_steps_list = zeros(N_total_steps, 2);
        azimuth_rad = deg2rad(sensor_data(:, 8));
        last_step_azimuth = azimuth_rad(step_indices(1));
        
        for i = 1:N_total_steps
            current_step_azimuth = azimuth_rad(step_indices(i));
            d_theta = current_step_azimuth - last_step_azimuth;
            pdr_steps_list(i, 1) = step_lengths_list(i); % [V4] 步长
            pdr_steps_list(i, 2) = atan2(sin(d_theta), cos(d_theta));
            last_step_azimuth = current_step_azimuth;
        end
        return; % 在这里返回，使用嘈杂的角度
    end
    
    % [V3/V4/V5 正常逻辑] 继续读取路径点
    num_segments_str = strtrim(line);
    m = sscanf(num_segments_str, '<%d>');
    segment_data = zeros(m, 6);
    for i = 1:m
        line = fgetl(fid);
        segment_data(i, :) = sscanf(line, '%f')';
    end
    fclose(fid);

    % --- 4. 计算每个走廊的"完美"角度 ---
    segment_headings = zeros(m, 1);
    for i = 1:m
        segment_headings(i) = calculate_bearing(segment_data(i, 1), segment_data(i, 2), segment_data(i, 3), segment_data(i, 4));
    end

    % --- 5. 合成 "完美" PDR 步骤列表 ---
    pdr_steps_list = zeros(N_total_steps, 2);
    last_heading = segment_headings(1);
    last_segment_idx = 1;

    for i = 1:N_total_steps % 遍历检测到的每一步
        current_step_sample_index = step_indices(i);
        current_segment_idx = find(current_step_sample_index >= (segment_data(:, 5)+1) & ...
                                   current_step_sample_index <= (segment_data(:, 6)+1), 1, 'first');
        if isempty(current_segment_idx)
            current_segment_idx = last_segment_idx; 
        end
        
        current_heading = segment_headings(current_segment_idx);
        current_step_length = step_lengths_list(i); % [V4]
        
        if current_segment_idx ~= last_segment_idx
            d_theta = current_heading - last_heading;
            d_theta = atan2(sin(d_theta), cos(d_theta)); 
        else
            d_theta = 0;
        end
        
        pdr_steps_list(i, 1) = current_step_length;
        pdr_steps_list(i, 2) = d_theta;
        
        last_heading = current_heading;
        last_segment_idx = current_segment_idx;
    end
end


function bearing_rad = calculate_bearing(lat1, lon1, lat2, lon2)
    lat1r = deg2rad(lat1); lon1r = deg2rad(lon1);
    lat2r = deg2rad(lat2); lon2r = deg2rad(lon2);
    dLon = lon2r - lon1r;
    y = sin(dLon) * cos(lat2r);
    x = cos(lat1r) * sin(lat2r) - sin(lat1r) * cos(lat2r) * cos(dLon);
    bearing_rad = atan2(y, x);
end