%==========================================================================
% SCRIPT: Build_UJI_Hybrid_Map.m (V2 - 你的高斯噪声方案)
%
% 描述:
%   "建图"脚本 V2。
%   此版本实现了你的新逻辑:
%   1. 正常插值 "岛屿" (走廊) 数据。
%   2. 找到 "海洋" (NaN) 区域。
%   3. 不再用一个定值填充, 而是用高斯噪声 (randn) 填充 "海洋" 区域。
%
%   最后, 它会将构建好的地图保存为 'my_hybrid_map.mat' 文件。
%
% 用法:
%   1. 将此脚本与 'lines' 和 'curves' 文件夹放在同一目录。
%   2. 运行此脚本。
%
%==========================================================================
clc;
clear;
close all;

%% 1. 设置参数
fprintf('正在自动检测脚本路径...\n');
try
    script_path = mfilename('fullpath');
    if isempty(script_path)
        script_path = pwd; 
        script_dir = pwd;
    else
        [script_dir, ~, ~] = fileparts(script_path);
    end
catch ME
    warning('获取脚本路径失败, 将使用当前工作目录。');
    script_dir = pwd; % Fallback
end

TRAIN_DATA_DIR = '.\UJIIndoorLoc-Mag\UJIIndoorLoc-Mag\'; 
OUTPUT_MAP_FILENAME = 'my_hybrid_map.mat'; % [新] 保存为新文件
GRID_RESOLUTION = 0.5; 
SMOOTHING_FACTOR = 2.0; 

fprintf('地图将保存为: %s\n', fullfile(script_dir, OUTPUT_MAP_FILENAME));

%% 2. 递归查找所有训练文件
fprintf('正在递归查找所有训练文件...\n');
line_files = dir(fullfile(TRAIN_DATA_DIR, 'lines', '**', 'l*.txt'));
curve_files = dir(fullfile(TRAIN_DATA_DIR, 'curves', '**', 'c*.txt'));
all_files = [line_files; curve_files];
if isempty(all_files)
    error('在 "%s" 的 "lines" 或 "curves" 子文件夹中未找到任何 .txt 文件。', TRAIN_DATA_DIR);
end
fprintf('... 找到了 %d 个训练文件。\n', length(all_files));

%% 3. 处理文件并构建 (X, Y, Mag) 点云
fprintf('正在处理所有文件以构建磁场点云...\n');
all_lats = [];
all_lons = [];
all_mags = [];
h_waitbar = waitbar(0, '正在处理训练文件...');

for i = 1:length(all_files)
    waitbar(i/length(all_files), h_waitbar, sprintf('处理中 (%d/%d)', i, length(all_files)));
    filepath = fullfile(all_files(i).folder, all_files(i).name);
    try
        [sensor_data, gt_segments] = parse_UJI_file(filepath);
    catch ME
        warning('跳过文件 (解析失败): %s. 错误: %s', all_files(i).name, ME.message);
        continue;
    end
    if isempty(gt_segments)
        continue;
    end
    
    for j = 1:size(gt_segments, 1)
        seg = gt_segments(j, :);
        lat1 = seg(1); lon1 = seg(2); lat2 = seg(3); lon2 = seg(4);
        fs = seg(5) + 1; ls = seg(6) + 1;
        num_samples = ls - fs + 1;
        
        if num_samples <= 0 || fs > height(sensor_data) || ls > height(sensor_data)
            continue;
        end
        
        seg_lats = linspace(lat1, lat2, num_samples)';
        seg_lons = linspace(lon1, lon2, num_samples)';
        seg_data = sensor_data(fs:ls, :);
        seg_mags = sqrt(seg_data.mx.^2 + seg_data.my.^2 + seg_data.mz.^2);
        
        all_lats = [all_lats; seg_lats];
        all_lons = [all_lons; seg_lons];
        all_mags = [all_mags; seg_mags];
    end
end
close(h_waitbar);
fprintf('... 点云构建完毕。总共 %d 个指纹点。\n', length(all_lats));
if isempty(all_lats)
    error('未能从任何训练文件中提取有效的指纹点。');
end

%% 4. 插值并创建 *混合* 地图矩阵
fprintf('正在将点云插值为 %s 的网格地图...\n', OUTPUT_MAP_FILENAME);

origin_latlon = [min(all_lats), min(all_lons)];
[map_x, map_y] = latlon_to_xy_helper(all_lats, all_lons, origin_latlon(1), origin_latlon(2));

x_min = min(map_x); x_max = max(map_x);
y_min = min(map_y); y_max = max(map_y);
[X_grid, Y_grid] = meshgrid(x_min:GRID_RESOLUTION:x_max, ...
                            y_min:GRID_RESOLUTION:y_max);
                        
fprintf('... 正在插值 "岛屿" (网格尺寸: %d x %d)...\n', size(X_grid, 1), size(X_grid, 2));

interpolant = scatteredInterpolant(map_x, map_y, all_mags, 'natural', 'none');
Mag_map = interpolant(X_grid, Y_grid); % [!!] 这张图现在包含了 NaN

% --- [!! 你的新逻辑从这里开始 !!] ---
fprintf('... 正在填充 "海洋" 区域...\n');

nan_mask = isnan(Mag_map);
num_nans = sum(nan_mask, 'all');

if num_nans > 0
    fprintf('... 找到了 %d 个 "海洋" (NaN) 点。\n', num_nans);
    
    % 1. 计算 "岛屿" (真实数据) 的统计特征
    island_data = Mag_map(~nan_mask);
    data_mean = mean(island_data);
    data_std = 5;%std(island_data); % 使用真实数据的标准差
    
    fprintf('... "岛屿" 数据的 均值=%.2f, 标准差=%.2f\n', data_mean, data_std);
    
    % 2. [你的要求] 用高斯噪声填充 "海洋"
    %    我们用 "岛屿" 的平均值和方差来生成噪声, 这样最真实
    fprintf('... 正在生成 %d 个高斯噪声点 (均值=%.2f, 标准差=%.2f)\n', ...
            num_nans, data_mean, data_std);
            
    noise = data_mean + randn(num_nans, 1) * data_std;
    
    % 3. 填充
    Mag_map(nan_mask) = noise;
else
    fprintf('... 地图已插值, 没有 "海洋" (NaN) 点。\n');
end
% --- [!! 你的新逻辑结束 !!] ---


% 4f. 平滑 *整个* 地图 (岛屿和海洋)
Mag_map = imgaussfilt(Mag_map, SMOOTHING_FACTOR); 
fprintf('... 整个混合地图已平滑。\n');


%% 5. 保存最终的地图文件
fprintf('正在封装并保存混合地图到 %s ...\n', OUTPUT_MAP_FILENAME);

geo_map_cpu.X_grid = X_grid;
geo_map_cpu.Y_grid = Y_grid;
geo_map_cpu.Mag_map = Mag_map; % 这就是你的 "二维大矩阵" (混合版)

save(OUTPUT_MAP_FILENAME, 'geo_map_cpu', 'origin_latlon', '-v7.3');

fprintf('\n=== 混合建图完成! ===\n');
fprintf('地图已保存到: %s\n', fullfile(script_dir, OUTPUT_MAP_FILENAME));

%% 6. (可选) 可视化检查
fprintf('正在生成地图预览图...\n');
figure('Name', '混合地磁图预览 (Build_UJI_Hybrid_Map.m)', 'Position', [100, 100, 800, 600]);
imagesc(geo_map_cpu.X_grid(1,:), geo_map_cpu.Y_grid(:,1), geo_map_cpu.Mag_map);
hold on;
axis xy; colormap('jet'); colorbar;
title(sprintf('构建的 *混合* 地磁图 (%d x %d 网格)', size(X_grid, 1), size(X_grid, 2)));
xlabel('X 位置 (米)'); ylabel('Y 位置 (米)');
hold off;
fprintf('预览图已生成。\n');


%% ========================================================================
%  == 辅助函数定义区 (与之前相同) ==
%  ========================================================================

function [sensor_data, gt_segments] = parse_UJI_file(filename)
    try
        raw_lines = readlines(filename);
    catch ME
        error('无法读取文件: %s. 错误: %s', filename, ME.message);
    end
    separator_idx = find(startsWith(raw_lines, '<'), 1);
    if isempty(separator_idx)
        error('未在文件中找到分隔符 "<m>": %s', filename);
    end
    sensor_lines = raw_lines(1 : separator_idx-1);
    sensor_lines = sensor_lines(strlength(sensor_lines) > 0);
    if isempty(sensor_lines)
         error('在文件中未找到传感器数据: %s', filename);
    end
    C_sensor = textscan(strjoin(sensor_lines, '\n'), '%f %f %f %f %f %f %f %f %f %f');
    if isempty(C_sensor) || isempty(C_sensor{1})
        error('无法解析文件中的传感器数据: %s', filename);
    end
    sensor_matrix = cell2mat(C_sensor);
    sensor_data = array2table(sensor_matrix, 'VariableNames', ...
        {'ts', 'mx', 'my', 'mz', 'ax', 'ay', 'az', 'ox', 'oy', 'oz'});
    gt_lines = raw_lines(separator_idx+1 : end);
    gt_lines = gt_lines(strlength(gt_lines) > 0); 
    if isempty(gt_lines)
        gt_segments = [];
    else
        C_gt = textscan(strjoin(gt_lines, '\n'), '%f %f %f %f %f %f');
        if isempty(C_gt) || isempty(C_gt{1})
             error('无法解析文件中的GT段数据: %s', filename);
        end
        gt_segments = cell2mat(C_gt);
    end
end
function [x, y] = latlon_to_xy_helper(lat, lon, origin_lat, origin_lon)
    R = 6371000; % 地球半径 (米)
    lat_rad = deg2rad(lat);
    lon_rad = deg2rad(lon);
    origin_lat_rad = deg2rad(origin_lat);
    origin_lon_rad = deg2rad(origin_lon);
    x = R * (lon_rad - origin_lon_rad) .* cos(origin_lat_rad);
    y = R * (lat_rad - origin_lat_rad);
end