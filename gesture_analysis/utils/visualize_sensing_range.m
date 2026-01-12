% =========================================================================
% visualize_sensing_range.m (函数版)
% 功能: 可视化手势平面的个体感知范围与整体整合包络 (v3.1)
% 描述:
%   该函数计算并展示在指定高度 (如 0.4m) 的手势平面上，各颗卫星信号的
%   投影点及其有效的感知范围 (Sensing Radius)。
%   同时计算这些感知区域的整体凸包 (Convex Hull)，量化总感知面积。
%
% [调用格式]:
%   visualize_sensing_range(obs_data, nav_data);
%
% [输入参数]:
%   1. obs_data (struct数组): 原始观测数据。
%   2. nav_data (struct结构体): 导航星历数据。
%
% [核心逻辑]:
%   1. 最佳时刻扫描: 遍历所有历元，自动锁定可见卫星数量最多的时刻作为快照。
%   2. 几何投影: 将视线矢量投影到指定高度平面。
%   3. 范围绘制:
%      - 个体范围: 以投影点为中心，绘制物理半径 (0.1m) 的圆。
%      - 整体范围: 计算所有圆周点的凸包 (Convex Hull) 并填充。
% =========================================================================

function visualize_sensing_range(obs_data, nav_data)

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动感知范围可视化 (Function版 v3.2: All Satellites)...\n');

%% 1. 参数设置
% --- 投影与感知参数 ---
TRAJ.gesture_height    = 0.3;   % [物理] 手势平面高度 (米)
TRAJ.min_elevation     = 0;     % [物理] 最低仰角过滤 (度)
TRAJ.sensing_radius    = 0.5;   % [物理] 个体感知半径 (米)

%% 2. 寻找最佳观测时刻 (快照)
fprintf('--> [计算] 扫描卫星数量最多的最佳时刻...\n');

num_epochs = length(obs_data);
best_epoch_idx = -1; 
max_sat_count = -1;

% 提取有效卫星列表
all_sat_ids = {}; 
for i=1:min(100,length(obs_data))
    if ~isempty(obs_data(i).data)
        all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; 
    end
end
unique_sat_ids = unique(all_sat_ids); 

% 【关键修改】不再进行 G/C/E/J 过滤，直接使用所有扫描到的卫星
valid_sats_list = unique_sat_ids; 

% 遍历寻找最佳时刻
for t_idx = 1:num_epochs 
    try 
        [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, t_idx); 
    catch
        continue; 
    end
    
    if isempty(rec_pos) || all(isnan(rec_pos)), continue; end
    [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
    
    current_count = 0;
    for s = 1:length(valid_sats_list)
        sid = valid_sats_list{s}; 
        if ~isfield(sat_states, sid), continue; end
        
        sat_p = sat_states.(sid).position;
        [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
        vec_u = [e, n, u]/norm([e, n, u]); 
        zen_deg = acosd(vec_u(3));
        
        if vec_u(3)>0 && (90-zen_deg)>=TRAJ.min_elevation
            current_count = current_count + 1; 
        end
    end
    
    if current_count > max_sat_count
        max_sat_count = current_count; 
        best_epoch_idx = t_idx; 
    end
end

if best_epoch_idx == -1
    error('未找到有效时刻，无法进行可视化。');
end

fprintf('✅ 锁定最佳时刻: %s (可见卫星: %d)\n', datestr(obs_data(best_epoch_idx).time), max_sat_count);

%% 3. 计算个体范围和整合包络
fprintf('--> [计算] 生成个体圆形范围和整体凸包...\n');

[rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, best_epoch_idx);
[lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));

proj_centers = []; 
all_circle_points = []; 

% 辅助角度
theta = linspace(0, 2*pi, 64); % 增加点数让圆更圆滑
circle_x_base = TRAJ.sensing_radius * cos(theta);
circle_y_base = TRAJ.sensing_radius * sin(theta);

for s = 1:length(valid_sats_list)
    sid = valid_sats_list{s}; 
    if ~isfield(sat_states, sid), continue; end
    
    sat_p = sat_states.(sid).position;
    [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
    vec_u = [e, n, u]/norm([e, n, u]); 
    zen_deg = acosd(vec_u(3));
    
    if vec_u(3) <= 0 || (90-zen_deg) < TRAJ.min_elevation, continue; end
    
    t_int = TRAJ.gesture_height / vec_u(3);
    pt_int = t_int * vec_u;
    
    if norm(pt_int(1:2)) > 5.0, continue; end 
    
    center_e = pt_int(1); center_n = pt_int(2);
    proj_centers = [proj_centers; center_e, center_n];
    
    this_circle_x = center_e + circle_x_base;
    this_circle_y = center_n + circle_y_base;
    all_circle_points = [all_circle_points; this_circle_x', this_circle_y'];
end

if ~isempty(all_circle_points)
    k_hull = convhull(all_circle_points(:,1), all_circle_points(:,2));
    hull_x = all_circle_points(k_hull, 1);
    hull_y = all_circle_points(k_hull, 2);
    area_val = polyarea(hull_x, hull_y);
else
    hull_x = []; hull_y = []; area_val = 0;
end

%% 4. 绘图 (Custom Aesthetics + Bottom Label)
fprintf('--> [绘图] 生成最终可视化...\n');

f = figure('Name', 'Sensing Range Analysis v3.2', 'Position', [400, 150, 750, 750], 'Color', 'w');
ax = axes('Parent', f); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
xlabel('East (m)'); ylabel('North (m)');

% --- 配色方案 ---
col_recv = 'k';                     % 接收机: 纯黑
col_sat  = [0.9290, 0.4940, 0.1250]; % 卫星点: 活力橙
col_ind  = [0.6, 0.6, 0.6];         % 个体范围: 中性灰
col_hull_line = [0.0, 0.5, 0.8];    % 整体包络线: 柔和亮蓝
col_hull_fill = [0.0, 0.6, 0.9];    % 整体填充: 天蓝色

% A. 画接收机
plot(ax, 0, 0, '^', 'MarkerSize', 11, 'MarkerFaceColor', col_recv, 'MarkerEdgeColor', 'none', 'DisplayName', 'Receiver');

% B. 画个体感知范围
if ~isempty(proj_centers)
    % 1. 画个体圆圈 (虚线，灰色)
    for i = 1:size(proj_centers, 1)
        cx = proj_centers(i,1) + circle_x_base;
        cy = proj_centers(i,2) + circle_y_base;
        if i==1
             plot(ax, cx, cy, '--', 'Color', col_ind, 'LineWidth', 1, 'DisplayName', sprintf('Individual Range (R=%.2fm)', TRAJ.sensing_radius));
        else
             plot(ax, cx, cy, '--', 'Color', col_ind, 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
    
    % 2. 画卫星中心点 (橙色实心点，最后画以防遮挡)
    plot(ax, proj_centers(:,1), proj_centers(:,2), '.', 'Color', col_sat, 'MarkerSize', 18, 'DisplayName', 'Projected Satellites');
end

% C. 画整合感知范围 (最外层)
if ~isempty(hull_x)
    % 填充 (极淡的蓝色，增加范围感)
    fill(ax, hull_x, hull_y, col_hull_fill, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    % 边界线 (柔和亮蓝，加粗)
    plot(ax, hull_x, hull_y, '-', 'Color', col_hull_line, 'LineWidth', 2.5, 'DisplayName', 'Total Sensing Scope');
end

% 标题与图例优化
title_str = {
    '\bf\fontsize{12}Sensing Scope Analysis (All Systems)',
    sprintf('\\rm\\fontsize{10}Height: %.2fm', TRAJ.gesture_height),
    sprintf('\\rm\\fontsize{10}Visible Satellites: %d | Total Area: %.2f m^2', size(proj_centers, 1), area_val)
};
title(ax, title_str);
legend(ax, 'Location', 'best');

% 视野微调
if ~isempty(hull_x)
    max_range = max(max(abs(hull_x)), max(abs(hull_y))) * 1.15;
    xlim(ax, [-max_range, max_range]); ylim(ax, [-max_range, max_range]);
end

% --- 在图像底部添加高度标号 ---
% 使用 annotation 创建一个不依赖于坐标轴的文本框
label_str = sprintf('\\bf (c) Sensing Plane Height: %.2f meters above receiver', TRAJ.gesture_height);
annotation(f, 'textbox',...
    [0.1, 0.01, 0.8, 0.06],... 
    'String', label_str,...
    'FontSize', 15,...         
    'FontWeight', 'bold',...
    'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'middle',...
    'LineStyle', 'none',...    
    'BackgroundColor', 'w');   

fprintf('✅ 可视化完成。总包络面积: %.2f m²\n', area_val);
end