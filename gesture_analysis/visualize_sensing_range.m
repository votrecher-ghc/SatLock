% =========================================================================
% visualize_sensing_range_v3.m
% 功能: 可视化手势平面的个体感知范围与整体整合包络 (美化版 + 底部大号标号)
% 风格: 现代科研风格，柔和配色，底部文本框字号加大
% 核心逻辑:
%   1. 寻找卫星最多的最佳时刻。
%   2. [个体范围]: 以每个投影点为中心，画一个物理半径的圆。
%   3. [整合范围]: 计算所有个体圆周点的凸包 (Convex Hull)。
%   4. [底部标号]: 在坐标轴下方添加文本框显示感知平面的高度。
% =========================================================================
%% 1. 环境检查与参数设置
clearvars -except obs_data nav_data;
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end
addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));
fprintf('--> 启动感知范围可视化 (Aesthetic Mode + Large Label)...\n');
% --- 投影与感知参数 ---
TRAJ.gesture_height    = 0.4;  % 手势平面高度 (米)
TRAJ.min_elevation     = 25;    % 最低仰角过滤 (度)
TRAJ.sensing_radius    = 0.1;  % 个体感知半径 (米)
%% 2. 寻找最佳观测时刻 (快照)
fprintf('--> [计算] 扫描卫星数量最多的最佳时刻...\n');
num_epochs = length(obs_data);
best_epoch_idx = -1; max_sat_count = -1;
all_sat_ids = {}; for i=1:min(100,length(obs_data)), if~isempty(obs_data(i).data), all_sat_ids=[all_sat_ids,fieldnames(obs_data(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids); valid_sats_list = {};
for i=1:length(unique_sat_ids), sid=unique_sat_ids{i}; if ismember(sid(1),['G','C','E','J']), valid_sats_list{end+1}=sid; end; end
for t_idx = 1:num_epochs 
    try [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, t_idx); catch, continue; end
    if isempty(rec_pos) || all(isnan(rec_pos)), continue; end
    [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
    current_count = 0;
    for s = 1:length(valid_sats_list)
        sid = valid_sats_list{s}; if ~isfield(sat_states, sid), continue; end
        sat_p = sat_states.(sid).position;
        [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
        vec_u = [e, n, u]/norm([e, n, u]); zen_deg = acosd(vec_u(3));
        if vec_u(3)>0 && (90-zen_deg)>=TRAJ.min_elevation, current_count = current_count + 1; end
    end
    if current_count > max_sat_count, max_sat_count = current_count; best_epoch_idx = t_idx; end
end
if best_epoch_idx == -1, error('未找到有效时刻'); end
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
    sid = valid_sats_list{s}; if ~isfield(sat_states, sid), continue; end
    sat_p = sat_states.(sid).position;
    [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
    vec_u = [e, n, u]/norm([e, n, u]); zen_deg = acosd(vec_u(3));
    
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
f = figure('Name', 'Sensing Range Analysis v3.1', 'Position', [400, 150, 750, 750], 'Color', 'w');
ax = axes('Parent', f); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
xlabel('East (m)'); ylabel('North (m)');
% --- 配色方案 (Color Palette) ---
col_recv = 'k';                     % 接收机: 纯黑
col_sat  = [0.9290, 0.4940, 0.1250]; % 卫星点: 活力橙 (Matlab Default Orange)
col_ind  = [0.6, 0.6, 0.6];         % 个体范围: 中性灰
col_hull_line = [0.0, 0.5, 0.8];    % 整体包络线: 柔和亮蓝 (Cerulean Blue) - 不深
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
    '\bf\fontsize{12}Sensing Scope Analysis',
    sprintf('\\rm\\fontsize{10}Height: %.2fm', TRAJ.gesture_height),
    sprintf('\\rm\\fontsize{10}Visible Satellites: %d | Total Area: %.2f m^2', size(proj_centers, 1), area_val)
};
title(ax, title_str);
legend(ax, 'Location', 'best'); % 自动寻找最佳位置，避免遮挡
% 视野微调
if ~isempty(hull_x)
    max_range = max(max(abs(hull_x)), max(abs(hull_y))) * 1.15;
    xlim(ax, [-max_range, max_range]); ylim(ax, [-max_range, max_range]);
end

% --- 新增: 在图像底部添加高度标号 ---
% 使用 annotation 创建一个不依赖于坐标轴的文本框
% position 坐标是归一化的 [x_begin, y_begin, width, height] (相对于整个figure窗口)
label_str = sprintf('\\bf (c) Sensing Plane Height: %.2f meters above receiver', TRAJ.gesture_height);
annotation(f, 'textbox',...
    [0.1, 0.01, 0.8, 0.06],... % <--- 修改：宽度从0.6加大到0.8，高度从0.05加大到0.06，以容纳大字体
    'String', label_str,...
    'FontSize', 15,...         % <--- 修改：字号从 11 加大到 15
    'FontWeight', 'bold',...
    'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'middle',...
    'LineStyle', 'none',...    % 无边框
    'BackgroundColor', 'w');   % 白色背景

fprintf('✅ 可视化完成。总包络面积: %.2f m²，已在底部添加大号高度标号。\n', area_val);