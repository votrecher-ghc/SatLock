% =========================================================================
% step2_direction_estimation.m (函数版)
% 功能: 手势分析第二步 - 3D 轨迹方向估算 (v9 3D Rigorous)
%
% [调用格式]:
%   step2_direction_estimation(obs_data, nav_data, segments, volatility_matrix, t_grid, valid_sats);
%
% [输入参数]:
%   1. obs_data, nav_data: 基础数据。
%   2. segments: Step 1 返回的分段结构体。
%   3. volatility_matrix: Step 1 返回的波动矩阵。
%   4. t_grid: Step 1 返回的时间轴。
%   5. valid_sats: Step 1 返回的卫星列表。
%
% [核心算法]:
%   - 3D 视距切割模型 (LoS Cutting Model)。
%   - 建立局部 ENU 坐标系。
%   - 计算射线与 Z=h 平面的 3D 交点。
% =========================================================================

function step2_direction_estimation(obs_data, nav_data, segments, volatility_matrix, t_grid, valid_sats)

if isempty(segments)
    fprintf('⚠️ 警告: 未收到有效的手势分段，Step 2 终止。\n');
    return;
end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> [Step 2] 启动 3D 空间轨迹分析...\n');

%% 1. 初始化 3D 画布
fig_handle = figure('Name', '3D Gesture Trajectory Model', 'Position', [100, 100, 1200, 900], 'Color', 'w');
ax = axes('Parent', fig_handle);
hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal'); view(ax, 3);
xlabel(ax, 'East (m)'); ylabel(ax, 'North (m)'); zlabel(ax, 'Up (m)');
title(ax, '基于三维视距路径切割的手势轨迹推演');

% --- 参数设置 ---
TRAJ_PARA.gesture_height       = 0.30;  % 手势物理高度 (米)
TRAJ_PARA.energy_threshold_ratio = 0.4; % 能量阈值
TRAJ_PARA.min_hit_sats         = 2;     % 最少卫星数
TRAJ_PARA.miss_conflict_dist   = 0.15;  % 冲突距离 (米)
TRAJ_PARA.min_elevation        = 15;    % 最低仰角

% 绘制场景基础
plot3(ax, 0, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
plane_range = 1.5; 
patch(ax, [-plane_range plane_range plane_range -plane_range], ...
          [-plane_range -plane_range plane_range plane_range], ...
          [TRAJ_PARA.gesture_height TRAJ_PARA.gesture_height TRAJ_PARA.gesture_height TRAJ_PARA.gesture_height], ...
          [0.9 0.9 1.0], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Gesture Plane');

%% 2. 核心循环
colors = lines(length(segments)); 
fprintf('--> 开始 3D 轨迹估算 (共 %d 个片段)...\n', length(segments));

for i = 1:length(segments)
    seg = segments(i);
    idx_range = seg.start_idx : seg.end_idx;
    seg_times = t_grid(idx_range);
    
    sub_volatility = volatility_matrix(idx_range, :);
    
    % 获取接收机位置
    [~, epoch_idx] = min(abs([obs_data.time] - seg.peak_time));
    try [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, epoch_idx); catch, continue; end
    if isempty(rec_pos), continue; end
    [rec_lat, rec_lon, rec_alt] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));

    % 存储 3D 交点
    sat_points = struct('id', {}, 'pos_3d', {}, 'energy', {}, 'time_offset', {}, 'vec_u', {});
    num_pts = 0;
    
    % 计算视线交点
    for s = 1:length(valid_sats)
        s_id = valid_sats{s}; if ~isfield(sat_states, s_id), continue; end
        
        s_energy = sum(sub_volatility(:, s), 'omitnan');
        [~, local_max_idx] = max(sub_volatility(:, s));
        t_offset = seconds(seg_times(local_max_idx) - seg_times(1));
        
        sat_pos = sat_states.(s_id).position;
        [e, n, u] = ecef2enu(sat_pos(1)-rec_pos(1), sat_pos(2)-rec_pos(2), sat_pos(3)-rec_pos(3), rec_lat, rec_lon, rec_alt);
        dist = norm([e, n, u]); vec_u = [e, n, u] / dist; 
        
        if asind(vec_u(3)) < TRAJ_PARA.min_elevation || vec_u(3) <= 0, continue; end 
        
        % 射线方程求解
        t_intersect = TRAJ_PARA.gesture_height / vec_u(3);
        intersect_point = t_intersect * vec_u; 
        if norm(intersect_point(1:2)) > 2.0, continue; end
        
        num_pts = num_pts + 1;
        sat_points(num_pts).id = s_id;
        sat_points(num_pts).pos_3d = intersect_point;
        sat_points(num_pts).energy = s_energy;
        sat_points(num_pts).time_offset = t_offset;
    end
    
    if num_pts < 2, continue; end
    
    % Hit/Miss 分类
    all_energies = [sat_points.energy];
    threshold = max(all_energies) * TRAJ_PARA.energy_threshold_ratio;
    hits = sat_points(all_energies > threshold);
    misses = sat_points(all_energies <= threshold);
    
    if length(hits) < TRAJ_PARA.min_hit_sats, continue; end
    
    % PCA 拟合 (在 Z=h 平面)
    coords_3d = vertcat(hits.pos_3d); 
    P_xy = coords_3d(:, 1:2); mean_P = mean(P_xy);
    [coeff, ~, ~] = pca(P_xy - mean_P);
    dir_vec_2d = coeff(:, 1)'; 
    
    % 时序定向
    projections = (P_xy - mean_P) * dir_vec_2d';
    times = [hits.time_offset]';
    corr_val = corr(projections, times);
    if ~isnan(corr_val) && corr_val < 0, dir_vec_2d = -dir_vec_2d; end
    
    % 计算端点
    proj_final = (P_xy - mean_P) * dir_vec_2d';
    start_3d = [mean_P + (min(proj_final) - 0.1) * dir_vec_2d, TRAJ_PARA.gesture_height];
    end_3d   = [mean_P + (max(proj_final) + 0.1) * dir_vec_2d, TRAJ_PARA.gesture_height];
    
    % 冲突检测
    conflict_count = 0;
    seg_vec = end_3d - start_3d; len_sq = dot(seg_vec, seg_vec);
    for m = 1:length(misses)
        m_pt = misses(m).pos_3d;
        if len_sq > 0
            t = dot(m_pt - start_3d, seg_vec) / len_sq;
            if t > 0 && t < 1 && norm(m_pt - (start_3d + t * seg_vec)) < TRAJ_PARA.miss_conflict_dist
                conflict_count = conflict_count + 1;
                plot3(ax, m_pt(1), m_pt(2), m_pt(3), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
            end
        end
    end
    
    traj_az = atan2d(dir_vec_2d(1), dir_vec_2d(2));
    if traj_az < 0, traj_az = traj_az + 360; end
    fprintf('   Seg #%d: 3D推演方向 %.1f 度 (相关性: %.2f)\n', i, traj_az, abs(corr_val));

    % 绘图
    draw_color = colors(mod(i-1, size(colors,1)) + 1, :);
    for k = 1:length(hits)
        pt = hits(k).pos_3d;
        plot3(ax, [0, pt(1)], [0, pt(2)], [0, pt(3)], '-', 'Color', [draw_color, 0.3], 'LineWidth', 1);
        plot3(ax, pt(1), pt(2), pt(3), 'o', 'MarkerFaceColor', draw_color, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    if ~isempty(misses)
        m_pts = vertcat(misses.pos_3d);
        plot3(ax, m_pts(:,1), m_pts(:,2), m_pts(:,3), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 5);
    end
    quiver3(ax, start_3d(1), start_3d(2), start_3d(3), ...
            end_3d(1)-start_3d(1), end_3d(2)-start_3d(2), end_3d(3)-start_3d(3), ...
            0, 'Color', draw_color, 'LineWidth', 3, 'MaxHeadSize', 0.5);
    text(ax, end_3d(1), end_3d(2), end_3d(3) + 0.05, sprintf('#%d: %.0f^o', i, traj_az), ...
        'Color', draw_color, 'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', 'w');
end
hold(ax, 'off');
fprintf('✅ Step 2 分析完成 (3D 视图已生成)。\n');
end