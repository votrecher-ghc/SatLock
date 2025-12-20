% =========================================================================
% run_gesture_analysis_3d_separate_v3.m
% 功能: 手势感知全流程 (独立窗口 + 几何约束可视化版)
% 改进:
%   1. [Figure 1] 3D 视图：展示射线、手势平面和空间几何关系。
%   2. [Figure 2] 2D 俯视图：(新增) 专注于手势平面上的点云拟合与排斥约束验证。
%   3. [优化] 颜色加深，视觉效果更清晰。
% =========================================================================

%% 1. 环境检查与参数设置
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end

% --- [Step 1] 分段参数 ---
PARA.smooth_window_sec = 1.5;   
PARA.gvi_threshold     = 4;     
PARA.sampling_rate     = 10;    
PARA.merge_gap_sec     = 0.5;   
PARA.min_duration_sec  = 0.4;   

% --- [Step 2] 3D 轨迹参数 ---
TRAJ.gesture_height    = 0.30;  % 手势平面高度 (0.3m)
TRAJ.sky_height        = 1.0;   % 射线绘制的顶部高度 (视觉更开阔)
TRAJ.energy_threshold  = 0.4;   
TRAJ.min_hit_sats      = 2;     
TRAJ.miss_conflict_dist= 0.15;  % 排斥半径 r_eff
TRAJ.min_elevation     = 15;    

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动手势感知全流程分析 (几何约束双视图版)...\n');

%% ================= [Step 1] 数据提取、滤波与分段 =================
% (核心逻辑保持不变，确保数据处理的一致性)
fprintf('--> [Step 1] 提取全星座数据...\n');
all_sat_ids = {};
for i = 1:min(100, length(obs_data)), if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids);
valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

raw_times = [obs_data.time];
mean_dt = mean(diff(raw_times));
if mean_dt < seconds(0.05), PARA.sampling_rate = 20; elseif mean_dt < seconds(0.2), PARA.sampling_rate = 10; else, PARA.sampling_rate = 1; end

t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
t_grid_plot = t_grid + hours(8) - seconds(20); 
num_samples = length(t_grid);
num_sats = length(valid_sats);

cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    sys = sat_id(1);
    if sys=='C', codes={'S2I','S2X','S1I','S6I','S7I'}; else, codes={'S1C','S1X','S2C','S2X','S5X'}; end
    for k = 1:min(50, length(obs_data))
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            avail_codes = fieldnames(obs_data(k).data.(sat_id).snr);
            for c = 1:length(codes), if ismember(codes{c}, avail_codes), target_snr_code = codes{c}; break; end; end
            if ~isempty(target_snr_code), break; end
        end
    end
    if isempty(target_snr_code), continue; end
    s_times = []; s_cn0 = [];
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id).snr, target_snr_code)
            val = obs_data(k).data.(sat_id).snr.(target_snr_code);
            if ~isnan(val) && val > 10, s_times = [s_times; obs_data(k).time]; s_cn0 = [s_cn0; val]; end
        end
    end
    if length(s_times) > 20
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
    end
end

% 1.2 预滤波
fprintf('--> [Step 1] 执行 Savitzky-Golay 预滤波去噪...\n');
sg_order = 2; sg_len = 7;
for s = 1:num_sats
    col = cn0_matrix(:, s);
    valid = ~isnan(col);
    if sum(valid) > sg_len*2
        idx = 1:length(col);
        filled = interp1(idx(valid), col(valid), idx, 'linear', 'extrap')';
        cn0_matrix(:, s) = sgolayfilt(filled, sg_order, sg_len);
    end
end

% 1.3 分段
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth);
gvi_curve = sum(volatility_matrix, 2, 'omitnan');
gvi_curve_clean = movmean(gvi_curve, 5);
is_active = gvi_curve_clean > PARA.gvi_threshold;
gap_pts = round(PARA.merge_gap_sec * PARA.sampling_rate);
pad_act = [1; is_active; 1];
g_starts = find(diff(pad_act)==-1); g_ends = find(diff(pad_act)==1)-1;
for i=1:length(g_starts)
    if (g_ends(i)-g_starts(i)+1) < gap_pts && g_starts(i)>1
        is_active(g_starts(i):g_ends(i)-1) = 1;
    end
end
edges = diff([0; is_active; 0]);
s_idxs = find(edges==1); e_idxs = find(edges==-1)-1;
min_dur = round(PARA.min_duration_sec * PARA.sampling_rate);
segments = struct('id', {}, 'start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {});
cnt = 0;
for i=1:length(s_idxs)
    if (e_idxs(i)-s_idxs(i)) >= min_dur
        cnt = cnt + 1;
        [m_v, m_i] = max(gvi_curve_clean(s_idxs(i):e_idxs(i)));
        peak_idx = s_idxs(i) + m_i - 1;
        segments(cnt).id = cnt;
        segments(cnt).start_idx = s_idxs(i);
        segments(cnt).end_idx = e_idxs(i);
        segments(cnt).peak_idx = peak_idx;
        segments(cnt).peak_time = t_grid(peak_idx);
        segments(cnt).peak_gvi = m_v;
    end
end
fprintf('✅ [Step 1] 完成。识别到 %d 个有效手势片段。\n', cnt);

% 1.4 可视化分段
figure('Name', 'Step 1: Segmentation Summary', 'Position', [50, 500, 800, 300]);
plot(t_grid_plot, gvi_curve_clean, 'k', 'LineWidth', 1); hold on;
yline(PARA.gvi_threshold, 'b--', 'Threshold');
yl = ylim;
for i=1:length(segments)
    idx = segments(i).start_idx : segments(i).end_idx;
    plot(t_grid_plot(idx), gvi_curve_clean(idx), 'r', 'LineWidth', 2);
end
title(sprintf('手势分段概览 (阈值=%d)', PARA.gvi_threshold));
xlabel('Time (BJT)'); ylabel('GVI');
datetick('x','HH:MM:ss','keepticks','keeplimits'); grid on;

%% ================= [Step 2] 3D 轨迹推演 (独立窗口 - 几何约束可视化版) =================
fprintf('--> [Step 2] 开始 3D 空间轨迹推演 (几何反演理论版)...\n');
traj_colors = lines(length(segments));
for i = 1:length(segments)
    seg = segments(i);
    idx_range = seg.start_idx : seg.end_idx;
    seg_times = t_grid(idx_range);
    
    sub_vol = volatility_matrix(idx_range, :);
    
    [~, ep_idx] = min(abs([obs_data.time] - seg.peak_time));
    [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, ep_idx);
    if isempty(rec_pos), continue; end
    [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
    
    % 存储 3D 点及其方向向量 (vec_u)
    sat_pts = struct('pos', {}, 'vec_u', {}, 'energy', {}, 't_off', {});
    pt_cnt = 0;
    
    for s = 1:length(valid_sats)
        sid = valid_sats{s};
        if ~isfield(sat_states, sid), continue; end
        
        sat_p = sat_states.(sid).position;
        [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
        vec = [e, n, u]; dist = norm(vec); vec_u = vec/dist;
        
        if asind(vec_u(3)) < TRAJ.min_elevation, continue; end
        
        if vec_u(3) > 0
            t_int = TRAJ.gesture_height / vec_u(3);
            pt_int = t_int * vec_u;
            
            % 稍微放宽范围，避免边缘卫星被切掉
            if norm(pt_int(1:2)) < 3.0 
                pt_cnt = pt_cnt + 1;
                sat_pts(pt_cnt).pos = pt_int;
                sat_pts(pt_cnt).vec_u = vec_u; % 关键：保存单位向量用于画射线
                sat_pts(pt_cnt).energy = sum(sub_vol(:, s), 'omitnan');
                [~, mx_i] = max(sub_vol(:, s));
                sat_pts(pt_cnt).t_off = seconds(seg_times(mx_i) - seg_times(1));
            end
        end
    end
    
    if pt_cnt < 2, continue; end
    
    eners = [sat_pts.energy];
    th_e = max(eners) * TRAJ.energy_threshold;
    hits = sat_pts(eners > th_e);
    misses = sat_pts(eners <= th_e);
    
    if length(hits) < TRAJ.min_hit_sats
        fprintf('   Seg #%d: 有效卫星不足 (%d)，跳过。\n', i, length(hits));
        continue;
    end
    
    % PCA & 时序
    coords = vertcat(hits.pos);
    pts_xy = coords(:, 1:2);
    mean_xy = mean(pts_xy);
    centered = pts_xy - mean_xy;
    [coeff, ~, ~] = pca(centered);
    dir_xy = coeff(:, 1)';
    
    proj = centered * dir_xy';
    times = [hits.t_off]';
    corr_v = corr(proj, times);
    if ~isnan(corr_v) && corr_v < 0, dir_xy = -dir_xy; end
    
    proj_vals = (pts_xy - mean_xy) * dir_xy';
    p_start_2d = mean_xy + (min(proj_vals)-0.1) * dir_xy;
    p_end_2d   = mean_xy + (max(proj_vals)+0.1) * dir_xy;
    p_start = [p_start_2d, TRAJ.gesture_height];
    p_end   = [p_end_2d,   TRAJ.gesture_height];
    
    conflict = false;
    vec_seg = p_end - p_start; len_sq = dot(vec_seg, vec_seg);
    for m=1:length(misses)
        mp = misses(m).pos;
        if len_sq>0
            t = dot(mp-p_start, vec_seg)/len_sq;
            if t>0 && t<1 && norm(mp-(p_start+t*vec_seg)) < TRAJ.miss_conflict_dist
                conflict = true; break;
            end
        end
    end
    
    % ================= [绘图 1：基于“视距切割”理论的 3D 可视化] =================
    fig_name = sprintf('Gesture #%d 3D View (T=%s)', i, datestr(seg.peak_time + hours(8)-seconds(20), 'HH:MM:SS'));
    f = figure('Name', fig_name, 'Position', [100 + (i-1)*2, 100 + (i-1)*2, 1000, 800], 'Color', 'w');
    ax3d = axes('Parent', f);
    hold(ax3d, 'on'); grid(ax3d, 'on'); axis(ax3d, 'equal'); view(ax3d, 3);
    xlabel(ax3d, 'East (m)'); ylabel(ax3d, 'North (m)'); zlabel(ax3d, 'Up (m)');
    
    % --- [修复]：明确定义颜色变量 col ---
    col = traj_colors(i, :); 
    
    az = atan2d(dir_xy(1), dir_xy(2)); if az<0, az=az+360; end
    str_conflict = ""; if conflict, str_conflict = "[Miss冲突]"; end
    title(ax3d, {sprintf('手势 #%d: 拟合方向 %.1f° (时序相关性 %.2f) %s', i, az, abs(corr_v), str_conflict), ...
                 '红色实心: Hit触点(吸引) | 灰色空心: Miss触点(排斥) | 蓝色虚线: 轨迹轴线'});
    % --- 1. 设置视图范围 ---
    view_range = 3.5; 
    xlim(ax3d, [-view_range, view_range]);
    ylim(ax3d, [-view_range, view_range]);
    zlim(ax3d, [0, TRAJ.sky_height + 0.2]); 
    % --- 2. 绘制“手势投影平面” (Theory: Projection Plane) ---
    % 绘制一个半透明的平面，代表理论中的“手势高度平面 H”
    plane_z = TRAJ.gesture_height;
    patch(ax3d, [-view_range view_range view_range -view_range], ...
                [-view_range -view_range view_range view_range], ...
          [plane_z plane_z plane_z plane_z], ...
          [0.4 0.5 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Gesture Plane');
    
    % 在平面上画网格，增加空间感
    grid_step = 1.0;
    for g = -view_range:grid_step:view_range
        plot3(ax3d, [-view_range, view_range], [g, g], [plane_z, plane_z], '-', 'Color', [0.8 0.9 1], 'LineWidth', 0.5);
        plot3(ax3d, [g, g], [-view_range, view_range], [plane_z, plane_z], '-', 'Color', [0.8 0.9 1], 'LineWidth', 0.5);
    end
    text(ax3d, -view_range+0.5, -view_range+0.5, plane_z, 'Gesture Projection Plane', 'Color', 'b', 'FontSize', 8);
    % --- 3. 绘制接收机 ---
    plot3(ax3d, 0,0,0, '^', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    plot3(ax3d, [0 0], [0 0], [0 plane_z], 'k:', 'LineWidth', 0.5); % 接收机到平面的垂线
    % --- 4. 绘制“轨迹轴线” (Theory: Spatial Axis T) ---
    % 这是根据 Hit 点集拟合出的无限长直线，体现“空间共线性”
    line_len = view_range * 1.5;
    p_axis_start = [mean_xy - line_len * dir_xy, plane_z];
    p_axis_end   = [mean_xy + line_len * dir_xy, plane_z];
    plot3(ax3d, [p_axis_start(1), p_axis_end(1)], ...
                [p_axis_start(2), p_axis_end(2)], ...
                [p_axis_start(3), p_axis_end(3)], ...
                '--b', 'LineWidth', 1.5, 'DisplayName', 'Trajectory Axis');
    % --- 5. 绘制 Miss 卫星 (Theory: Exclusion Constraints) ---
    if ~isempty(misses)
        for m = 1:length(misses)
            mp = misses(m).pos; % 这是在手势平面上的投影点 P_j
            vec = misses(m).vec_u;
            
            % 画射线: 接收机 -> 投影点 -> 天空
            if vec(3) > 0
                p_top = (TRAJ.sky_height / vec(3)) * vec;
                plot3(ax3d, [0 p_top(1)], [0 p_top(2)], [0 p_top(3)], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.5);
            end
            
            % 画投影触点 (空心灰圈，表示排斥)
            plot3(ax3d, mp(1), mp(2), mp(3), 'o', ...
                  'Color', [0.6 0.6 0.6], 'MarkerSize', 4, 'DisplayName', 'Miss Contact');
        end
    end
    % --- 6. 绘制 Hit 卫星 (Theory: Attraction Constraints) ---
    for k=1:length(hits)
        hp = hits(k).pos; % 这是在手势平面上的投影点 P_k
        vec = hits(k).vec_u;
        
        if vec(3) > 0
            p_top = (TRAJ.sky_height / vec(3)) * vec;
            % --- [修复]：使用 patch 实现半透明射线 ---
            patch(ax3d, 'XData', [0 p_top(1)], 'YData', [0 p_top(2)], 'ZData', [0 p_top(3)], ...
                  'EdgeColor', col, 'EdgeAlpha', 0.6, 'LineWidth', 1.5, 'FaceColor', 'none');
        end
        
        % 画拟合误差线 (触点到轴线的距离)
        pt_on_line = mean_xy + dot(hp(1:2)-mean_xy, dir_xy) * dir_xy;
        plot3(ax3d, [hp(1), pt_on_line(1)], [hp(2), pt_on_line(2)], [hp(3), plane_z], '-', 'Color', col, 'LineWidth', 0.5);
        
        % 画投影触点 (实心彩点，表示吸引/相容)
        plot3(ax3d, hp(1), hp(2), hp(3), 'o', ...
              'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Hit Contact');
    end
    
    % --- 7. 绘制手势运动向量 (Theory: Temporal Gradient) ---
    % 在轴线上画一个箭头，表示时序推断出的方向
    quiver3(ax3d, p_start(1), p_start(2), p_start(3), ...
            p_end(1)-p_start(1), p_end(2)-p_start(2), 0, ...
            0, 'Color', 'k', 'LineWidth', 3, 'MaxHeadSize', 0.5, 'DisplayName', 'Motion Vector');
    
%% ================= [图 2: 2D 投影与排斥约束视图 (新增)] =================
    fig_name_2d = sprintf('Gesture #%d 2D Geometry Analysis', i);
    f2 = figure('Name', fig_name_2d, 'Position', [1000 + (i-1)*2, 100 + (i-1)*2, 600, 630], 'Color', 'w');
    
    % 调整坐标轴位置，为底部标签留出空间
    ax2d = axes('Parent', f2, 'Position', [0.13, 0.15, 0.775, 0.75]); 
    
    hold(ax2d, 'on'); grid(ax2d, 'on'); axis(ax2d, 'equal');
    xlabel(ax2d, 'East (m)'); ylabel(ax2d, 'North (m)');
    
    % 1. 绘制手势线段
    plot(ax2d, [p_axis_start(1), p_axis_end(1)], [p_axis_start(2), p_axis_end(2)], ...
         '--b', 'LineWidth', 1.5, 'DisplayName', 'Spatial Axis');
    % 2. 绘制接收机位置
    plot(ax2d, 0, 0, '^', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    % 3. 绘制 Miss 点
    if ~isempty(misses)
        miss_coords = vertcat(misses.pos);
        plot(ax2d, miss_coords(:,1), miss_coords(:,2), 'o', ...
             'Color', [0.4 0.4 0.4], 'MarkerSize', 5, 'LineWidth', 1, 'DisplayName', 'Miss Pts');
        theta = linspace(0, 2*pi, 30);
        for m=1:length(misses)
            cx = misses(m).pos(1) + TRAJ.miss_conflict_dist * cos(theta);
            cy = misses(m).pos(2) + TRAJ.miss_conflict_dist * sin(theta);
            plot(ax2d, cx, cy, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end
    % 4. 绘制 Hit 点
    hit_coords = vertcat(hits.pos);
    plot(ax2d, hit_coords(:,1), hit_coords(:,2), 'o', ...
         'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Hit Pts');
    % 5. 绘制手势向量
    quiver(ax2d, p_start_2d(1), p_start_2d(2), ...
           p_end_2d(1)-p_start_2d(1), p_end_2d(2)-p_start_2d(2), ...
           'Color', 'k', 'LineWidth', 2.5, 'MaxHeadSize', 0.5, 'DisplayName', 'Motion Vector', 'AutoScale', 'off');
    legend(ax2d, 'Location', 'bestoutside');
    xlim(ax2d, [-view_range, view_range]);
    ylim(ax2d, [-view_range, view_range]);
    
    % ================= [新增: 自动判断英文方向标签] =================
    % 逻辑：将角度归一化到 0-360，根据区间判断东南西北
    % 假设坐标系为 ENU (East=0, North=90)
    az_norm = mod(az, 360); 
    if (az_norm >= 315 || az_norm < 45)
        dir_str = 'West';
    elseif (az_norm >= 45 && az_norm < 135)
        dir_str = 'South';
    elseif (az_norm >= 135 && az_norm < 225)
        dir_str = 'East';
    else
        dir_str = 'North';
    end
    
    % 生成标签字符串: (a) Moving East
    fig_label_str = sprintf('\\bf(a) Moving %s', dir_str);
    
    annotation(f2, 'textbox',...
        [0.1, 0.01, 0.8, 0.05],...
        'String', fig_label_str,... % <--- 使用自动生成的英文标签
        'FontSize', 14,...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle',...
        'LineStyle', 'none',...
        'BackgroundColor', 'w');
    % =======================================================

    fprintf('   Seg #%d: 方向 %.1f 度 (Moving %s) (相关性 %.2f) %s\n', i, az, dir_str, abs(corr_v), str_conflict);
end
fprintf('✅ 全流程分析完成！\n');