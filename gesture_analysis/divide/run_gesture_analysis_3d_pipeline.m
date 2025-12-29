% =========================================================================
% run_gesture_analysis_3d_pipeline.m (函数版)
% 功能: 手势感知全流程 - 3D 几何反演与可视化版
% 核心逻辑:
%   1. [分段]: 基于 GVI 能量波动提取手势动作区间。
%   2. [几何反演]: "Hit/Miss" 视距切割理论。
%      - Hit (触点): 能量下降的卫星，产生"吸引"约束。
%      - Miss (未触点): 能量稳定的卫星，产生"排斥"约束。
%   3. [PCA 拟合]: 在手势平面上拟合点云主方向 (轴线)。
%   4. [双视图可视化]: 
%      - Figure 1: 3D 空间视图 (射线、平面、轴线)。
%      - Figure 2: 2D 平面视图 (排斥约束验证、方向标签)。
%
% [调用格式]:
%   [analysis_results, segments] = run_gesture_analysis_3d_pipeline(obs_data, nav_data);
%
% [返回值说明]:
%   1. analysis_results (struct数组): 几何分析详情。
%      - .id: 手势序号
%      - .p_start / .p_end: 拟合出的 3D 起终点坐标 (米)
%      - .azimuth: 运动方位角 (度)
%      - .direction: 文字描述 (如 "Moving East")
%      - .conflict: 是否违反排斥约束 (布尔值)
%   2. segments (struct数组): GVI 分段信息。
% =========================================================================

function [analysis_results, segments] = run_gesture_analysis_3d_pipeline(obs_data, nav_data)

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动手势感知全流程分析 (Function版: 3D Geometry Pipeline)...\n');

%% 1. 参数设置
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

%% ================= [Step 1] 数据提取、滤波与分段 =================
fprintf('--> [Step 1] 提取全星座数据...\n');

% 1.1 提取有效卫星
all_sat_ids = {};
for i = 1:min(100, length(obs_data)), if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids);
valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

% 自动调整采样率
raw_times = [obs_data.time];
mean_dt = mean(diff(raw_times));
if mean_dt < seconds(0.05), PARA.sampling_rate = 20; elseif mean_dt < seconds(0.2), PARA.sampling_rate = 10; else, PARA.sampling_rate = 1; end
fprintf('    自动匹配采样率: %d Hz\n', PARA.sampling_rate);

t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
t_grid_plot = t_grid + hours(8) - seconds(20); 
num_samples = length(t_grid);
num_sats = length(valid_sats);

% 提取 CN0 数据
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

% 1.2 预滤波 (SG Filter)
fprintf('--> [Step 1] 执行 Savitzky-Golay 预滤波去噪...\n');
sg_order = 2; sg_len = 7;
for s = 1:num_sats
    col = cn0_matrix(:, s); valid = ~isnan(col);
    if sum(valid) > sg_len*2
        idx = 1:length(col);
        filled = interp1(idx(valid), col(valid), idx, 'linear', 'extrap')';
        cn0_matrix(:, s) = sgolayfilt(filled, sg_order, sg_len);
    end
end

% 1.3 GVI 计算与分段
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth);
gvi_curve = sum(volatility_matrix, 2, 'omitnan');
gvi_curve_clean = movmean(gvi_curve, 5);
is_active = gvi_curve_clean > PARA.gvi_threshold;

% 缝合间隙
gap_pts = round(PARA.merge_gap_sec * PARA.sampling_rate);
pad_act = [1; is_active; 1];
g_starts = find(diff(pad_act)==-1); g_ends = find(diff(pad_act)==1)-1;
for i=1:length(g_starts)
    if (g_ends(i)-g_starts(i)+1) < gap_pts && g_starts(i)>1
        is_active(g_starts(i):g_ends(i)-1) = 1;
    end
end

% 提取有效片段
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

% 1.4 可视化分段概览
if cnt > 0
    figure('Name', 'Step 1: Segmentation Summary', 'Position', [50, 500, 800, 300], 'Color', 'w');
    plot(t_grid_plot, gvi_curve_clean, 'k', 'LineWidth', 1); hold on;
    yline(PARA.gvi_threshold, 'b--', 'Threshold');
    for i=1:length(segments)
        idx = segments(i).start_idx : segments(i).end_idx;
        plot(t_grid_plot(idx), gvi_curve_clean(idx), 'r', 'LineWidth', 2);
    end
    title(sprintf('手势分段概览 (阈值=%d)', PARA.gvi_threshold));
    xlabel('Time (BJT)'); ylabel('GVI');
    datetick('x','HH:MM:ss','keepticks','keeplimits'); grid on;
end

%% ================= [Step 2] 3D 轨迹推演 (独立窗口 - 几何约束可视化版) =================
fprintf('--> [Step 2] 开始 3D 空间轨迹推演 (几何反演理论版)...\n');

analysis_results = struct('id', {}, 'p_start', {}, 'p_end', {}, 'azimuth', {}, 'direction', {}, 'conflict', {});
traj_colors = lines(length(segments));

for i = 1:length(segments)
    seg = segments(i);
    idx_range = seg.start_idx : seg.end_idx;
    seg_times = t_grid(idx_range);
    sub_vol = volatility_matrix(idx_range, :);
    
    % 计算参考接收机位置
    [~, ep_idx] = min(abs([obs_data.time] - seg.peak_time));
    try [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, ep_idx); catch, continue; end
    if isempty(rec_pos), continue; end
    [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
    
    % 存储 3D 点及其方向向量
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
                sat_pts(pt_cnt).vec_u = vec_u; 
                sat_pts(pt_cnt).energy = sum(sub_vol(:, s), 'omitnan');
                [~, mx_i] = max(sub_vol(:, s));
                sat_pts(pt_cnt).t_off = seconds(seg_times(mx_i) - seg_times(1));
            end
        end
    end
    
    if pt_cnt < 2, continue; end
    
    % 区分 Hit 和 Miss
    eners = [sat_pts.energy];
    th_e = max(eners) * TRAJ.energy_threshold;
    hits = sat_pts(eners > th_e);
    misses = sat_pts(eners <= th_e);
    
    if length(hits) < TRAJ.min_hit_sats
        fprintf('   Seg #%d: 有效卫星不足 (%d)，跳过。\n', i, length(hits));
        continue;
    end
    
    % PCA & 时序相关性分析
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
    
    % 计算轴线端点
    proj_vals = (pts_xy - mean_xy) * dir_xy';
    p_start_2d = mean_xy + (min(proj_vals)-0.1) * dir_xy;
    p_end_2d   = mean_xy + (max(proj_vals)+0.1) * dir_xy;
    p_start = [p_start_2d, TRAJ.gesture_height];
    p_end   = [p_end_2d,   TRAJ.gesture_height];
    
    % 排斥冲突检测
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
    
    % 计算方位角与方向描述
    az = atan2d(dir_xy(1), dir_xy(2)); if az<0, az=az+360; end
    az_norm = mod(az, 360); 
    if (az_norm >= 315 || az_norm < 45), dir_str = 'West'; % ENU系: x=East, y=North. 0度为East? 不, atan2(x,y) 0度为North(y轴), 90度为East(x轴). 
    % 修正: atan2d(x,y): x是East, y是North.
    % 0度(x=0,y=1) -> North. 90度(x=1,y=0) -> East. 
    % 45-135 -> East. 135-225 -> South. 225-315 -> West. 315-45 -> North.
        if (az_norm >= 315 || az_norm < 45), dir_str = 'North';
        elseif (az_norm >= 45 && az_norm < 135), dir_str = 'East';
        elseif (az_norm >= 135 && az_norm < 225), dir_str = 'South';
        else, dir_str = 'West';
        end
    else % 使用原始代码逻辑
         if (az_norm >= 315 || az_norm < 45), dir_str = 'North'; % 原代码可能有误，此处按标准ENU北向为0修正? 
         % 暂保持与您提供的脚本逻辑一致，假设 azimuth 定义为 standard geographic heading? 
         % 若 atan2d(e, n), 则 0=North, 90=East.
         % 您的代码: atan2d(dir_xy(1), dir_xy(2)) -> atan2d(e, n). Correct.
         end
    end
    
    % 记录结果
    analysis_results(i).id = i;
    analysis_results(i).p_start = p_start;
    analysis_results(i).p_end = p_end;
    analysis_results(i).azimuth = az;
    analysis_results(i).direction = dir_str;
    analysis_results(i).conflict = conflict;
    
    str_conflict = ""; if conflict, str_conflict = "[Miss冲突]"; end
    fprintf('   Seg #%d: 方向 %.1f 度 (Moving %s) (相关性 %.2f) %s\n', i, az, dir_str, abs(corr_v), str_conflict);

    % ================= [绘图 1：3D 视图] =================
    fig_name = sprintf('Gesture #%d 3D View (T=%s)', i, datestr(seg.peak_time + hours(8)-seconds(20), 'HH:MM:SS'));
    f = figure('Name', fig_name, 'Position', [100 + (i-1)*2, 100 + (i-1)*2, 1000, 800], 'Color', 'w');
    ax3d = axes('Parent', f);
    hold(ax3d, 'on'); grid(ax3d, 'on'); axis(ax3d, 'equal'); view(ax3d, 3);
    xlabel(ax3d, 'East (m)'); ylabel(ax3d, 'North (m)'); zlabel(ax3d, 'Up (m)');
    
    col = traj_colors(i, :); 
    title(ax3d, {sprintf('手势 #%d: 拟合方向 %.1f° (时序相关性 %.2f) %s', i, az, abs(corr_v), str_conflict), ...
                 '红色实心: Hit触点(吸引) | 灰色空心: Miss触点(排斥) | 蓝色虚线: 轨迹轴线'});
    
    % 视图范围
    view_range = 3.5; 
    xlim(ax3d, [-view_range, view_range]); ylim(ax3d, [-view_range, view_range]); zlim(ax3d, [0, TRAJ.sky_height + 0.2]); 
    
    % 绘制平面
    plane_z = TRAJ.gesture_height;
    patch(ax3d, [-view_range view_range view_range -view_range], ...
                [-view_range -view_range view_range view_range], ...
          [plane_z plane_z plane_z plane_z], ...
          [0.4 0.5 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Gesture Plane');
    
    % 绘制轴线
    line_len = view_range * 1.5;
    p_axis_start = [mean_xy - line_len * dir_xy, plane_z];
    p_axis_end   = [mean_xy + line_len * dir_xy, plane_z];
    plot3(ax3d, [p_axis_start(1), p_axis_end(1)], [p_axis_start(2), p_axis_end(2)], [p_axis_start(3), p_axis_end(3)], ...
                '--b', 'LineWidth', 1.5, 'DisplayName', 'Trajectory Axis');
    
    % 绘制接收机
    plot3(ax3d, 0,0,0, '^', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    
    % 绘制 Miss 点
    if ~isempty(misses)
        for m = 1:length(misses)
            mp = misses(m).pos; vec = misses(m).vec_u;
            if vec(3) > 0, p_top = (TRAJ.sky_height / vec(3)) * vec; plot3(ax3d, [0 p_top(1)], [0 p_top(2)], [0 p_top(3)], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.5); end
            plot3(ax3d, mp(1), mp(2), mp(3), 'o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 4);
        end
    end
    
    % 绘制 Hit 点
    for k=1:length(hits)
        hp = hits(k).pos; vec = hits(k).vec_u;
        if vec(3) > 0, p_top = (TRAJ.sky_height / vec(3)) * vec; patch(ax3d, 'XData', [0 p_top(1)], 'YData', [0 p_top(2)], 'ZData', [0 p_top(3)], 'EdgeColor', col, 'EdgeAlpha', 0.6, 'LineWidth', 1.5, 'FaceColor', 'none'); end
        plot3(ax3d, hp(1), hp(2), hp(3), 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
    end
    
    % 绘制运动向量
    quiver3(ax3d, p_start(1), p_start(2), p_start(3), p_end(1)-p_start(1), p_end(2)-p_start(2), 0, 0, 'Color', 'k', 'LineWidth', 3, 'MaxHeadSize', 0.5);
    
    % ================= [绘图 2：2D 视图] =================
    fig_name_2d = sprintf('Gesture #%d 2D Geometry Analysis', i);
    f2 = figure('Name', fig_name_2d, 'Position', [1000 + (i-1)*2, 100 + (i-1)*2, 600, 630], 'Color', 'w');
    ax2d = axes('Parent', f2, 'Position', [0.13, 0.15, 0.775, 0.75]); 
    hold(ax2d, 'on'); grid(ax2d, 'on'); axis(ax2d, 'equal');
    ylabel(ax2d, 'North (m)');
    
    plot(ax2d, [p_axis_start(1), p_axis_end(1)], [p_axis_start(2), p_axis_end(2)], '--b', 'LineWidth', 1.5);
    plot(ax2d, 0, 0, '^', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    
    if ~isempty(misses)
        miss_coords = vertcat(misses.pos);
        plot(ax2d, miss_coords(:,1), miss_coords(:,2), 'o', 'Color', [0.4 0.4 0.4], 'MarkerSize', 5);
        theta = linspace(0, 2*pi, 30);
        for m=1:length(misses)
            cx = misses(m).pos(1) + TRAJ.miss_conflict_dist * cos(theta);
            cy = misses(m).pos(2) + TRAJ.miss_conflict_dist * sin(theta);
            plot(ax2d, cx, cy, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
        end
    end
    
    hit_coords = vertcat(hits.pos);
    plot(ax2d, hit_coords(:,1), hit_coords(:,2), 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    quiver(ax2d, p_start_2d(1), p_start_2d(2), p_end_2d(1)-p_start_2d(1), p_end_2d(2)-p_start_2d(2), 'Color', 'k', 'LineWidth', 2.5, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
    
    xlim(ax2d, [-view_range, view_range]); ylim(ax2d, [-view_range, view_range]);
    
    fig_label_str = sprintf('\\bf(a) Moving %s', dir_str);
    title(ax2d, sprintf('Gesture #%d 2D Geometry', i), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel(ax2d, fig_label_str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    text(ax2d, max(xlim), min(ylim), 'East (m) ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    plotedit(f2, 'on');
end

fprintf('✅ 全流程分析完成 (已返回几何分析数据)。\n');
end