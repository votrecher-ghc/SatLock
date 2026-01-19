% =========================================================================
% run_gesture_analysis_linear_boundary.m (重构版)
% 功能: 鲁棒手势感知 Step 2 - 边界追踪 + 线性约束 (Boundary V3 + Linear Constraint)
% 描述:
%   该函数结合了 "V3 边界追踪" 的抗遮挡能力与 "Line 线性约束" 的几何整形能力。
%   
%   核心流程:
%   1. [Step 2 - 抓点]: 采用 Boundary V3 逻辑。
%      - 计算所有波动点(Active Points)离身体中心的距离。
%      - 仅保留距离最远的 Top-K% (如 50%) 点的重心。
%      - 目的: 锁定"指尖"位置，抵抗手臂遮挡带来的重心后移。
%
%   2. [Step 2.5 - 整形]: 采用 Continuous Track Line 逻辑。
%      - 识别独立笔画 (Stroke Detection)。
%      - 对每一段轨迹进行 PCA 主成分分析，提取主轴方向。
%      - 将轨迹点强制投影到直线上。
%      - 目的: 消除手部微小抖动，输出平直、美观的几何轨迹。
%
% [调用格式]:
%   [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_linear_boundary(obs_waveform, nav_data, step1_res_shaped);
%
% [输入参数]:
%   1. obs_waveform (struct): 经过 waveform_reshaping 处理后的方波观测数据。
%   2. nav_data (struct): 导航星历数据。
%   3. step1_res_shaped (struct): Step 1 结果包。
%
% [返回值说明]:
%   1. traj_x / traj_y: 最终线性化后的轨迹坐标。
%   2. traj_t: 时间索引。
%   3. traj_e: 能量/置信度。
% =========================================================================

function [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_linear_boundary(obs_waveform, nav_data, step1_res_shaped)

% --- 1. 数据解包 ---
segments = step1_res_shaped.segments;
volatility_matrix = step1_res_shaped.volatility_matrix;
t_grid = step1_res_shaped.t_grid;
valid_sats = step1_res_shaped.valid_sats;
% PARA = step1_res_shaped.PARA;

if isempty(segments)
    fprintf('⚠️ 警告: step1_res 中无有效分段，Step 2 (Linear Boundary) 终止。\n');
    traj_x=[]; traj_y=[]; traj_t=[]; traj_e=[];
    return;
end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> [Step 2] 启动线性边界追踪 (Boundary V3 + PCA Constraint)...\n');

%% ================= [Part 1] 参数设置 =================

% 1. 几何与物理参数
TRAJ.gesture_height    = 0.20;   % [物理] 手势平面高度 (米)
TRAJ.min_action_dist   = 0.05;   % [触发] 动作死区 (米)

% 2. 边界追踪参数 (V3 Logic)
ALG.zenith_safe_deg    = 15;     % [掩膜] 最低安全仰角 (度)，保留头顶
ALG.top_k_percent      = 0.5;    % [筛选] 仅取最远的 50% 点 (指尖假设)

% 3. 线性约束参数 (Line Logic)
ALG.stroke_gap_sec     = 1.0;    % [分段] 笔画间隔时间阈值 (秒)

% 4. 时域聚类与平滑
TRAJ.time_cluster_k    = 3;      % [聚类] 时间窗口 (Step 2)
TRAJ.traj_smooth_m     = 2;      % [平滑] 最终平滑 (Step 3)

%% ================= [Part 2] 核心追踪 (Boundary V3) =================

% 1. 构建激活掩膜
num_samples = length(t_grid);
is_active_mask = false(num_samples, 1);
for k = 1:length(segments)
    is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
end

% 2. 预计算接收机位置
rec_pos_cache = NaN(num_samples, 3);
for t = 1:num_samples
    if is_active_mask(t)
        [~, epoch_idx] = min(abs([obs_waveform.time] - t_grid(t)));
        try [rp, ~, ~] = calculate_receiver_position(obs_waveform, nav_data, epoch_idx); rec_pos_cache(t,:) = rp; catch, end
    end
end

% 3. 主循环: 初始轨迹生成
K_cluster = TRAJ.time_cluster_k;
raw_track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
track_cnt = 0;
is_tracking_started = false; 

fprintf('--> 执行边界追踪 (Top %.0f%%)...\n', ALG.top_k_percent*100);

for t_start = 1 : K_cluster : num_samples
    
    t_end = min(t_start + K_cluster - 1, num_samples);
    range_indices = t_start : t_end;
    
    if ~any(is_active_mask(range_indices))
        is_tracking_started = false; 
        continue;
    end
    
    % --- A. 候选点收集 ---
    window_candidates = struct('x', {}, 'y', {}, 'dist', {}, 'weight', {});
    wc_cnt = 0;
    sum_window_energy = 0;
    
    for t = range_indices
        if ~is_active_mask(t), continue; end
        
        rec_pos = rec_pos_cache(t, :);
        if any(isnan(rec_pos)), continue; end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
        
        current_vols = volatility_matrix(t, :);
        % 方波数据: >1 即为激活 (通常为 10)
        valid_energy_idx = find(current_vols > 1.0); 
        
        if isempty(valid_energy_idx), continue; end
        
        [~, epoch_idx] = min(abs([obs_waveform.time] - t_grid(t)));
        try [~, ~, sat_states] = calculate_receiver_position(obs_waveform, nav_data, epoch_idx); catch, continue; end
        
        for k = 1:length(valid_energy_idx)
            s_idx = valid_energy_idx(k); 
            sid = valid_sats{s_idx};
            if ~isfield(sat_states, sid), continue; end
            
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
            dist = norm([e, n, u]);
            vec_u = [e, n, u] / dist; 
            
            % 仰角过滤 (保留头顶)
            el_deg = asind(vec_u(3));
            if el_deg < ALG.zenith_safe_deg, continue; end
            
            % 投影
            t_int = TRAJ.gesture_height / vec_u(3); 
            pt_int = t_int * vec_u;
            
            % 半径门控
            dist_from_center = norm(pt_int(1:2));
            if dist_from_center > 5.0, continue; end
            
            wc_cnt = wc_cnt + 1;
            window_candidates(wc_cnt).x = pt_int(1);
            window_candidates(wc_cnt).y = pt_int(2);
            window_candidates(wc_cnt).dist = dist_from_center;
            
            w = current_vols(s_idx);
            window_candidates(wc_cnt).weight = w;
            sum_window_energy = sum_window_energy + w;
        end
    end
    
    if wc_cnt == 0, continue; end
    
    % --- B. 边界筛选 (Top-K Logic) ---
    all_dists = [window_candidates.dist];
    [~, sort_idx] = sort(all_dists, 'descend'); % 从远到近排序
    
    num_to_pick = max(1, ceil(wc_cnt * ALG.top_k_percent));
    selected_idx = sort_idx(1:num_to_pick);
    
    sum_w = 0; sum_wx = 0; sum_wy = 0;
    for k = 1:length(selected_idx)
        idx = selected_idx(k);
        w = window_candidates(idx).weight;
        sum_w = sum_w + w;
        sum_wx = sum_wx + window_candidates(idx).x * w;
        sum_wy = sum_wy + window_candidates(idx).y * w;
    end
    
    if sum_w == 0, continue; end
    center_x = sum_wx / sum_w;
    center_y = sum_wy / sum_w;
    
    % --- C. 动作触发锁 ---
    if ~is_tracking_started
        if norm([center_x, center_y]) > TRAJ.min_action_dist
            is_tracking_started = true; 
        else
            continue; 
        end
    end
    
    track_cnt = track_cnt + 1;
    raw_track_results(track_cnt).t_idx = mean(range_indices);
    raw_track_results(track_cnt).x = center_x;
    raw_track_results(track_cnt).y = center_y;
    raw_track_results(track_cnt).total_energy = sum_window_energy;
end

% 提取原始轨迹
if track_cnt > 0
    traj_x = [raw_track_results.x]'; 
    traj_y = [raw_track_results.y]'; 
    traj_t = [raw_track_results.t_idx]'; 
    traj_e = [raw_track_results.total_energy]';
else
    traj_x = []; traj_y = []; traj_t = []; traj_e = [];
    fprintf('⚠️ Step 2 未生成有效轨迹。\n');
    return;
end

%% ================= [Part 3] 线性约束 (PCA Line) =================
fprintf('--> [Step 2.5] 执行轨迹线性化约束 (PCA Linear Constraint)...\n');

if ~isempty(traj_x)
    % 1. 识别笔画分段 (根据时间跳变)
    % 这里的 sampling_rate 需要估算或从 step1_res 获取
    % 假设 step1_res.t_grid 均匀，计算平均 dt
    if length(t_grid) > 1
        mean_dt = seconds(mean(diff(t_grid)));
        gap_idx_th = round(ALG.stroke_gap_sec / mean_dt);
    else
        gap_idx_th = 25; % 默认
    end
    
    time_diffs = diff(traj_t);
    break_indices = find(time_diffs > gap_idx_th);
    
    seg_starts = [1; break_indices + 1];
    seg_ends   = [break_indices; length(traj_t)];
    
    % 2. 对每一段分别进行 PCA 线性投影
    for k = 1:length(seg_starts)
        idx_s = seg_starts(k);
        idx_e = seg_ends(k);
        
        num_pts = idx_e - idx_s + 1;
        if num_pts < 3, continue; end % 点太少不拟合
        
        pts_segment = [traj_x(idx_s : idx_e), traj_y(idx_s : idx_e)];
        
        % A. PCA 主成分提取
        mu = mean(pts_segment, 1);
        pts_centered = pts_segment - mu;
        
        % SVD 分解协方差矩阵
        if size(pts_centered, 1) >= 2
            [~, ~, V] = svd(cov(pts_centered));
            main_direction = V(:, 1); % 第一主成分 (最大方差方向)
        else
            continue;
        end
        
        % B. 线性投影 (Projection)
        % Proj = (P * v) * v^T
        pts_projected = (pts_centered * main_direction) * main_direction';
        pts_final = pts_projected + mu;
        
        % C. 回填
        traj_x(idx_s : idx_e) = pts_final(:, 1);
        traj_y(idx_s : idx_e) = pts_final(:, 2);
    end
    fprintf('    已对 %d 个独立笔画进行了线性矫正。\n', length(seg_starts));
end

% 最终平滑
if TRAJ.traj_smooth_m > 1
    traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
    traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
end

%% ================= [Part 4] 绘图流程 =================
fprintf('\n--> 开始生成图表...\n');

if ~isempty(traj_x)
    figure('Name', 'Linear Boundary Track (V3 + Line)', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    % 画接收机
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--');
    
    % 画轨迹
    plot(ax, traj_x, traj_y, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Linear Path');
    
    % 标记起终点
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'End');
    
    title({'Linear Boundary Fitting', ...
           sprintf('Top-%.0f%% Boundary + PCA Line', ALG.top_k_percent*100)});
    legend('Location', 'best');
    
    max_range = max([max(abs(traj_x)), max(abs(traj_y)), 0.5]);
    xlim([-max_range*1.2, max_range*1.2]); 
    ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end

fprintf('✅ 线性边界追踪分析完成。\n');
end