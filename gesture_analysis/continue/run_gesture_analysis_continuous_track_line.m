% =========================================================================
% run_gesture_analysis_continuous_track_line.m (函数版)
% 功能: 鲁棒手势感知 Step 2 - 线性约束增强版 (v2.5 Struct Input)
% 描述:
%   该函数是手势分析流水线的 Step 2 变体，专门针对规则几何形状（如 A, Z, N）优化。
%   它接收 Step 1 处理后的干净数据，在"强仰角加权追踪"的基础上，增加了
%   "PCA 线性约束" (PCA Constraint) 后处理步骤。
%   
%   它能自动识别独立的笔画分段，并利用主成分分析 (SVD/PCA) 提取每一笔的
%   主方向，强制将轨迹点投影到直线上，从而消除手部悬空书写时的微小抖动，
%   使绘制出的线条更加平直、美观。
%
% [调用格式]:
%   [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_continuous_track_line(obs_clean, nav_data, step1_res);
%
% [输入参数]:
%   1. obs_clean (struct): 
%      Step 1 返回的清洗后观测数据。用于提供伪距以计算接收机位置。
%   2. nav_data (struct): 
%      导航星历数据。用于计算卫星仰角和方位角。
%   3. step1_res (struct): 
%      Step 1 返回的结果包。包含核心的 .volatility_matrix (能量权重) 
%      以及 .segments (分段信息) 和 .t_grid (时间轴)。
%
% [返回值说明]:
%   1. traj_x / traj_y (double列向量): 
%      [线性矫正后] 的轨迹坐标 (East/North)，单位: 米。
%   2. traj_t (double列向量): 
%      轨迹点对应的时间索引。
%   3. traj_e (double列向量): 
%      轨迹点的能量权重总和 (置信度)。
%
% [核心算法]:
%   1. Zenith Focus: Weight = Energy * (sin(Elevation))^N，仅信任头顶卫星。
%   2. PCA Straightening: 对每一笔画点云做奇异值分解，提取第一主成分作为轴线进行投影。
% =========================================================================

function [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_continuous_track_line(obs_clean, nav_data, step1_res)

% --- 1. 数据解包 (Unpack Data) ---
segments = step1_res.segments;
volatility_matrix = step1_res.volatility_matrix;
t_grid = step1_res.t_grid;
valid_sats = step1_res.valid_sats;
PARA = step1_res.PARA; 

if isempty(segments)
    fprintf('⚠️ 警告: step1_res 中无有效分段，Step 2 (Linear) 终止。\n');
    traj_x=[]; traj_y=[]; traj_t=[]; traj_e=[];
    return;
end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> [Step 2] 启动线性约束追踪 (Function版 v2.5: PCA Linear Constraint)...\n');

%% ================= [Part 1] 参数设置 =================

% 1. 轨迹重建参数
TRAJ.gesture_height    = 0.30;   % [物理] 手势平面的假设高度 (米)
TRAJ.min_elevation     = 0;      % [物理] 最低仰角门限 (度): 靠加权压制
TRAJ.min_action_dist   = 0.05;   % [触发] 动作死区 (米)

% 2. 过滤参数
PARA.min_sat_vol       = 1.0;    % [过滤] 单星波动门限 (dB)

% 3. 时域聚类参数
TRAJ.time_cluster_k    = 5;      % [聚类] 时间窗口: 每 5 个原始采样点合并为一个轨迹点
TRAJ.traj_smooth_m     = 1;      % [平滑] 最终轨迹的平滑窗口大小

% 4. 仰角加权策略 (Strategy: Zenith Focus)
TRAJ.elevation_power   = 2;      % [加权] 仰角权重指数。越高越只信头顶信号。

%% ================= [Part 2] 核心计算流程 =================

% 1. 构建激活掩膜 (Mask)
%    为了保持轨迹连续性，我们将所有分段的区间标记为 Active
num_samples = length(t_grid);
is_active_mask = false(num_samples, 1);
for k = 1:length(segments)
    is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
end

% 2. 预计算接收机位置缓存 (加速)
fprintf('    预计算接收机坐标...\n');
rec_pos_cache = NaN(num_samples, 3);
for t = 1:num_samples
    if is_active_mask(t)
        % 使用 t_grid 对齐 obs_clean 时间戳
        [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
        try [rp,~,~]=calculate_receiver_position(obs_clean, nav_data, epoch_idx); rec_pos_cache(t,:)=rp; catch, end
    end
end

% 3. 主循环: 时间窗口步进
track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
track_cnt = 0;
is_tracking_started = false; 

K_cluster = TRAJ.time_cluster_k;
fprintf('--> 执行强仰角加权追踪 (Power=%d)...\n', TRAJ.elevation_power);

for t_start = 1 : K_cluster : num_samples
    
    t_end = min(t_start + K_cluster - 1, num_samples);
    range_indices = t_start : t_end;
    
    if ~any(is_active_mask(range_indices))
        is_tracking_started = false; 
        continue;
    end
    
    % --- A. 收集阶段 ---
    cluster_pts_x = [];
    cluster_pts_y = [];
    cluster_weights = [];
    
    for t = range_indices
        if ~is_active_mask(t), continue; end
        
        rec_pos = rec_pos_cache(t, :); 
        if any(isnan(rec_pos)), continue; end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
        
        % [关键] 直接使用 Step 1 传入的能量矩阵
        current_vols = volatility_matrix(t, :);
        valid_energy_idx = find(current_vols > PARA.min_sat_vol); 
        
        if isempty(valid_energy_idx), continue; end
        
        % 计算卫星位置
        [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
        try [~, ~, sat_states] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); catch, continue; end
        
        for k = 1:length(valid_energy_idx)
            s_idx = valid_energy_idx(k); 
            sid = valid_sats{s_idx};
            if ~isfield(sat_states, sid), continue; end
            
            % 几何投影
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
            vec_u = [e, n, u] / norm([e, n, u]); 
            
            % 仰角检查
            el_rad = asin(vec_u(3)); el_deg = rad2deg(el_rad);
            if vec_u(3) <= 0 || el_deg < TRAJ.min_elevation, continue; end 
            
            t_int = TRAJ.gesture_height / vec_u(3); 
            pt = t_int * vec_u;
            if norm(pt(1:2)) > 5.0, continue; end 
            
            % === 强仰角加权计算 ===
            base_energy = current_vols(s_idx);
            w_elevation = (sin(el_rad)) ^ TRAJ.elevation_power;
            final_w = base_energy * w_elevation;
            
            cluster_pts_x(end+1) = pt(1);
            cluster_pts_y(end+1) = pt(2);
            cluster_weights(end+1) = final_w;
        end
    end
    
    % --- B. 聚合阶段 ---
    if isempty(cluster_weights), continue; end
    sum_w = sum(cluster_weights);
    if sum_w == 0, continue; end
    
    center_x = sum(cluster_pts_x .* cluster_weights) / sum_w;
    center_y = sum(cluster_pts_y .* cluster_weights) / sum_w;
    
    % --- C. 动作触发锁逻辑 ---
    dist_from_origin = norm([center_x, center_y]);
    if ~is_tracking_started
        if dist_from_origin > TRAJ.min_action_dist, is_tracking_started = true; else, continue; end
    end
    
    track_cnt = track_cnt + 1;
    track_results(track_cnt).t_idx = mean(range_indices);
    track_results(track_cnt).x = center_x;
    track_results(track_cnt).y = center_y;
    track_results(track_cnt).total_energy = sum_w; 
end

% 轨迹提取
if track_cnt > 0
    traj_x = [track_results.x]'; 
    traj_y = [track_results.y]'; 
    traj_t = [track_results.t_idx]'; 
    traj_e = [track_results.total_energy]';
else
    traj_x = []; traj_y = []; traj_t = []; traj_e = [];
end

% ================= [Step 2.5] 轨迹线性化约束 (PCA Straightening) =================
if ~isempty(traj_x)
    fprintf('--> [Step 2.5] 执行轨迹线性化约束 (PCA Linear Constraint)...\n');
    
    % 1. 识别笔画分段 (根据时间跳变)
    stroke_gap_threshold = 1.0 * PARA.sampling_rate; % 1秒间隔
    
    time_diffs = diff(traj_t);
    break_indices = find(time_diffs > stroke_gap_threshold);
    
    seg_starts = [1; break_indices + 1];
    seg_ends   = [break_indices; length(traj_t)];
    
    % 2. 对每一段分别进行 PCA 线性投影
    for k = 1:length(seg_starts)
        idx_s = seg_starts(k);
        idx_e = seg_ends(k);
        
        num_pts = idx_e - idx_s + 1;
        if num_pts < 3, continue; end 
        
        pts_segment = [traj_x(idx_s : idx_e), traj_y(idx_s : idx_e)];
        
        % A. PCA 主成分提取
        mu = mean(pts_segment, 1);
        pts_centered = pts_segment - mu;
        [~, ~, V] = svd(cov(pts_centered));
        main_direction = V(:, 1); 
        
        % B. 线性投影约束 (Projection)
        pts_projected = (pts_centered * main_direction) * main_direction';
        pts_final = pts_projected + mu;
        
        % C. 回填
        traj_x(idx_s : idx_e) = pts_final(:, 1);
        traj_y(idx_s : idx_e) = pts_final(:, 2);
    end
    fprintf('    已对 %d 个独立笔画进行了线性矫正。\n', length(seg_starts));
end

% 轨迹平滑
if ~isempty(traj_x) && TRAJ.traj_smooth_m > 1
    traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
    traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
end

%% ================= [Part 3] 绘图流程 =================
fprintf('\n--> 开始生成图表...\n');

if ~isempty(traj_x)
    figure('Name', 'Trajectory v2.5 (Linear Constraint)', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    % 画接收机位置
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--');
    
    % 画轨迹
    plot(ax, traj_x, traj_y, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'PCA Path');
    
    % 标记起终点
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'End');
    
    title({sprintf('v2.5 线性约束增强版 (PCA Constraint)'), ...
           sprintf('El Power=%d, MinVol=%.1f', TRAJ.elevation_power, PARA.min_sat_vol)});
    legend('Location', 'best');
    
    max_range = max([abs(traj_x); abs(traj_y); 0.5]);
    xlim([-max_range*1.2, max_range*1.2]); 
    ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end

fprintf('✅ v2.5 分析完成 (已返回线性矫正后的轨迹数据)。\n');
end