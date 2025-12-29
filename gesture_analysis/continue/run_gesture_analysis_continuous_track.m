% =========================================================================
% run_gesture_analysis_continuous_track_line.m (函数版)
% 功能: 鲁棒手势感知 Step 2 - 线性约束增强版 (v2.5)
% 描述:
%   该函数是 Step 2 的一种实现，接收 Step 1 处理好的干净数据。
%   它在强仰角加权追踪 (Zenith Focus) 的基础上，增加了后处理线性化约束。
%   利用 PCA (主成分分析) 识别每一笔画的主方向，并强制将轨迹点投影到
%   该主轴上，专门用于消除手部微小抖动，绘制规则几何图形 (如 A, N, Z)。
%
% [调用格式]:
%   [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_continuous_track_line(obs_clean, nav_data, step1_res);
%
% [输入参数]:
%   1. obs_clean (struct): Step 1 返回的清洗后观测数据 (用于计算接收机位置)。
%   2. nav_data (struct): 导航星历数据 (用于计算卫星仰角/方位角)。
%   3. step1_res (struct): Step 1 返回的结果包 (包含能量矩阵、分段信息等)。
%
% [返回值说明]:
%   1. traj_x / traj_y (double列向量): [线性矫正后] 的轨迹坐标 (米)。
%   2. traj_t (double列向量): 轨迹点对应的时间索引。
%   3. traj_e (double列向量): 轨迹点的能量权重。
%
% [核心算法]:
%   1. Elevation Weighting: Weight = Energy * (sin(El))^4，强力压制低仰角干扰。
%   2. PCA Projection: 对提取出的独立笔画进行 SVD 分解，投影至第一主成分方向。
% =========================================================================

function [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_continuous_track(obs_clean, nav_data, step1_res)

% --- 1. 解包数据 ---
segments = step1_res.segments;
volatility_matrix = step1_res.volatility_matrix;
t_grid = step1_res.t_grid;
valid_sats = step1_res.valid_sats;
% PARA = step1_res.PARA; 

if isempty(segments)
    fprintf('⚠️ 警告: 未收到有效手势分段，Step 2 终止。\n');
    traj_x=[]; traj_y=[]; traj_t=[]; traj_e=[];
    return;
end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> [Step 2] 启动强仰角加权追踪 (Continuous Track)...\n');

%% 2. 追踪参数
TRAJ.gesture_height    = 0.30;   
TRAJ.min_elevation     = 0;      % 不设硬阈值，靠加权压制
TRAJ.min_action_dist   = 0.05;   
TRAJ.elevation_power   = 2;      % [核心] 仰角权重指数 (Zenith Focus)
TRAJ.time_cluster_k    = 5;      
TRAJ.traj_smooth_m     = 2; 
PARA.min_sat_vol       = 1.0;    % 最小波动门限

%% 3. 核心追踪循环
% 由于 Step 1 已经给了分段，我们可以选择只计算分段内的，也可以计算全量的。
% 这里的 continuous_track 通常是针对整个时间轴或者激活区间的。
% 为了保持连续性，我们利用 step1 计算的 GVI 掩膜或者分段来决定计算范围。

% 重新构建一个全局 mask，或者直接遍历所有分段
num_samples = length(t_grid);
is_active_mask = false(num_samples, 1);
for k=1:length(segments)
    is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
end

track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
track_cnt = 0;
is_tracking_started = false;

% 预计算接收机位置 (使用 obs_clean)
rec_pos_cache = NaN(num_samples, 3);
for t = 1:num_samples
    if is_active_mask(t)
        % 使用 t_grid 对应的时间去 obs_clean 里找最近的观测
        [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
        try [rp,~,~]=calculate_receiver_position(obs_clean, nav_data, epoch_idx); rec_pos_cache(t,:)=rp; catch, end
    end
end

% 聚类循环
K_cluster = TRAJ.time_cluster_k;
for t_start = 1 : K_cluster : num_samples
    
    t_end = min(t_start + K_cluster - 1, num_samples);
    range_indices = t_start : t_end;
    
    if ~any(is_active_mask(range_indices)), is_tracking_started = false; continue; end
    
    cluster_pts_x = []; cluster_pts_y = []; cluster_weights = [];
    
    for t = range_indices
        if ~is_active_mask(t), continue; end
        
        rec_pos = rec_pos_cache(t, :); 
        if any(isnan(rec_pos)), continue; end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
        
        % 直接使用 step1 传来的能量矩阵
        current_vols = volatility_matrix(t, :);
        valid_energy_idx = find(current_vols > PARA.min_sat_vol);
        
        if isempty(valid_energy_idx), continue; end
        
        [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
        try [~, ~, sat_states] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); catch, continue; end
        
        for k = 1:length(valid_energy_idx)
            s_idx = valid_energy_idx(k); 
            sid = valid_sats{s_idx};
            if ~isfield(sat_states, sid), continue; end
            
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
            vec_u = [e, n, u] / norm([e, n, u]); 
            
            % 仰角计算
            el_rad = asin(vec_u(3)); el_deg = rad2deg(el_rad);
            if vec_u(3) <= 0 || el_deg < TRAJ.min_elevation, continue; end
            
            t_int = TRAJ.gesture_height / vec_u(3); 
            pt = t_int * vec_u;
            if norm(pt(1:2)) > 5.0, continue; end
            
            % === 强仰角加权 ===
            base_energy = current_vols(s_idx);
            w_elevation = (sin(el_rad)) ^ TRAJ.elevation_power;
            final_w = base_energy * w_elevation;
            
            cluster_pts_x(end+1) = pt(1);
            cluster_pts_y(end+1) = pt(2);
            cluster_weights(end+1) = final_w;
        end
    end
    
    if isempty(cluster_weights), continue; end
    sum_w = sum(cluster_weights);
    if sum_w == 0, continue; end
    
    center_x = sum(cluster_pts_x .* cluster_weights) / sum_w;
    center_y = sum(cluster_pts_y .* cluster_weights) / sum_w;
    
    % 动作触发锁
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

% 提取结果
if track_cnt > 0
    traj_x = [track_results.x]'; 
    traj_y = [track_results.y]'; 
    traj_t = [track_results.t_idx]'; 
    traj_e = [track_results.total_energy]';
    
    if TRAJ.traj_smooth_m > 1
        traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
        traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
    end
else
    traj_x = []; traj_y = []; traj_t = []; traj_e = [];
end

%% 4. 绘图 (保留原版绘图逻辑)
if ~isempty(traj_x)
    figure('Name', 'Continuous Track (Step 2 Result)', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--');
    plot(ax, traj_x, traj_y, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r');
    
    title(sprintf('Continuous Track (Power=%d)', TRAJ.elevation_power));
    
    max_range = max([abs(traj_x); abs(traj_y); 0.5]);
    xlim([-max_range*1.2, max_range*1.2]); 
    ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ Step 2 未生成有效轨迹。\n');
end

fprintf('✅ Step 2 追踪完成。\n');
end