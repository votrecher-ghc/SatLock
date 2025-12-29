% =========================================================================
% run_gesture_analysis_boundary_track.m (函数版)
% 功能: 鲁棒手势感知 Step 2 - 稳定对抗边界追踪版 (v5.1 Struct Input)
% 描述:
%   该函数是 Step 2 的一种变体，采用"稳定-波动对抗" (Stable-Active Opposition) 算法。
%   它将可见卫星分为两组：
%     1. 波动组 (Active): 能量高，代表被遮挡(通常是手臂/手掌)。
%     2. 稳定组 (Stable): 能量低，代表背景或未遮挡区域。
%   算法计算从"波动重心"指向"稳定重心"的对抗矢量，并在此方向上筛选
%   最前沿的波动点 (Top-K%)，以此来抵抗手臂对重心的拖拽效应。
%
% [调用格式]:
%   [traj_x, traj_y, traj_t, traj_counts] = run_gesture_analysis_boundary_track(obs_clean, nav_data, step1_res);
%
% [输入参数]:
%   1. obs_clean (struct): Step 1 返回的清洗后观测数据 (用于定位)。
%   2. nav_data (struct): 导航星历数据。
%   3. step1_res (struct): Step 1 返回的结果包。
%
% [返回值说明]:
%   1. traj_x / traj_y (double): 轨迹坐标 (米)。
%   2. traj_t (double): 轨迹点时间索引。
%   3. traj_counts (double): 
%      "波动卫星计数" (Active Sat Count)。
%      该值越大，代表遮挡越显著，点越可信。
%
% [核心算法]:
%   1. Classification: 基于瞬时能量阈值将卫星分为 Active/Stable 两类。
%   2. Opposition Vector: Vec = Centroid(Stable) - Centroid(Active)。
%   3. Top-K Projection: 在 Vec 方向上投影并提取前沿点。
% =========================================================================

function [traj_x, traj_y, traj_t, traj_active_counts] = run_gesture_analysis_boundary_track(obs_clean, nav_data, step1_res)

% --- 1. 数据解包 ---
segments = step1_res.segments;
volatility_matrix = step1_res.volatility_matrix;
t_grid = step1_res.t_grid;
valid_sats = step1_res.valid_sats;
% PARA = step1_res.PARA;

if isempty(segments)
    fprintf('⚠️ 警告: step1_res 中无有效分段，Step 2 (Opposition) 终止。\n');
    traj_x=[]; traj_y=[]; traj_t=[]; traj_active_counts=[];
    return;
end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> [Step 2] 启动稳定对抗追踪 (Function版 v5.1: Stable-Active Opposition)...\n');

%% ================= [Part 1] 参数设置 =================

% 1. 轨迹与几何参数
TRAJ.gesture_height    = 0.20;  % [物理] 手势高度 (米)
TRAJ.min_elevation     = 15;    % [物理] 最低仰角 (度)
TRAJ.min_action_dist   = 0.05;  % [触发] 动作死区 (米)

% 2. 算法特有参数
% [核心] 单星波动判定阈值: 高于此值视为 Active，低于视为 Stable
% 注意: 这里的能量是经过 SG 滤波后的纯净波动值
PARA.active_th         = 1.0;   

% 3. 对抗与边界参数
ALG.zenith_safe_deg    = 30;    % [掩膜] 安全仰角门限 (度)
ALG.top_k_percent      = 0.2;   % [核心] 前沿比例: 取对抗方向上最远的前 20%

% 4. 时域聚类与平滑
TRAJ.time_cluster_k    = 5;     % [聚类] 时间窗口
TRAJ.traj_smooth_m     = 5;     % [平滑] 平滑窗口

%% ================= [Part 2] 核心计算流程 =================

% 1. 构建激活掩膜
num_samples = length(t_grid);
is_active_mask = false(num_samples, 1);
for k = 1:length(segments)
    is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
end

% 2. 预计算接收机位置缓存
rec_pos_cache = NaN(num_samples, 3);
for t = 1:num_samples
    if is_active_mask(t)
        [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
        try [rp, ~, ~] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); rec_pos_cache(t,:) = rp; catch, end
    end
end

% 3. 主循环: 稳定对抗追踪
K_cluster = TRAJ.time_cluster_k;
track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'active_count', {});
track_cnt = 0;
is_tracking_started = false; 

fprintf('--> 执行对抗追踪 (Cluster K=%d)...\n', K_cluster);

for t_start = 1 : K_cluster : num_samples
    
    t_end = min(t_start + K_cluster - 1, num_samples);
    range_indices = t_start : t_end;
    
    if ~any(is_active_mask(range_indices))
        is_tracking_started = false; 
        continue;
    end
    
    % --- 定义两组候选点 ---
    active_points  = []; % 波动组 (代表遮挡/手臂)
    stable_points  = []; % 稳定组 (代表背景/未遮挡)
    
    % --- A. 收集点云 ---
    for t = range_indices
        if ~is_active_mask(t), continue; end
        
        rec_pos = rec_pos_cache(t, :); 
        if any(isnan(rec_pos)), continue; end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
        
        % 获取 Step 1 传入的能量
        current_vols = volatility_matrix(t, :);
        
        [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
        try [~, ~, sat_states] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); catch, continue; end
        
        for s = 1:length(valid_sats)
            sid = valid_sats{s}; 
            if ~isfield(sat_states, sid), continue; end
            
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
            dist = norm([e, n, u]);
            vec_u = [e, n, u] / dist; 
            zen_deg = acosd(vec_u(3));
            
            if vec_u(3) <= 0, continue; end
            
            % [分类] 波动 vs 稳定
            is_active_sat = current_vols(s) > PARA.active_th; 
            
            % 投影
            t_int = TRAJ.gesture_height / vec_u(3); 
            pt_int = t_int * vec_u;
            
            if norm(pt_int(1:2)) > 5.0, continue; end
            
            % 收集逻辑 (应用安全仰角过滤)
            if is_active_sat
                % 波动组: 剔除过低仰角 (可能被身体挡住)
                if zen_deg >= ALG.zenith_safe_deg
                    active_points(end+1, :) = [pt_int(1), pt_int(2)];
                end
            else
                % 稳定组: 同样只参考有效工作区的卫星
                if zen_deg >= ALG.zenith_safe_deg
                    stable_points(end+1, :) = [pt_int(1), pt_int(2)];
                end
            end
        end
    end
    
    if isempty(active_points), continue; end
    
    % --- B. 对抗矢量计算 ---
    center_x = 0; center_y = 0;
    
    if size(stable_points, 1) > 2 && size(active_points, 1) > 2
        % 计算重心
        centroid_stable = mean(stable_points, 1);
        centroid_active = mean(active_points, 1);
        
        % 对抗矢量: 修正逻辑 (Arm occlusion -> Hand front)
        % 假设稳定区域在指尖前方，波动区域在手臂后方
        direction_vec = centroid_stable - centroid_active;
        
        if norm(direction_vec) > 1e-3
            direction_vec = direction_vec / norm(direction_vec);
        else
            direction_vec = [0, 0];
        end
        
        % 投影排序: 找在方向上投影最大的 Active 点 (最靠前的波动点)
        scores = active_points * direction_vec';
        [~, sort_idx] = sort(scores, 'descend'); 
        
        % 取 Top K%
        num_to_pick = max(1, ceil(size(active_points, 1) * ALG.top_k_percent));
        selected_pts = active_points(sort_idx(1:num_to_pick), :);
        
        center_x = mean(selected_pts(:,1));
        center_y = mean(selected_pts(:,2));
        
    else
        % 兜底: 退化为几何中心
        center_x = mean(active_points(:,1));
        center_y = mean(active_points(:,2));
    end
    
    % --- C. 动作触发锁 ---
    if ~is_tracking_started
        if norm([center_x, center_y]) > TRAJ.min_action_dist
            is_tracking_started = true; 
        else
            continue; 
        end
    end
    
    track_cnt = track_cnt + 1;
    track_results(track_cnt).t_idx = mean(range_indices);
    track_results(track_cnt).x = center_x;
    track_results(track_cnt).y = center_y;
    track_results(track_cnt).active_count = size(active_points, 1);
end

% 轨迹提取
if track_cnt > 0
    traj_x = [track_results.x]'; 
    traj_y = [track_results.y]'; 
    traj_t = [track_results.t_idx]'; 
    traj_active_counts = [track_results.active_count]';
    
    if TRAJ.traj_smooth_m > 1
        traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
        traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
    end
else
    traj_x = []; traj_y = []; traj_t = []; traj_active_counts = [];
end


%% ================= [Part 3] 绘图流程 =================
fprintf('\n--> 开始生成图表...\n');

if ~isempty(traj_x)
    figure('Name', 'Reconstructed Opposition Trajectory v5.1', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--'); 
    
    % 轨迹线 (用时间颜色映射)
    scatter(ax, traj_x, traj_y, 40, traj_t, 'filled', 'DisplayName', 'Finger Path');
    plot(ax, traj_x, traj_y, 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    c = colorbar; c.Label.String = 'Time Index'; colormap(ax, 'turbo');
    
    % 起终点
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'End');
    
    title({'Opposition Track v5.1', ...
           sprintf('Top %.0f%% Leading Edge | Cluster K=%d', ALG.top_k_percent*100, TRAJ.time_cluster_k)});
    legend('Location', 'best');
    
    max_range = max(max(abs(traj_x)), max(abs(traj_y)));
    if max_range < 0.5, max_range = 0.5; end
    xlim([-max_range*1.2, max_range*1.2]); ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end

fprintf('✅ v5.1 对抗分析完成 (已返回前沿轨迹数据)。\n');
end