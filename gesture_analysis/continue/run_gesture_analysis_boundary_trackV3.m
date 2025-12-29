% =========================================================================
% run_gesture_analysis_boundary_trackV3.m (函数版)
% 功能: 鲁棒手势感知 Step 2 - 边界/前沿追踪版 (v3.3 Struct Input)
% 描述:
%   该函数是 Step 2 的一种实现，针对"远端手势" (如伸出手画图) 设计。
%   它接收 Step 1 处理后的干净数据，不再计算所有遮挡点的重心，而是筛选
%   距离身体中心"最远"的一批点 (Top-K%)，以此模拟追踪"指尖"的位置，
%   有效避免了手掌或手臂遮挡信号将轨迹重心向身体侧拉扯的问题。
%
% [调用格式]:
%   [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_boundary_trackV3(obs_clean, nav_data, step1_res);
%
% [输入参数]:
%   1. obs_clean (struct): Step 1 返回的清洗后观测数据 (用于定位)。
%   2. nav_data (struct): 导航星历数据。
%   3. step1_res (struct): Step 1 返回的结果包 (含能量矩阵、分段信息等)。
%
% [返回值说明]:
%   1. traj_x / traj_y: [指尖/前沿] 的轨迹坐标 (米)。
%   2. traj_t: 轨迹点时间索引。
%   3. traj_e: 前沿聚类点的总能量。
%
% [核心算法]:
%   1. Top-K Filtering: 仅保留距离圆心最远的 K% (如50%) 的点参与计算。
%   2. Zenith Safe Mask: 强制剔除低仰角区域，防止身体核心区遮挡干扰。
% =========================================================================

function [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_boundary_trackV3(obs_clean, nav_data, step1_res)

% --- 1. 数据解包 (Unpack Data) ---
segments = step1_res.segments;
volatility_matrix = step1_res.volatility_matrix;
t_grid = step1_res.t_grid;
valid_sats = step1_res.valid_sats;
PARA = step1_res.PARA; % 获取 Step 1 的参数

if isempty(segments)
    fprintf('⚠️ 警告: step1_res 中无有效分段，Step 2 (Boundary) 终止。\n');
    traj_x=[]; traj_y=[]; traj_t=[]; traj_e=[];
    return;
end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> [Step 2] 启动前沿边界追踪 (Function版 v3.3: Boundary/Top-K)...\n');

%% ================= [Part 1] 参数设置 =================

% 1. 轨迹与几何参数
TRAJ.gesture_height    = 0.30;  % [物理] 手势平面高度 (米)
TRAJ.min_elevation     = 10;    % [物理] 最低仰角门限 (度)
TRAJ.min_action_dist   = 0.05;  % [触发] 动作死区 (米)

% 2. 过滤参数
% [修复] 使用 Step 1 传入的 PARA 进行兼容，或直接定义局部参数
PARA.min_sat_vol       = 1.0;   % [过滤] 单星波动门限 (dB)

% 3. 抗干扰与边界算法参数
ALG.zenith_safe_deg    = 10;    % [掩膜] 安全仰角门限 (度): 剔除身体遮挡
ALG.top_k_percent      = 0.5;   % [核心] 前沿追踪比例: 仅取最远的 50% 点

% 4. 时域聚类
TRAJ.time_cluster_k    = 5;     % [聚类] 时间窗口
TRAJ.traj_smooth_m     = 2;     % [平滑] 轨迹平滑窗口

%% ================= [Part 2] 核心计算流程 =================

% 1. 构建激活掩膜 (Mask)
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

% 3. 主循环: 边界追踪
K_cluster = TRAJ.time_cluster_k;
fprintf('--> 执行边界追踪 (Zenith<%d°, Top%.0f%%)...\n', ALG.zenith_safe_deg, ALG.top_k_percent*100);

track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
track_cnt = 0;
is_tracking_started = false; 

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
        
        % 读取能量矩阵
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
            dist = norm([e, n, u]);
            vec_u = [e, n, u] / dist; 
            
            zen_deg = acosd(vec_u(3)); 
            el_deg  = 90 - zen_deg;
            
            if vec_u(3) <= 0, continue; end 
            
            % --- [安全掩膜] 剔除身体遮挡 ---
            if el_deg < ALG.zenith_safe_deg
                continue; 
            end
            
            % 计算投影点
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
    
    % --- B. 边界提取 (Top-K Logic) ---
    % 核心思想: 离圆心最远的点通常对应"指尖"位置
    all_dists = [window_candidates.dist];
    [~, sort_idx] = sort(all_dists, 'descend'); 
    
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
    track_results(track_cnt).t_idx = mean(range_indices);
    track_results(track_cnt).x = center_x;
    track_results(track_cnt).y = center_y;
    track_results(track_cnt).total_energy = sum_window_energy;
end

% 轨迹提取
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

%% ================= [Part 3] 绘图流程 =================
fprintf('\n--> 开始生成图表...\n');

if ~isempty(traj_x)
    figure('Name', 'Reconstructed Boundary Trajectory v3.3', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    % 画接收机
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--'); 
    
    % 画轨迹 (带能量颜色)
    plot(ax, traj_x, traj_y, 'b-', 'LineWidth', 2, 'DisplayName', 'Finger Path');
    scatter(ax, traj_x, traj_y, 40, traj_e, 'filled', 'DisplayName', 'Cluster Point');
    c = colorbar; c.Label.String = 'Cluster Energy'; colormap(ax, 'parula');
    
    % 起终点
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'End');
    
    title({'Boundary Track v3.3', ...
           sprintf('Elev > %d° | Top %.0f%% | Cluster K=%d', ALG.zenith_safe_deg, ALG.top_k_percent*100, TRAJ.time_cluster_k)});
    legend('Location', 'best');
    
    max_range = max(max(abs(traj_x)), max(abs(traj_y)));
    if max_range < 0.5, max_range = 0.5; end
    xlim([-max_range*1.2, max_range*1.2]); ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end

fprintf('✅ v3.3 边界追踪分析完成 (已返回轨迹数据)。\n');
end