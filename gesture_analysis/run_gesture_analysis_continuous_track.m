% =========================================================================
% run_gesture_analysis_continuous_track_v2.m
% 功能: 鲁棒手势感知 - 连续轨迹追踪版 (v2.0 - 时域聚类增强)
% 核心改进: 
%   1. [时域聚类]: 每 K 个采样点聚合计算一个重心，降低采样率，提升稳定性。
%   2. [轨迹平滑]: 对聚合后的轨迹点进行 M 点滑动平均。
% =========================================================================

%% ================= [Part 1] 环境检查与参数设置 =================
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动连续手势轨迹追踪 (v2.0: Temporal Clustering)...\n');

% --- 参数设置 ---
% [Step 0: 信号预处理参数 (Baseline Algorithm)]
PARA.diff_lag_N         = 10;     % 趋势窗口
PARA.noise_cutoff_db    = 1;     % 噪声阈值
PARA.spike_th_db        = 2;     % 毛刺阈值
PARA.spike_max_duration = 3;     % 毛刺宽度

% [Step 1: 分段与检测参数]
PARA.smooth_window_sec = 1.5;   
PARA.gvi_threshold     = 5;     
PARA.sampling_rate     = 25;    
PARA.merge_gap_sec     = 0.01;  
PARA.min_duration_sec  = 0.4;   

% [Step 2: 3D 轨迹投影参数]
TRAJ.gesture_height    = 0.20;  % 手势物理高度
TRAJ.min_elevation     = 15;    % 最低仰角
TRAJ.energy_th_ratio   = 0.2;   % 能量筛选比例
TRAJ.min_action_dist   = 0.05;  % 动作触发距离 (5cm)

% [新增: 聚类与平滑参数]
TRAJ.time_cluster_k    = 10;     % 1. 聚类窗口: 每 5 个原始点 -> 聚合成 1 个轨迹点
TRAJ.traj_smooth_m     = 5;     % 2. 轨迹平滑: 对聚合后的点进行 3 点滑动平均

%% ================= [Part 2] 核心计算流程 =================

% ----------------- [Step 0 计算] 数据提取与基准线清洗 -----------------
fprintf('--> [Step 0] 提取数据并执行基准线清洗...\n');

% 1. 提取卫星ID
all_sat_ids = {};
for i = 1:min(100, length(obs_data)), if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids);
valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

% 2. 构建时间轴
raw_times = [obs_data.time];
t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
num_samples = length(t_grid);
num_sats = length(valid_sats);
t_grid_plot = t_grid + hours(8) - seconds(18); 

% 3. 提取原始数据
cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    for k = 1:min(50, length(obs_data))
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            fields = fieldnames(obs_data(k).data.(sat_id).snr);
            if ismember('S1C', fields), target_snr_code = 'S1C'; elseif ismember('S2I', fields), target_snr_code = 'S2I'; else, if ~isempty(fields), target_snr_code = fields{1}; end; end
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

% 4. 全局基准线清洗
N = PARA.diff_lag_N; NoiseTh = PARA.noise_cutoff_db; SpikeTh = PARA.spike_th_db; SpikeDur = PARA.spike_max_duration;
for s = 1:num_sats
    raw_col = cn0_matrix(:, s); col = raw_col;
    valid_data = raw_col(~isnan(raw_col)); if isempty(valid_data), continue; end
    baseline = mode(round(valid_data)); 
    for t = 1:num_samples
        curr_val = raw_col(t); if isnan(curr_val), continue; end
        if abs(curr_val - baseline) > SpikeTh
            is_spike = false;
            for k = 1:SpikeDur
                if t + k > num_samples, break; end
                if abs(raw_col(t+k) - baseline) <= NoiseTh, is_spike = true; break; end
            end
            if is_spike, col(t) = baseline; continue; end
        end
        win_end = min(t + N - 1, num_samples);
        diffs = raw_col(t : win_end) - baseline;
        sig_diffs = diffs(abs(diffs) > NoiseTh);
        if isempty(sig_diffs)
            col(t) = baseline;
        else
            if all(sig_diffs > 0) || all(sig_diffs < 0), col(t) = curr_val; else, col(t) = baseline; end
        end
    end
    cn0_matrix(:, s) = col;
end


% ----------------- [Step 1 计算] SG滤波与波动矩阵计算 -----------------
fprintf('--> [Step 1] SG滤波与波动能量计算...\n');

for s = 1:num_sats
    col = cn0_matrix(:, s); valid = ~isnan(col);
    if sum(valid) > 14
        idx = 1:length(col); filled = interp1(idx(valid), col(valid), idx, 'linear', 'extrap')';
        cn0_matrix(:, s) = sgolayfilt(filled, 2, 7); 
    end
end

cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth); 
gvi_curve_clean = movmean(sum(volatility_matrix, 2, 'omitnan'), 5);
is_active_mask = gvi_curve_clean > PARA.gvi_threshold; 


% ----------------- [Step 2 计算] 时域聚类轨迹追踪 -----------------
K_cluster = TRAJ.time_cluster_k;
fprintf('--> [Step 2] 时域聚类追踪 (Cluster Window: %d pts)...\n', K_cluster);

track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {}, 'active_sats', {});
track_cnt = 0;
is_tracking_started = false; 

% 预先计算所有时刻的接收机位置 (缓存以加速)
fprintf('    预计算接收机位置缓存...\n');
rec_pos_cache = NaN(num_samples, 3);
for t = 1:num_samples
    if ~is_active_mask(t), continue; end
    try [rp, ~, ~] = calculate_receiver_position(obs_data, nav_data, t); rec_pos_cache(t,:) = rp; catch, end
end

% [修改] 循环步长改为 K_cluster，进行降采样聚类
for t_start = 1 : K_cluster : num_samples
    
    t_end = min(t_start + K_cluster - 1, num_samples);
    range_indices = t_start : t_end;
    
    % 1. 检查该窗口内是否有任意一刻激活
    if ~any(is_active_mask(range_indices))
        is_tracking_started = false; % 能量断裂，重置
        continue;
    end
    
    % 2. 收集窗口内【所有时刻、所有卫星】的加权点
    cluster_points_x = [];
    cluster_points_y = [];
    cluster_weights  = [];
    
    for t = range_indices
        if ~is_active_mask(t), continue; end
        
        % 读取缓存的接收机位置
        rec_pos = rec_pos_cache(t, :);
        if any(isnan(rec_pos)), continue; end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
        
        % 提取波动能量
        current_vols = volatility_matrix(t, :);
        valid_energy_idx = find(current_vols > 0.1); 
        if isempty(valid_energy_idx), continue; end
        
        % 获取卫星状态
        try [~, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, t); catch, continue; end
        
        for k = 1:length(valid_energy_idx)
            s_idx = valid_energy_idx(k); sid = valid_sats{s_idx};
            if ~isfield(sat_states, sid), continue; end
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
            vec_u = [e, n, u] / norm([e, n, u]); zen_deg = acosd(vec_u(3));
            
            if vec_u(3) <= 0 || (90-zen_deg) < TRAJ.min_elevation, continue; end
            t_int = TRAJ.gesture_height / vec_u(3); pt_int = t_int * vec_u;
            if norm(pt_int(1:2)) > 5.0, continue; end
            
            % 累积点
            w = current_vols(s_idx); 
            cluster_points_x(end+1) = pt_int(1); 
            cluster_points_y(end+1) = pt_int(2); 
            cluster_weights(end+1)  = w;
        end
    end
    
    % 3. 计算该时间窗口的总重心
    if isempty(cluster_weights), continue; end
    sum_w = sum(cluster_weights); if sum_w == 0, continue; end
    
    center_x = sum(cluster_points_x .* cluster_weights) / sum_w;
    center_y = sum(cluster_points_y .* cluster_weights) / sum_w;
    
    % 4. 动作触发锁逻辑
    dist_from_origin = norm([center_x, center_y]);
    if ~is_tracking_started
        if dist_from_origin > TRAJ.min_action_dist
            is_tracking_started = true; 
        else
            continue; 
        end
    end
    
    % 5. 记录聚类后的点 (时间戳取窗口中心)
    track_cnt = track_cnt + 1;
    track_results(track_cnt).t_idx = mean(range_indices);
    track_results(track_cnt).x = center_x;
    track_results(track_cnt).y = center_y;
    track_results(track_cnt).total_energy = sum_w;
end

% 轨迹提取与后处理
if track_cnt > 0
    traj_x = [track_results.x]'; 
    traj_y = [track_results.y]'; 
    traj_t = [track_results.t_idx]'; 
    traj_e = [track_results.total_energy]';
    
    % [新增] 对聚合后的轨迹进行平滑
    if TRAJ.traj_smooth_m > 1
        fprintf('    执行轨迹平滑 (Window: %d pts)...\n', TRAJ.traj_smooth_m);
        traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
        traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
    end
else
    traj_x = []; traj_y = []; traj_t = [];
end


%% ================= [Part 3] 统一绘图流程 =================
fprintf('\n--> 开始生成图表...\n');

% [图表 1] GVI
figure('Name', 'Continuous Tracking: GVI Overview', 'Position', [50, 100, 1000, 400], 'Color', 'w');
plot(t_grid_plot, gvi_curve_clean, 'k-', 'LineWidth', 1); hold on;
yline(PARA.gvi_threshold, 'b--', 'Activation Threshold');
if ~isempty(traj_t)
    area_x = [t_grid_plot(1); t_grid_plot; t_grid_plot(end)];
    area_y = [0; double(is_active_mask) * max(gvi_curve_clean)*0.1; 0];
    fill(area_x, area_y, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'GVI Active Zone');
end
title('GVI 能量与激活区'); ylabel('GVI'); xlabel('Time'); grid on; axis tight;

% [图表 2] 轨迹图
if ~isempty(traj_x)
    figure('Name', 'Reconstructed Continuous Trajectory', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    % 画接收机和触发圈
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--'); 
    
    % 画轨迹 (点更稀疏，线更平滑)
    plot(ax, traj_x, traj_y, 'b-', 'LineWidth', 2, 'DisplayName', 'Aggregated Path');
    scatter(ax, traj_x, traj_y, 40, traj_e, 'filled', 'DisplayName', 'Cluster Center');
    c = colorbar; c.Label.String = 'Cluster Energy'; colormap(ax, 'parula');
    
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'End');
    
    title({sprintf('时域聚类轨迹 (Cluster K=%d, Smooth M=%d)', TRAJ.time_cluster_k, TRAJ.traj_smooth_m), ...
           sprintf('Trigger Dist: %.2fm', TRAJ.min_action_dist)});
    legend('Location', 'best');
    
    max_range = max(max(abs(traj_x)), max(abs(traj_y)));
    if max_range < 0.5, max_range = 0.5; end
    xlim([-max_range*1.2, max_range*1.2]); ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end

fprintf('✅ 所有分析完成。\n');