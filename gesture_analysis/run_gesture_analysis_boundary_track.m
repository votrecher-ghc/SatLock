% =========================================================================
% run_gesture_analysis_boundary_track.m
% 功能: 鲁棒手势感知 - 稳定对抗边界追踪版 (v5.0)
% 核心逻辑: 
%   1. [预处理]: 基准线清洗 (Baseline Mode)。
%   2. [核心算法]: 稳定-波动对抗 (Stable-Active Opposition)。
%      利用"未波动卫星"的位置，推算出指尖指向，从而锁定波动区域的"最前沿边界"。
%      解决手臂遮挡导致重心后移的问题。
%   3. [时空优化]: 引入时域聚类 (Cluster K) 和轨迹平滑 (Smooth M)。
% =========================================================================

%% ================= [Part 1] 环境检查与参数设置 =================
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动稳定对抗边界追踪 (v5.0: Stable-Active Opposition)...\n');

% --- 参数设置 ---
% [Step 0: 信号预处理 (Baseline)]
PARA.diff_lag_N         = 5;     % 趋势窗口
PARA.noise_cutoff_db    = 1;     % 噪声阈值
PARA.spike_th_db        = 2;     % 毛刺阈值
PARA.spike_max_duration = 3;     % 毛刺点数

% [Step 1: 能量门控]
PARA.smooth_window_sec = 1.5;   % GVI平滑窗口
PARA.gvi_threshold     = 5;     % GVI能量阈值 (判定是否有动作)
PARA.sampling_rate     = 25;    
PARA.merge_gap_sec     = 0.2;   
PARA.min_duration_sec  = 0.4;   

% [Step 2: 几何与追踪参数]
TRAJ.gesture_height    = 0.20;  % 手势高度
TRAJ.min_elevation     = 15;    % 最低仰角
TRAJ.min_action_dist   = 0.05;  % 动作触发距离 (>5cm)
ALG.zenith_safe_deg    = 30;    % 天顶角掩膜 (剔除身体核心区遮挡)

% [Step 2: 边界锁定参数]
ALG.top_k_percent      = 0.2;   % 前沿比例: 取对抗方向上最远的前 20% 点

% [Step 2: 聚类与平滑 (用户指定)]
TRAJ.time_cluster_k    = 5;    % 时域聚类: 每 10 个采样点聚类出一个坐标
TRAJ.traj_smooth_m     = 5;     % 轨迹平滑: 对最终轨迹做 5 点平滑

%% ================= [Part 2] 核心计算流程 =================

% ----------------- [Step 0] 数据提取与基准线清洗 -----------------
fprintf('--> [Step 0] 数据提取与基准线清洗...\n');

% 1. 提取卫星
all_sat_ids = {}; for i=1:min(100,length(obs_data)), if~isempty(obs_data(i).data), all_sat_ids=[all_sat_ids,fieldnames(obs_data(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids); valid_sats = {};
for i=1:length(unique_sat_ids), sid=unique_sat_ids{i}; if ismember(sid(1),['G','C','E','J']), valid_sats{end+1}=sid; end; end

raw_times = [obs_data.time];
t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
num_samples = length(t_grid); num_sats = length(valid_sats);
t_grid_plot = t_grid + hours(8) - seconds(18); 

% 2. 提取原始数据
cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx}; target_snr_code = '';
    for k = 1:min(50, length(obs_data))
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            fields = fieldnames(obs_data(k).data.(sat_id).snr);
            if ~isempty(fields), target_snr_code = fields{1}; break; end
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

% 3. 基准线清洗 (Baseline Algorithm)
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


% ----------------- [Step 1 计算] 生成掩膜 -----------------
fprintf('--> [Step 1] 生成智能追踪掩膜...\n');

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

% 生成二值掩膜
raw_mask = gvi_curve_clean > PARA.gvi_threshold;
is_active_mask = raw_mask;

% 缝合间隙
gap_pts = round(PARA.merge_gap_sec * PARA.sampling_rate);
pad_mask = [0; raw_mask; 0]; edges = diff(pad_mask);
starts = find(edges == 1); ends = find(edges == -1) - 1;
for k = 1 : length(starts)-1
    gap_len = starts(k+1) - ends(k) - 1;
    if gap_len > 0 && gap_len <= gap_pts, is_active_mask(ends(k)+1 : starts(k+1)-1) = 1; end
end

% 剔除短时噪声
min_dur_pts = round(PARA.min_duration_sec * PARA.sampling_rate);
pad_mask = [0; is_active_mask; 0]; edges = diff(pad_mask);
starts = find(edges == 1); ends = find(edges == -1) - 1;
for k = 1 : length(starts)
    if (ends(k) - starts(k) + 1) < min_dur_pts, is_active_mask(starts(k) : ends(k)) = 0; end
end


% ----------------- [Step 2 计算] 稳定对抗边界追踪 -----------------
K_cluster = TRAJ.time_cluster_k;
fprintf('--> [Step 2] 执行稳定对抗追踪 (Stable-Active Opposition, Cluster K=%d)...\n', K_cluster);

track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'active_count', {});
track_cnt = 0;
is_tracking_started = false; 

rec_pos_cache = NaN(num_samples, 3);
for t = 1:num_samples
    if ~is_active_mask(t), continue; end
    try [rp, ~, ~] = calculate_receiver_position(obs_data, nav_data, t); rec_pos_cache(t,:) = rp; catch, end
end

% 时域聚类循环
for t_start = 1 : K_cluster : num_samples
    
    t_end = min(t_start + K_cluster - 1, num_samples);
    range_indices = t_start : t_end;
    
    if ~any(is_active_mask(range_indices))
        is_tracking_started = false; 
        continue;
    end
    
    % 定义两组候选点
    active_points  = []; % 波动组 (代表遮挡)
    stable_points  = []; % 稳定组 (代表未遮挡背景)
    
    % 在时间窗内收集所有点
    for t = range_indices
        if ~is_active_mask(t), continue; end
        rec_pos = rec_pos_cache(t, :); if any(isnan(rec_pos)), continue; end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
        
        current_vols = volatility_matrix(t, :);
        
        try [~, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, t); catch, continue; end
        
        for s = 1:num_sats
            sid = valid_sats{s}; 
            if ~isfield(sat_states, sid), continue; end
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
            vec_u = [e, n, u] / norm([e, n, u]); zen_deg = acosd(vec_u(3));
            
            if vec_u(3) <= 0, continue; end
            
            % [分类] 波动 vs 稳定
            % 波动判据: GVI > 阈值 (这里的阈值可以用全局的，也可以用局部自适应，暂用0.5倍GVI阈值作为单星判据)
            % 注意：GVI是总能量，单星波动阈值一般较小，这里取 1.0 (baseline清洗后的波动)
            is_active_sat = current_vols(s) > 1.0; 
            
            % 投影计算
            t_int = TRAJ.gesture_height / vec_u(3); pt_int = t_int * vec_u;
            dist_from_center = norm(pt_int(1:2));
            if dist_from_center > 5.0, continue; end
            
            if is_active_sat
                % 波动组: 需要剔除身体遮挡 (天顶角过小)
                if zen_deg >= ALG.zenith_safe_deg
                    active_points(end+1, :) = [pt_int(1), pt_int(2)];
                end
            else
                % 稳定组: 同样剔除身体遮挡区，只看有效工作区的稳定卫星
                if zen_deg >= ALG.zenith_safe_deg
                    stable_points(end+1, :) = [pt_int(1), pt_int(2)];
                end
            end
        end
    end
    
    if isempty(active_points), continue; end
    
    % [核心算法] 对抗矢量计算
    % 如果有足够多的稳定点，计算稳定重心 -> 波动重心 的反向矢量
    center_x = 0; center_y = 0;
    
    if size(stable_points, 1) > 2 && size(active_points, 1) > 2
        % 计算重心
        centroid_stable = mean(stable_points, 1);
        centroid_active = mean(active_points, 1);
        
        % 对抗矢量: 从 波动中心 指向 稳定中心 (即指向未遮挡区域)
        % 修正: 手指在前，手臂在后。手臂遮挡区(Active) -> 手指前方未遮挡区(Stable)
        direction_vec = centroid_stable - centroid_active;
        if norm(direction_vec) > 1e-3
            direction_vec = direction_vec / norm(direction_vec);
        else
            direction_vec = [0, 0];
        end
        
        % 投影排序: 找在 direction_vec 方向上投影最大的 Active 点
        % Proj = P dot Dir
        scores = active_points * direction_vec';
        [~, sort_idx] = sort(scores, 'descend'); % 分数越高，越靠近稳定区边界
        
        % 取前 K% (前沿)
        num_to_pick = max(1, ceil(size(active_points, 1) * ALG.top_k_percent));
        selected_pts = active_points(sort_idx(1:num_to_pick), :);
        
        center_x = mean(selected_pts(:,1));
        center_y = mean(selected_pts(:,2));
        
    else
        % [兜底逻辑] 如果没有稳定点(全遮挡)或波动点太少，退化为波动点几何中心
        % 这里的几何中心是指 Convex Hull 的中心，或者简单平均
        % 为了稳健，使用简单平均
        center_x = mean(active_points(:,1));
        center_y = mean(active_points(:,2));
    end
    
    % [动作触发锁]
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

if track_cnt > 0
    traj_x = [track_results.x]'; traj_y = [track_results.y]'; 
    traj_t = [track_results.t_idx]'; 
    
    % [轨迹平滑]
    if TRAJ.traj_smooth_m > 1
        traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
        traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
    end
else
    traj_x = []; traj_y = []; traj_t = [];
end


%% ================= [Part 3] 统一绘图流程 =================
fprintf('\n--> 开始生成图表...\n');

% [图表 1] GVI
figure('Name', 'Opposition Tracking: GVI Overview', 'Position', [50, 100, 1000, 400], 'Color', 'w');
plot(t_grid_plot, gvi_curve_clean, 'k-', 'LineWidth', 1); hold on;
yline(PARA.gvi_threshold, 'b--', 'Threshold');
if ~isempty(traj_t)
    area_x = [t_grid_plot(1); t_grid_plot; t_grid_plot(end)];
    area_y = [0; double(is_active_mask) * max(gvi_curve_clean)*0.1; 0];
    fill(area_x, area_y, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Active Zone');
end
title('GVI 能量与激活区'); ylabel('GVI'); xlabel('Time'); grid on; axis tight;

% [图表 2] 轨迹图
if ~isempty(traj_x)
    figure('Name', 'Reconstructed Opposition Trajectory', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--'); 
    
    % 轨迹线 (用时间颜色映射)
    scatter(ax, traj_x, traj_y, 40, traj_t, 'filled', 'DisplayName', 'Finger Path');
    plot(ax, traj_x, traj_y, 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    c = colorbar; c.Label.String = 'Time Index'; colormap(ax, 'turbo');
    
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'End');
    
    title({'稳定对抗前沿轨迹 (Stable-Active Opposition)', ...
           sprintf('Top %.0f%% Leading Edge | Cluster K=%d', ALG.top_k_percent*100, TRAJ.time_cluster_k)});
    legend('Location', 'best');
    
    max_range = max(max(abs(traj_x)), max(abs(traj_y)));
    if max_range < 0.5, max_range = 0.5; end
    xlim([-max_range*1.2, max_range*1.2]); ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end

fprintf('✅ 所有分析完成。\n');










