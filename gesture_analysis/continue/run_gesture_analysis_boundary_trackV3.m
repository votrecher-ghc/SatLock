% =========================================================================
% run_gesture_analysis_continuous_track_v2_3.m
% 功能: 鲁棒手势感知 - 边界追踪修正版 (v3.2 - 逻辑修正 + 参数化过滤)
% 改进: 
%   1. [参数化]: 将单星波动阈值 (0.1) 提取为可配置参数 PARA.min_sat_vol。
%   2. [注释]: 严格保留 v3.1 原版详细注释，确保物理含义清晰。
%   3. [逻辑]: 保持 v3.1 的天顶角剔除逻辑和前沿追踪策略。
% =========================================================================

%% ================= [Part 1] 环境检查与参数设置 =================
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动前沿边界轨迹追踪 (v3.2: Detailed Comments & Parametrized)...\n');

% --- 参数设置 (完全恢复详细注释) ---

% [Step 0: 信号预处理 (Baseline Algorithm)]
% 作用: 利用众数锁定信号基准线，去除环境慢漂移和非手势的震荡
PARA.diff_lag_N         = 5;     % 趋势判定窗口 (点数): 只有连续N点偏离基准线才被视为有效趋势，否则视为震荡被抹平
PARA.noise_cutoff_db    = 1;     % 噪声容忍度 (dB): 信号波动小于等于此值时，被强制归零(视为平稳)
PARA.spike_th_db        = 2;     % 毛刺检测阈值 (dB): 瞬间跳变超过此值，将启动毛刺剔除逻辑
PARA.spike_max_duration = 3;     % 毛刺最大宽度 (点数): 如果跳变持续时间极短(<=3点)，则判定为观测噪声并剔除

% [Step 1: 能量门控与分段 (Mask Generation)]
% 作用: 生成时域掩膜，决定哪些时间段需要进行轨迹追踪
PARA.smooth_window_sec = 1;     % GVI基线窗口 (秒): 用于计算波动能量(GVI)的滑动平均窗口大小
PARA.gvi_threshold     = 2;     % GVI触发阈值 (dB): 全局波动能量超过此值，系统才进入"追踪模式"
PARA.sampling_rate     = 25;    % 数据采样率 (Hz): 必须与观测数据的实际采样率匹配
PARA.merge_gap_sec     = 0.4;   % 动作缝合阈值 (秒): 两个动作间隔小于0.2s时，视为同一个连续动作(防止画圆时中间断开)
PARA.min_duration_sec  = 0.2;   % 最小动作时长 (秒): 持续时间小于0.4s的波动被视为偶发干扰，直接丢弃

% [Step 2: 几何反演与轨迹追踪 (Geometry & Tracking)]
% 作用: 将一维的信号波动反演为二维的空间轨迹
TRAJ.gesture_height    = 0.30;  % 手势物理高度 (米): 假设手在接收机上方20cm处运动(用于视线投影)
TRAJ.min_elevation     = 10;    % 最低仰角门限 (度): 剔除地平线附近(仰角<15°)的多径严重卫星
TRAJ.min_action_dist   = 0.05;  % 动作触发死区 (米): 轨迹重心偏移超过5cm才开始记录，防止静止时的零漂

% [新增参数 v3.2]
PARA.min_sat_vol       = 1;   % [过滤] 单星波动门限 (dB): 只有波动超过此值的卫星才参与位置计算

% [关键抗干扰参数]
ALG.zenith_safe_deg    = 10;    % 安全仰角门限 (度): 用于剔除容易受身体遮挡的低仰角区域
ALG.top_k_percent      = 0.5;   % 前沿追踪比例 (0.2=20%): 仅利用距离圆心最远的20%的卫星点来计算重心(模拟指尖位置，抗手臂干扰)
TRAJ.time_cluster_k    = 5;     % 时域聚类窗口 (点数): 每5个采样点合并计算一个轨迹坐标(降低轨迹抖动)
TRAJ.traj_smooth_m     = 2;     % 轨迹平滑窗口 (点数): 对生成的最终轨迹进行移动平均平滑


%% ================= [Part 2] 核心计算流程 =================

% ----------------- [Step 0] 数据提取与基准线清洗 -----------------
fprintf('--> [Step 0] 数据提取与基准线清洗...\n');
% 1. 提取卫星与时间轴
all_sat_ids = {}; 
for i=1:min(100,length(obs_data))
    if ~isempty(obs_data(i).data), all_sat_ids=[all_sat_ids,fieldnames(obs_data(i).data)']; end
end
unique_sat_ids = unique(all_sat_ids); 
valid_sats = {};
for i=1:length(unique_sat_ids)
    sid=unique_sat_ids{i}; 
    if ismember(sid(1),['G','C','E','J']), valid_sats{end+1}=sid; end
end
raw_times = [obs_data.time];
t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
num_samples = length(t_grid); 
num_sats = length(valid_sats);
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
            if ~isnan(val) && val > 10
                s_times = [s_times; obs_data(k).time]; 
                s_cn0 = [s_cn0; val]; 
            end
        end
    end
    if length(s_times) > 20
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
    end
end

% 3. 全局基准线清洗 (Baseline Algorithm)
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

% ----------------- [Step 1 计算] 生成追踪掩膜 (Smart Mask) -----------------
fprintf('--> [Step 1] 生成智能追踪掩膜 (Mask Generation)...\n');
% SG 滤波
for s = 1:num_sats
    col = cn0_matrix(:, s); valid = ~isnan(col);
    if sum(valid) > 14
        idx = 1:length(col); filled = interp1(idx(valid), col(valid), idx, 'linear', 'extrap')';
        cn0_matrix(:, s) = sgolayfilt(filled, 2, 7); 
    end
end
% 计算波动能量 GVI
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth); 
gvi_curve_clean = movmean(sum(volatility_matrix, 2, 'omitnan'), 5);
% [核心优化] 智能掩膜生成 (替代旧的分段struct)
raw_mask = gvi_curve_clean > PARA.gvi_threshold;
is_active_mask = raw_mask;
% A. 缝合间隙 (Merge Gap)
gap_pts = round(PARA.merge_gap_sec * PARA.sampling_rate);
pad_mask = [0; raw_mask; 0];
edges = diff(pad_mask);
starts = find(edges == 1);
ends   = find(edges == -1) - 1;
for k = 1 : length(starts)-1
    gap_len = starts(k+1) - ends(k) - 1;
    if gap_len > 0 && gap_len <= gap_pts
        is_active_mask(ends(k)+1 : starts(k+1)-1) = 1; % 填补间隙
    end
end
% B. 剔除短时噪声 (Min Duration)
min_dur_pts = round(PARA.min_duration_sec * PARA.sampling_rate);
pad_mask = [0; is_active_mask; 0];
edges = diff(pad_mask);
starts = find(edges == 1);
ends   = find(edges == -1) - 1;
for k = 1 : length(starts)
    dur = ends(k) - starts(k) + 1;
    if dur < min_dur_pts
        is_active_mask(starts(k) : ends(k)) = 0; % 擦除短脉冲
    end
end
fprintf('   掩膜生成完毕。有效追踪时段占比: %.1f%%\n', sum(is_active_mask)/num_samples*100);

% ----------------- [Step 2 计算] 时空边界追踪 (Spatiotemporal Edge) -----------------
K_cluster = TRAJ.time_cluster_k;
fprintf('--> [Step 2] 执行边界追踪 (Zenith<%d°, Top%.0f%%)...\n', ALG.zenith_safe_deg, ALG.top_k_percent*100);
track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
track_cnt = 0;
is_tracking_started = false; 
% 预计算接收机位置缓存
rec_pos_cache = NaN(num_samples, 3);
for t = 1:num_samples
    if ~is_active_mask(t), continue; end
    try [rp, ~, ~] = calculate_receiver_position(obs_data, nav_data, t); rec_pos_cache(t,:) = rp; catch, end
end
% [时域聚类循环]
for t_start = 1 : K_cluster : num_samples
    
    t_end = min(t_start + K_cluster - 1, num_samples);
    range_indices = t_start : t_end;
    
    % 如果这一小段完全没激活，重置状态
    if ~any(is_active_mask(range_indices))
        is_tracking_started = false; 
        continue;
    end
    
    % 收集该时间窗内所有合格点
    window_candidates = struct('x', {}, 'y', {}, 'dist', {}, 'weight', {});
    wc_cnt = 0;
    sum_window_energy = 0;
    
    for t = range_indices
        if ~is_active_mask(t), continue; end
        
        rec_pos = rec_pos_cache(t, :);
        if any(isnan(rec_pos)), continue; end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
        
        current_vols = volatility_matrix(t, :);
        
        % [修改点 v3.2] 使用 PARA.min_sat_vol 替代硬编码 0.1
        valid_energy_idx = find(current_vols > PARA.min_sat_vol); 
        
        if isempty(valid_energy_idx), continue; end
        
        try [~, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, t); catch, continue; end
        
        for k = 1:length(valid_energy_idx)
            s_idx = valid_energy_idx(k); sid = valid_sats{s_idx};
            if ~isfield(sat_states, sid), continue; end
            
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
            dist = norm([e, n, u]);
            vec_u = [e, n, u] / dist; 
            zen_deg = acosd(vec_u(3)); % 计算天顶角 (0度=头顶)
            el_deg  = 90 - zen_deg;    % 计算仰角 (90度=头顶)
            
            if vec_u(3) <= 0, continue; end % 忽略地平线以下
            
            % --- [逻辑修复重点] 剔除身体遮挡 (通常在低仰角) ---
            % 新逻辑: 剔除仰角过低的卫星。如果仰角 < 25度，认为可能被身体挡住，跳过。
            if el_deg < ALG.zenith_safe_deg
                continue; 
            end
            
            % 计算投影点
            t_int = TRAJ.gesture_height / vec_u(3); 
            pt_int = t_int * vec_u;
            
            % --- [半径门控] Radius Gating ---
            % 物理意义: 限制有效感知的水平半径。
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
    
    % [边界提取] Top K% Furthest
    % 核心思想: 离圆心最远的点通常对应"指尖"位置 (前沿)
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
    track_results(track_cnt).total_energy = sum_window_energy;
end
if track_cnt > 0
    traj_x = [track_results.x]'; traj_y = [track_results.y]'; 
    traj_t = [track_results.t_idx]'; traj_e = [track_results.total_energy]';
    
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
figure('Name', 'Boundary Tracking v3.2: GVI Overview', 'Position', [50, 100, 1000, 400], 'Color', 'w');
plot(t_grid_plot, gvi_curve_clean, 'k-', 'LineWidth', 1); hold on;
yline(PARA.gvi_threshold, 'b--', 'Threshold');
% 画出最终的 mask 区域
area_x = [t_grid_plot(1); t_grid_plot; t_grid_plot(end)];
area_y = [0; double(is_active_mask) * max(gvi_curve_clean)*0.1; 0];
fill(area_x, area_y, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Tracking Active');
title('GVI 能量与激活区 (Masked)'); ylabel('GVI'); xlabel('Time'); grid on; axis tight;
% [图表 2] 轨迹图
if ~isempty(traj_x)
    figure('Name', 'Reconstructed Boundary Trajectory v3.2', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--'); 
    
    plot(ax, traj_x, traj_y, 'b-', 'LineWidth', 2, 'DisplayName', 'Finger Path');
    scatter(ax, traj_x, traj_y, 40, traj_e, 'filled', 'DisplayName', 'Cluster Point');
    c = colorbar; c.Label.String = 'Energy'; colormap(ax, 'parula');
    
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'End');
    
    title({'boundary trackV3', ...
           sprintf('Elev > %d° | Top %.0f%% | Cluster K=%d', ALG.zenith_safe_deg, ALG.top_k_percent*100, TRAJ.time_cluster_k)});
    legend('Location', 'best');
    
    max_range = max(max(abs(traj_x)), max(abs(traj_y)));
    if max_range < 0.5, max_range = 0.5; end
    xlim([-max_range*1.2, max_range*1.2]); ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end
fprintf('✅ 所有分析完成。\n');