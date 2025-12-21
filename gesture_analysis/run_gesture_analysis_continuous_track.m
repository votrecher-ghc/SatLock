% =========================================================================
% run_gesture_analysis_continuous_track_v2_3.m
% 功能: 鲁棒手势感知 - 强仰角加权版 (v2.3)
% 核心逻辑: 
%   1. [时域聚类]: 每 K 个采样点聚合计算一个重心。
%   2. [强仰角加权]: 使用 sin(el)^8 策略，极度压制非头顶区域（如手臂）的干扰。
%   3. [参数化过滤]: 单星波动阈值提取为可配置参数。
% =========================================================================

%% ================= [Part 1] 环境检查与参数设置 =================
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('错误: 请先加载 obs_data 到工作区!'); end
if ~exist('nav_data', 'var'), error('错误: 请先加载 nav_data 到工作区!'); end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动强仰角加权追踪 (v2.3: Power=8, Parametrized Volatility)...\n');

% --- [关键参数设置] ---
% 1. 信号预处理参数 (用于 Step 0)
PARA.diff_lag_N         = 5;     % 差分趋势窗口: 5个点
PARA.noise_cutoff_db    = 1;     % 噪声底噪阈值 (dB): 小于此波动视为静止
PARA.spike_th_db        = 2;     % 毛刺检测阈值 (dB)
PARA.spike_max_duration = 3;     % 毛刺最大持续点数: 超过则认为是真实信号

% 2. 能量检测与分段参数 (用于 Step 1)
PARA.smooth_window_sec = 1.5;    % 滑动平均窗口(秒): 用于计算长时间背景基线
PARA.gvi_threshold     = 2;      % GVI 激活阈值: 超过此值才开始计算轨迹
PARA.sampling_rate     = 25;     % 数据采样率 (Hz)
PARA.merge_gap_sec     = 0.2;    % 动作合并窗口: 断开小于0.2秒视为同一动作

% 3. 轨迹重建参数 (用于 Step 2)
TRAJ.gesture_height    = 0.20;   % [物理] 手势平面的假设高度 (米)
TRAJ.min_elevation     = 30;     % [物理] 最低仰角门限 (度): 低于此角度的卫星直接丢弃
TRAJ.min_action_dist   = 0.05;   % [触发] 动作死区 (米): 只有重心移出 5cm 圈外才开始记录

% [新增参数 v2.3]
PARA.min_sat_vol       = 0.5;    % [过滤] 单星波动门限 (dB): 只有波动超过此值的卫星才参与位置计算

% 4. 时域聚类参数 (V2 核心)
TRAJ.time_cluster_k    = 10;     % [聚类] 时间窗口: 每 10 个原始采样点(约0.4秒)合并为一个轨迹点
TRAJ.traj_smooth_m     = 5;      % [平滑] 最终轨迹的平滑窗口大小

% 5. [核心修改] 仰角加权策略 (Strategy: Zenith Focus)
% 权重公式: Weight = Energy * (sin(Elevation))^N
% N=1: 线性加权 (温和)
% N=2: 标准加权 (常用)
% N=4: 强力加权 (激进)
% N=8: 极强加权 (本次设定) -> 30度仰角权重仅为 0.39%，90度为 100%
TRAJ.elevation_power   = 8;      % [重要] 仰角权重指数。越高越只信头顶信号，抗手臂干扰越强。


%% ================= [Part 2] 核心计算流程 =================

% ----------------- [Step 0-1] 数据清洗与能量计算 (标准流程) -----------------
fprintf('--> [Step 0-1] 正在执行数据清洗与 GVI 计算...\n');

% 1. 提取有效卫星ID
all_sat_ids = {};
for i = 1:min(100, length(obs_data)), if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids);
valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

% 2. 构建时间网格
raw_times = [obs_data.time];
t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
num_samples = length(t_grid);
num_sats = length(valid_sats);
t_grid_plot = t_grid + hours(8) - seconds(18); % 北京时间用于绘图

% 3. 提取 CN0 数据 (插值对齐)
cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    % 智能查找可用频点
    for k = 1:min(50, length(obs_data))
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            fields = fieldnames(obs_data(k).data.(sat_id).snr);
            if ismember('S1C', fields), target_snr_code='S1C'; elseif ismember('S2I', fields), target_snr_code='S2I'; elseif ~isempty(fields), target_snr_code=fields{1}; end
            if ~isempty(target_snr_code), break; end
        end
    end
    if isempty(target_snr_code), continue; end
    
    % 提取与插值
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

% 4. 基准线清洗 (去噪)
for s = 1:num_sats
    raw_col = cn0_matrix(:, s); col = raw_col;
    valid_data = raw_col(~isnan(raw_col)); if isempty(valid_data), continue; end
    baseline = mode(round(valid_data)); 
    for t = 1:num_samples
        curr_val = raw_col(t); if isnan(curr_val), continue; end
        if abs(curr_val - baseline) > PARA.spike_th_db
             % 简单 Spike 抑制
             is_spike = true;
             for k=1:PARA.spike_max_duration, if t+k<=num_samples && abs(raw_col(t+k)-baseline)<=PARA.noise_cutoff_db, is_spike=true; break; end; end
        end
        if abs(curr_val - baseline) < PARA.noise_cutoff_db, col(t) = baseline; end
    end
    cn0_matrix(:, s) = col;
end

% 5. 计算波动能量 (GVI)
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth); 
gvi_curve_clean = movmean(sum(volatility_matrix, 2, 'omitnan'), 5);
is_active_mask = gvi_curve_clean > PARA.gvi_threshold; 

% ----------------- [Step 2] 强仰角加权聚类追踪 (v2.3 核心) -----------------
K_cluster = TRAJ.time_cluster_k;
fprintf('--> [Step 2] 执行聚类追踪 (N=%d 强仰角加权)...\n', TRAJ.elevation_power);

track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
track_cnt = 0;
is_tracking_started = false; 

% 预计算接收机位置缓存 (加速)
rec_pos_cache = NaN(num_samples, 3);
fprintf('    预计算接收机坐标...\n');
for t = 1:num_samples
    if is_active_mask(t), try [rp,~,~]=calculate_receiver_position(obs_data,nav_data,t); rec_pos_cache(t,:)=rp; catch, end; end
end

% === 主循环: 时间窗口步进 (每 K 个点处理一次) ===
for t_start = 1 : K_cluster : num_samples
    
    t_end = min(t_start + K_cluster - 1, num_samples);
    range_indices = t_start : t_end;
    
    % 若该窗口无任何激活信号，跳过
    if ~any(is_active_mask(range_indices))
        is_tracking_started = false; 
        continue;
    end
    
    % --- 1. 收集阶段: 收集窗口内所有点的(位置, 能量, 仰角) ---
    cluster_pts_x = [];
    cluster_pts_y = [];
    cluster_weights = [];
    
    for t = range_indices
        if ~is_active_mask(t), continue; end
        
        rec_pos = rec_pos_cache(t, :); 
        if any(isnan(rec_pos)), continue; end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
        
        % 获取当前时刻各卫星波动能量
        current_vols = volatility_matrix(t, :);
        
        % [修改点 v2.3] 使用 PARA.min_sat_vol 替代硬编码 0.1
        valid_energy_idx = find(current_vols > PARA.min_sat_vol); 
        
        if isempty(valid_energy_idx), continue; end
        
        % 计算卫星位置
        try [~, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, t); catch, continue; end
        
        for k = 1:length(valid_energy_idx)
            s_idx = valid_energy_idx(k); 
            sid = valid_sats{s_idx};
            if ~isfield(sat_states, sid), continue; end
            
            % 计算几何投影
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
            vec_u = [e, n, u] / norm([e, n, u]); 
            
            % 仰角检查
            zenith_rad = acos(vec_u(3)); % 天顶角(弧度)
            el_rad = pi/2 - zenith_rad;  % 仰角(弧度)
            el_deg = rad2deg(el_rad);
            
            if vec_u(3) <= 0 || el_deg < TRAJ.min_elevation, continue; end % 低于 15 度丢弃
            
            % 投影到手势平面
            t_int = TRAJ.gesture_height / vec_u(3); 
            pt = t_int * vec_u;
            if norm(pt(1:2)) > 5.0, continue; end % 丢弃异常远的投影
            
            % === [核心修改 v2.3] 强仰角加权计算 ===
            % 原始能量
            base_energy = current_vols(s_idx);
            
            % 仰角权重: (sin(el))^8
            % N=8: 30度时权重仅 0.39%，90度时权重 100%
            w_elevation = (sin(el_rad)) ^ TRAJ.elevation_power;
            
            % 最终合成权重 (只含仰角加权)
            final_w = base_energy * w_elevation;
            
            % 收集数据
            cluster_pts_x(end+1) = pt(1);
            cluster_pts_y(end+1) = pt(2);
            cluster_weights(end+1) = final_w;
        end
    end
    
    % --- 2. 聚合阶段: 计算加权质心 ---
    if isempty(cluster_weights), continue; end
    
    sum_w = sum(cluster_weights);
    if sum_w == 0, continue; end
    
    % 质心公式
    center_x = sum(cluster_pts_x .* cluster_weights) / sum_w;
    center_y = sum(cluster_pts_y .* cluster_weights) / sum_w;
    
    % --- 3. 动作触发锁逻辑 ---
    dist_from_origin = norm([center_x, center_y]);
    if ~is_tracking_started
        if dist_from_origin > TRAJ.min_action_dist
            is_tracking_started = true; 
        else
            continue; % 未触发，跳过
        end
    end
    
    % 记录结果
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
    
    % 轨迹平滑 (V2 标准步骤)
    if TRAJ.traj_smooth_m > 1
        fprintf('    执行轨迹平滑 (Window: %d)...\n', TRAJ.traj_smooth_m);
        traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
        traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
    end
else
    traj_x = []; traj_y = []; traj_t = [];
end


%% ================= [Part 3] 统一绘图流程 =================
fprintf('\n--> 开始生成图表...\n');

% [图表 1] GVI
figure('Name', 'GVI Overview', 'Position', [50, 100, 1000, 300], 'Color', 'w');
plot(t_grid_plot, gvi_curve_clean, 'b-', 'LineWidth', 1); hold on; % 默认蓝色
yline(PARA.gvi_threshold, 'k--', 'Activation Threshold');
title('GVI 能量波动图'); ylabel('GVI'); grid on; axis tight;

% [图表 2] 轨迹重建结果
if ~isempty(traj_x)
    figure('Name', 'Trajectory v2.3 (Extreme Elevation Weighted)', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    % 画接收机位置
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    % 画触发圈
    viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--');
    
    % 画轨迹 (默认蓝色)
    plot(ax, traj_x, traj_y, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Weighted Path');
    
    % 标记起终点
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'End');
    
    title({sprintf('continuous track (Power N=%d)', TRAJ.elevation_power), ...
           '使用参数化单星过滤'});
    legend('Location', 'best');
    
    % 自动缩放视角
    max_range = max([abs(traj_x); abs(traj_y); 0.5]);
    xlim([-max_range*1.2, max_range*1.2]); 
    ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end

fprintf('✅ v2.3 分析完成。\n');