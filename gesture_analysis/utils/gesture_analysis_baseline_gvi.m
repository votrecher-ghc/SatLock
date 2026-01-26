% =========================================================================
% gesture_analysis_baseline_gvi.m (函数版)
% 功能: 手势分析 Step 1 - 基准线清洗、GVI特征提取与动作分段 (v7.4 GLONASS Support)
% 描述:
%   [重要更新] 更改了 Volatility (波动能量) 的计算方式。
%   不再使用 "当前值 - 滑动平均"，改为 "当前值 - 全局基准线" (abs(Data - Baseline))。
%
%   [v7.4 Update]:
%   1. [修复] 增加了对 'R' (GLONASS) 卫星的支持。之前的版本只包含了 G/C/E/J。
%   2. 包含 Data Injection 模块：将清洗后的基准线数据回填至 obs_clean。
%   3. 修复了 containers.Map 报错问题。
%
% [调用格式]:
%   [obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data);
%
% [输入参数]:
%   1. obs_data (struct数组): 
%      原始观测数据。
% =========================================================================

function [obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data)

fprintf('--> [Step 1] 启动基准线清洗与GVI分段 (Deviation Mode v7.4 + GLONASS)...\n');

%% 1. 初始化与参数设置
obs_clean = obs_data; % 初始化副本

% --- 核心参数 ---
PARA.sampling_rate      = 25;    % [系统] 锁定 25Hz
% [Baseline] 清洗参数
PARA.diff_lag_N         = 5;     % 趋势窗口
PARA.noise_cutoff_db    = 1.0;   % 噪声阈值
PARA.spike_th_db        = 1.0;   % 毛刺阈值
PARA.spike_max_duration = 5;     % 毛刺最大持续点数
% [GVI & Segment] 参数
PARA.smooth_window_pts  = 25;    % [平滑窗口]
PARA.gvi_threshold      = 20;     % [激活阈值]
PARA.merge_gap_pts      = 10;    % [合并容差]
PARA.min_duration_pts   = 3;     % [最小力度]

fprintf('    算法模式: 绝对偏差 (Abs Deviation from Baseline)\n');

%% 2. 数据对齐 (插值) & 矩阵构建
% 2.1 提取卫星
all_sat_ids = {};
scan_range = unique([1:min(100, length(obs_data)), max(1, length(obs_data)-100):length(obs_data)]);
for i = scan_range
    if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end
end
unique_sat_ids = unique(all_sat_ids);
valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    % [修复] 在此处添加了 'R' 以支持 GLONASS
    if ismember(sid(1), ['G','R','C','E','J']) 
        valid_sats{end+1} = sid; 
    end
end

% 2.2 构建网格
raw_times = [obs_data.time];
t_grid  = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
num_samples = length(t_grid);
num_sats = length(valid_sats);

% 2.3 插值原始 SNR (构建 cn0_matrix)
cn0_matrix = NaN(num_samples, num_sats);

for s_idx = 1:num_sats
    sid = valid_sats{s_idx};
    target_code = '';
    % 寻找该卫星可用的 SNR 信号码
    for k = 1:min(50, length(obs_data)) 
        if isfield(obs_data(k).data, sid) && isfield(obs_data(k).data.(sid), 'snr')
            fds = fieldnames(obs_data(k).data.(sid).snr);
            if ~isempty(fds), target_code = fds{1}; break; end
        end
    end
    if isempty(target_code), continue; end
    
    % 提取离散数据
    s_times = []; s_vals = [];
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sid) && isfield(obs_data(k).data.(sid).snr, target_code)
            val = obs_data(k).data.(sid).snr.(target_code);
            if ~isnan(val) && val > 10
                s_times = [s_times; obs_data(k).time]; s_vals = [s_vals; val];
            end
        end
    end
    
    % 执行线性插值
    if length(s_times) > 5
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_vals(u_idx), t_grid, 'linear', NaN);
    end
end

%% 3. 基准线清洗 (清洗算法逻辑)
fprintf('--> [预处理] 执行基准线清洗...\n');
N = PARA.diff_lag_N;
NoiseTh = PARA.noise_cutoff_db;
SpikeTh = PARA.spike_th_db;
SpikeDur = PARA.spike_max_duration;

sat_baselines = NaN(1, num_sats);

for s = 1:num_sats
    raw_col = cn0_matrix(:, s);
    valid_mask = ~isnan(raw_col);
    if sum(valid_mask) < 10, continue; end
    
    % 计算全局基准线 (众数)
    base_val = mode(round(raw_col(valid_mask))); 
    sat_baselines(s) = base_val; 
    
    col = raw_col;
    
    % 清洗循环
    for t = 1:num_samples
        curr_val = raw_col(t);
        if isnan(curr_val), continue; end
        
        % A. Spike Check
        if abs(curr_val - base_val) > SpikeTh
            is_spike = false;
            for k = 1:SpikeDur
                if t + k > num_samples, break; end
                if abs(raw_col(t+k) - base_val) <= NoiseTh
                    is_spike = true; break;
                end
            end
            if is_spike, col(t) = base_val; continue; end
        end
        
        % B. Trend Check
        win_end = min(t + N - 1, num_samples);
        diffs = raw_col(t : win_end) - base_val;
        sig_diffs = diffs(abs(diffs) > NoiseTh);
        
        if isempty(sig_diffs)
            col(t) = base_val; 
        else
            if all(sig_diffs > 0) || all(sig_diffs < 0)
                col(t) = curr_val; 
            else
                col(t) = base_val; 
            end
        end
    end
    cn0_matrix(:, s) = col; 
end

%% 4. GVI 计算 (基于绝对偏差)
fprintf('--> [特征提取] 计算绝对偏差能量 (Deviation Energy)...\n');
volatility_matrix = zeros(size(cn0_matrix));
for s = 1:num_sats
    if isnan(sat_baselines(s)), continue; end
    
    clean_col = cn0_matrix(:, s);
    base_val  = sat_baselines(s);
    
    dev_col = abs(clean_col - base_val);
    dev_col(dev_col <= PARA.noise_cutoff_db) = 0;
    
    volatility_matrix(:, s) = dev_col;
end

% GVI Sum & Segmentation
gvi_curve = sum(volatility_matrix, 2, 'omitnan');
gvi_curve_clean = movmean(gvi_curve, 5);
is_active = gvi_curve_clean > PARA.gvi_threshold;

% Merge Gap
padded = [1; is_active; 1];
gap_starts = find(diff(padded) == -1); gap_ends = find(diff(padded) == 1) - 1; 
for i = 1:length(gap_starts)
    if (gap_ends(i) - gap_starts(i) + 1) < PARA.merge_gap_pts && gap_starts(i) > 1
        is_active(gap_starts(i):gap_ends(i)-1) = 1; 
    end
end

% Extract Segments
edges = diff([0; is_active; 0]);
s_indices = find(edges == 1); e_indices = find(edges == -1) - 1;
segments = struct('id', {}, 'start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {}, 'peak_idx', {});
num_segs = 0;
for i = 1:length(s_indices)
    if (e_indices(i) - s_indices(i)) >= PARA.min_duration_pts
        num_segs = num_segs + 1;
        [max_val, max_loc] = max(gvi_curve_clean(s_indices(i):e_indices(i)));
        segments(num_segs).id = num_segs;
        segments(num_segs).start_idx = s_indices(i);
        segments(num_segs).end_idx = e_indices(i);
        segments(num_segs).peak_idx = s_indices(i) + max_loc - 1;
        segments(num_segs).peak_time = t_grid(segments(num_segs).peak_idx);
        segments(num_segs).peak_gvi = max_val;
    end
end
fprintf('✅ [Step 1] 完成: 识别到 %d 个片段。\n', num_segs);

%% 5. [核心] 数据回填 (Data Injection)
fprintf('--> [Injection] 正在将清洗后的基准线回填至 obs_clean...\n');

sat_map = containers.Map(valid_sats, 1:num_sats);

for k = 1:length(obs_clean)
    t_now = obs_clean(k).time;
    [~, t_idx] = min(abs(t_grid - t_now));
    
    epoch_sats = fieldnames(obs_clean(k).data);
    for i = 1:length(epoch_sats)
        sid = epoch_sats{i};
        if isKey(sat_map, sid)
            col_idx = sat_map(sid);
            val_clean = cn0_matrix(t_idx, col_idx);
            
            snr_struct = obs_clean(k).data.(sid).snr;
            fds = fieldnames(snr_struct);
            if ~isempty(fds)
                target_code = fds{1};
                obs_clean(k).data.(sid).snr.(target_code) = val_clean;
            end
        end
    end
end
fprintf('✅ obs_clean 数据回填完成 (包含 R 系列卫星)。\n');

%% 7. 结果可视化 (Visual Debug)
% =========================================================================
% 说明: 此模块用于绘制 GVI 曲线和阈值线，方便调试参数。
%       如果你不需要绘图，可以将下方的 if true 改为 if false。
% =========================================================================

if  false
    figure('Name', 'GVI Segmentation Analysis (Step 1)', 'Position', [100, 100, 1000, 600], 'Color', 'w');
    
    % 为了更符合直觉，将绘图时间转为北京时间 (UTC+8)
    t_grid_plot = t_grid + hours(8);
    
    % --- 子图 1: GVI 曲线与阈值 ---
    subplot(2, 1, 1);
    plot(t_grid_plot, gvi_curve_clean, 'k-', 'LineWidth', 1); hold on;
    yline(PARA.gvi_threshold, 'b--', 'Threshold', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    
    title(['1. GVI Energy Curve (Threshold = ' num2str(PARA.gvi_threshold) ')']);
    ylabel('GVI Value');
    grid on; axis tight;
    if ~isempty(t_grid_plot), datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); end
    
    % --- 子图 2: 动作片段详情 (红色高亮) ---
    subplot(2, 1, 2);
    plot(t_grid_plot, gvi_curve_clean, 'Color', [0.8 0.8 0.8], 'LineWidth', 1); hold on;
    yline(PARA.gvi_threshold, 'b--', 'LineWidth', 1);
    
    title(['2. Detected Segments (Total: ' num2str(num_segs) ')']);
    ylabel('Segment Detail');
    grid on; axis tight;
    if ~isempty(t_grid_plot), datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); end
    
    % 绘制红色高亮区域
    yl = ylim;
    for i = 1:num_segs
        idx_r = segments(i).start_idx : segments(i).end_idx;
        if isempty(idx_r), continue; end
        
        t_s = t_grid_plot(segments(i).start_idx);
        t_e = t_grid_plot(segments(i).end_idx);
        
        % 绘制半透明红色背景块
        patch([t_s t_e t_e t_s], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        % 绘制该段的高亮曲线
        plot(t_grid_plot(idx_r), gvi_curve_clean(idx_r), 'r-', 'LineWidth', 1.5);
        
        % 标记 ID
        text(t_grid_plot(segments(i).peak_idx), segments(i).peak_gvi, sprintf('#%d', segments(i).id), ...
             'Color', 'r', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
    end
    
    fprintf('✅ [Plot] GVI 可视化图表已生成。\n');
end
% =========================================================================

%% 6. 结果打包
step1_res.segments = segments;
step1_res.volatility_matrix = volatility_matrix; 
step1_res.cn0_clean_matrix  = cn0_matrix;  
step1_res.t_grid = t_grid;
step1_res.valid_sats = valid_sats;
step1_res.PARA = PARA;

end