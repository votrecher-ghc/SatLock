% =========================================================================
% step1_segmentation_GVI.m (函数版)
% 功能: 手势分析第一步 - 信号预处理与分段 (v5.1 Fixed)
%
% [调用格式]:
%   [segments, volatility_matrix, t_grid, valid_sats, PARA] = step1_segmentation_GVI(obs_data);
%
% [输入参数]:
%   1. obs_data (struct): 原始观测数据。
%
% [返回值说明]:
%   1. segments (struct): 识别出的手势片段信息（起止索引、峰值时间等）。
%   2. volatility_matrix (double): 去噪后的信号波动矩阵 [Samples x Sats]。
%   3. t_grid (datetime): 统一的时间轴网格。
%   4. valid_sats (cell): 有效卫星ID列表，与矩阵列对应。
%   5. PARA (struct): 计算过程中确定的参数（如自动计算的 sampling_rate）。
% =========================================================================

function [segments, volatility_matrix, t_grid, valid_sats, PARA] = step1_segmentation_GVI(obs_data)

fprintf('--> [Step 1] 启动手势检测与分段 (Function版 v5.1)...\n');

%% 1. 参数设置
PARA.smooth_window_sec = 1.5;  % 基线平滑窗口(秒)
PARA.gvi_threshold     = 6;    % GVI 阈值 (dB)
PARA.sampling_rate     = 25;   % 默认采样率 (会被自动覆盖)
PARA.merge_gap_sec     = 0.5;  % 连通域合并窗口
PARA.min_duration_sec  = 0.4;  % 最小持续时间

%% 2. 数据提取与对齐
all_sat_ids = {};
% 扫描数据以获取所有卫星ID (取首尾部分扫描以加速)
scan_range = unique([1:min(100, length(obs_data)), max(1, length(obs_data)-100):length(obs_data)]);
for i = scan_range
    if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end
end
unique_sat_ids = unique(all_sat_ids);

valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

% === 采样率自动校准 ===
raw_times = [obs_data.time];
mean_dt = mean(diff(raw_times));
fprintf('    检测到平均时间间隔: %.4f 秒\n', seconds(mean_dt));

if mean_dt < seconds(0.03), PARA.sampling_rate = 50;
elseif mean_dt < seconds(0.05), PARA.sampling_rate = 25;
elseif mean_dt < seconds(0.12), PARA.sampling_rate = 10;
else, PARA.sampling_rate = 1; end
fprintf('    -> 自动匹配采样率为: %d Hz\n', PARA.sampling_rate);

% 时间网格构建
t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
t_grid_plot = t_grid + hours(8) - seconds(20); % 北京时间用于绘图

num_samples = length(t_grid);
num_sats = length(valid_sats);

% 提取 C/N0 矩阵
cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    for k = 1:min(50, length(obs_data)) 
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            fields = fieldnames(obs_data(k).data.(sat_id).snr);
            if ~isempty(fields)
                if ismember('S1C', fields), target_snr_code = 'S1C'; 
                elseif ismember('S2I', fields), target_snr_code = 'S2I';
                else, target_snr_code = fields{1}; end
                break;
            end
        end
    end
    if isempty(target_snr_code), continue; end
    
    s_times = []; s_cn0 = [];
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id).snr, target_snr_code)
            val = obs_data(k).data.(sat_id).snr.(target_snr_code);
            if ~isnan(val) && val > 10
                s_times = [s_times; obs_data(k).time]; s_cn0 = [s_cn0; val];
            end
        end
    end
    if length(s_times) > 20
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
    end
end

%% 3. 核心预处理 (SG滤波)
fprintf('--> [预处理] 执行 Savitzky-Golay 滤波...\n');
sg_order = 2; sg_len = 7; 

for s = 1:num_sats
    col_data = cn0_matrix(:, s); valid_mask = ~isnan(col_data);
    if sum(valid_mask) > sg_len * 2
        x = 1:length(col_data);
        filled_data = interp1(x(valid_mask), col_data(valid_mask), x, 'linear', 'extrap')';
        cn0_matrix(:, s) = sgolayfilt(filled_data, sg_order, sg_len);
    end
end

%% 4. 计算波动与分段
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth);
gvi_curve = sum(volatility_matrix, 2, 'omitnan');
gvi_curve_clean = movmean(gvi_curve, 5); 

fprintf('--> 正在执行连通域分段...\n');
is_active = gvi_curve_clean > PARA.gvi_threshold;

% Merge Gap
min_gap_idx = round(PARA.merge_gap_sec * PARA.sampling_rate);
padded_active = [1; is_active; 1];
gap_starts = find(diff(padded_active) == -1); gap_ends = find(diff(padded_active) == 1) - 1; 
for i = 1:length(gap_starts)
    if (gap_ends(i) - gap_starts(i) + 1) < min_gap_idx && gap_starts(i) > 1, is_active(gap_starts(i):gap_ends(i)-1) = 1; end
end

% Extract Segments
edges = diff([0; is_active; 0]);
s_indices = find(edges == 1); e_indices = find(edges == -1) - 1;
segments = struct('id', {}, 'start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {}, 'peak_idx', {});
num_segs = 0;
min_dur_idx = round(PARA.min_duration_sec * PARA.sampling_rate);

for i = 1:length(s_indices)
    if (e_indices(i) - s_indices(i)) >= min_dur_idx
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
fprintf('✅ 识别到 %d 个手势片段。\n', num_segs);

%% 5. 结果可视化 (保留绘图功能)
figure('Name', 'Gesture Detection Analysis (v5.1 Fixed)', 'Position', [50, 50, 1000, 800], 'Color', 'w');
subplot(3, 1, 1); plot(t_grid_plot, cn0_matrix); 
title('1. 滤波后的全星座 C/N0 数据'); ylabel('SNR'); datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); axis tight; grid on;

subplot(3, 1, 2); plot(t_grid_plot, gvi_curve_clean, 'k-'); hold on; yline(PARA.gvi_threshold, 'b--');
title('2. 波动指数 (GVI)'); ylabel('GVI'); datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); axis tight; grid on;

subplot(3, 1, 3); plot(t_grid_plot, gvi_curve_clean, 'Color', [0.8 0.8 0.8]); hold on;
yline(PARA.gvi_threshold, 'b--'); 
title('3. 动作片段详情'); ylabel('GVI Detail'); datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); grid on; axis tight;
yl3 = ylim; 
for i = 1:num_segs
    idx_r = segments(i).start_idx : segments(i).end_idx;
    t_s = t_grid_plot(segments(i).start_idx); t_e = t_grid_plot(segments(i).end_idx);
    patch([t_s t_e t_e t_s], [yl3(1) yl3(1) yl3(2) yl3(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    plot(t_grid_plot(idx_r), gvi_curve_clean(idx_r), 'r-', 'LineWidth', 2);
    text(t_grid_plot(segments(i).peak_idx), segments(i).peak_gvi, sprintf('#%d', i), 'Color', 'r', 'VerticalAlignment', 'bottom');
end

fprintf('✅ Step 1 分析完成 (已返回分段数据)。\n');
end