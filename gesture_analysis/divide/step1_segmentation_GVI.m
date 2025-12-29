% =========================================================================
% step1_segmentation_GVI.m (v5.1 - 修复采样率Bug版)
% 功能: 手势检测 Step 1 
% 改进记录:
%   v5.1: [修复] 修正采样率自动识别逻辑，正确支持 25Hz 数据 (dt=0.04s)。
%   v5.0: [新增] 引入 Savitzky-Golay 预滤波，消除 SNR 量化噪声。
%   v4.x: [保持] 连通域合并逻辑，防止动作断裂。
%   v4.x: [保持] 统一的北京时间轴 (UTC+8 -20s)。
% =========================================================================

%% 1. 准备工作与参数设置
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end

% --- 核心参数 ---
PARA.smooth_window_sec = 1.5;  % 基线平滑窗口(秒): 用于计算长时间背景基线
PARA.gvi_threshold     = 6;    % GVI 阈值 (dB): 超过此值的波动才会被检测

% [修正] 默认采样率设为 25，后续会根据实际数据覆盖
PARA.sampling_rate     = 25;   

% 连通域参数 (解决 M 型波断裂)
PARA.merge_gap_sec     = 0.5;  % 中间断开小于 0.5 秒视为不断开
PARA.min_duration_sec  = 0.4;  % 小于 0.4 秒视为噪声

fprintf('--> [Step 1] 开始手势检测 (v5.1 SG滤波 + 25Hz修复)...\n');

%% 2. 数据提取与对齐
all_sat_ids = {};
% 扫描数据以获取所有卫星ID
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

% === [核心修复] 采样率自动校准 ===
raw_times = [obs_data.time];
mean_dt = mean(diff(raw_times));
fprintf('    检测到数据平均时间间隔: %.4f 秒\n', seconds(mean_dt));

if mean_dt < seconds(0.03)      % < 30ms (如 50Hz, dt=0.02)
    PARA.sampling_rate = 50;
elseif mean_dt < seconds(0.05)  % < 50ms (如 25Hz, dt=0.04) -> 命中这里
    PARA.sampling_rate = 25;
elseif mean_dt < seconds(0.12)  % < 120ms (如 10Hz, dt=0.10)
    PARA.sampling_rate = 10;
else
    PARA.sampling_rate = 1; 
end
fprintf('    -> 自动匹配采样率为: %d Hz\n', PARA.sampling_rate);

% 原始时间网格 (UTC)
t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 

% 构建统一的绘图时间轴 (北京时间 + 20s跳秒修正)
t_grid_plot = t_grid + hours(8) - seconds(20); 

num_samples = length(t_grid);
num_sats = length(valid_sats);

% 提取 C/N0 矩阵
cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    % 智能搜索可用的 SNR 信号类型
    for k = 1:min(50, length(obs_data)) 
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            fields = fieldnames(obs_data(k).data.(sat_id).snr);
            if ~isempty(fields)
                if ismember('S1C', fields), target_snr_code = 'S1C'; 
                elseif ismember('S2I', fields), target_snr_code = 'S2I';
                else, target_snr_code = fields{1}; 
                end
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
                s_times = [s_times; obs_data(k).time]; %#ok<AGROW>
                s_cn0   = [s_cn0; val]; %#ok<AGROW>
            end
        end
    end
    if length(s_times) > 20
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
    end
end

%% =================== 2.5 核心预处理 (SG滤波) ===================
fprintf('--> [预处理] 正在对 C/N0 矩阵进行 Savitzky-Golay 滤波...\n');

% 备份原始数据 (用于绘图对比)
cn0_matrix_raw = cn0_matrix; 

% SG 滤波参数: 2阶多项式，7点窗口 (约0.28秒 @ 25Hz)
% 作用: 填平量化台阶，去除单点毛刺，保留手势波形特征
sg_order = 2;
sg_len = 7; 

for s = 1:num_sats
    col_data = cn0_matrix(:, s);
    valid_mask = ~isnan(col_data);
    
    if sum(valid_mask) > sg_len * 2
        % 临时填补 NaN 以进行滤波 (线性插值)
        x = 1:length(col_data);
        filled_data = interp1(x(valid_mask), col_data(valid_mask), x, 'linear', 'extrap')';
        
        % 滤波
        filtered_col = sgolayfilt(filled_data, sg_order, sg_len);
        
        % 将滤波结果写回
        cn0_matrix(:, s) = filtered_col;
    end
end
fprintf('    滤波完成。现在的 cn0_matrix 更加平滑，底噪更低。\n');

%% 3. 计算波动与分段 (连通域逻辑)
% 这里的 cn0_matrix 已经是去噪后的版本
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth);
gvi_curve = sum(volatility_matrix, 2, 'omitnan');
gvi_curve_clean = movmean(gvi_curve, 5); % 再做一次轻微平滑

fprintf('--> 正在执行连通域分段...\n');

is_active = gvi_curve_clean > PARA.gvi_threshold;

% 填补缝隙 (Merge Gap)
min_gap_idx = round(PARA.merge_gap_sec * PARA.sampling_rate);
padded_active = [1; is_active; 1];
gap_starts = find(diff(padded_active) == -1); 
gap_ends   = find(diff(padded_active) == 1) - 1; 

for i = 1:length(gap_starts)
    len = gap_ends(i) - gap_starts(i) + 1;
    if len < min_gap_idx && gap_starts(i) > 1 && gap_ends(i) < length(padded_active)
        idx_s = gap_starts(i); 
        idx_e = gap_ends(i) - 1; 
        if idx_s <= idx_e, is_active(idx_s:idx_e) = 1; end
    end
end

% 提取分段 (Segments)
edges = diff([0; is_active; 0]);
s_indices = find(edges == 1);
e_indices = find(edges == -1) - 1;

segments = struct('id', {}, 'start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {}, 'peak_idx', {});
num_segs = 0;
min_dur_idx = round(PARA.min_duration_sec * PARA.sampling_rate);

for i = 1:length(s_indices)
    s = s_indices(i);
    e = e_indices(i);
    if (e - s) >= min_dur_idx
        num_segs = num_segs + 1;
        [max_val, max_loc] = max(gvi_curve_clean(s:e));
        peak_idx = s + max_loc - 1;
        segments(num_segs).id = num_segs;
        segments(num_segs).start_idx = s;
        segments(num_segs).end_idx = e;
        segments(num_segs).peak_idx = peak_idx;
        segments(num_segs).peak_time = t_grid(peak_idx); % 存储 UTC 时间供索引
        segments(num_segs).peak_gvi = max_val;
    end
end

fprintf('✅ 识别到 %d 个手势片段。\n', num_segs);

%% 4. 结果可视化
figure('Name', 'Gesture Detection Analysis (v5.1 Fixed)', 'Position', [50, 50, 1000, 800], 'Color', 'w');

% --- 图1: C/N0 数据对比 ---
subplot(3, 1, 1);
plot(t_grid_plot, cn0_matrix); % 画滤波后的数据
title(sprintf('1. 滤波后的全星座 C/N0 数据 (SG Filtered, %d sats)', num_sats));
ylabel('SNR (Filtered)'); 
xlabel('时间 (北京时间)');
datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
axis tight; grid on;

yl = ylim;
for i=1:num_segs
    t_s = t_grid_plot(segments(i).start_idx);
    t_e = t_grid_plot(segments(i).end_idx);
    patch([t_s t_e t_e t_s], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

% --- 图2: 波动指数 GVI ---
subplot(3, 1, 2);
plot(t_grid_plot, gvi_curve_clean, 'k-', 'LineWidth', 1); hold on;
yline(PARA.gvi_threshold, 'b--', '阈值'); % 默认蓝色
title(sprintf('2. 波动指数 (基于滤波数据), 阈值:%.1f', PARA.gvi_threshold));
ylabel('GVI'); 
xlabel('时间 (北京时间)');
datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
axis tight; grid on;

yl2 = ylim;
for i=1:num_segs
    t_s = t_grid_plot(segments(i).start_idx);
    t_e = t_grid_plot(segments(i).end_idx);
    patch([t_s t_e t_e t_s], [yl2(1) yl2(1) yl2(2) yl2(2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% --- 图3：手势片段详情 ---
subplot(3, 1, 3);
plot(t_grid_plot, gvi_curve_clean, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); hold on;
yline(PARA.gvi_threshold, 'b--', '阈值'); 

title('3. 检测到的手势动作片段详情 (红色高亮)');
ylabel('GVI Detail'); 
xlabel('时间 (北京时间)');
datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
grid on; axis tight;

yl3 = ylim; 
for i = 1:num_segs
    idx_range = segments(i).start_idx : segments(i).end_idx;
    t_s = t_grid_plot(segments(i).start_idx);
    t_e = t_grid_plot(segments(i).end_idx);
    t_peak = t_grid_plot(segments(i).peak_idx);
    
    % 1. 红色背景块
    patch([t_s t_e t_e t_s], [yl3(1) yl3(1) yl3(2) yl3(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
      
    % 2. 红色波形线
    plot(t_grid_plot(idx_range), gvi_curve_clean(idx_range), 'r-', 'LineWidth', 2);
    
    % 3. 标签
    text(t_peak, segments(i).peak_gvi, sprintf('  #%d', i), ...
        'Color', 'r', 'FontWeight', 'bold', 'FontSize', 11, ...
        'VerticalAlignment', 'bottom');
end

%% 5. 计算并导出起止时间
fprintf('\n--> 正在生成结果表格...\n');

ids = [segments.id]';
start_idxs = [segments.start_idx]';
end_idxs   = [segments.end_idx]';
dur_samples = end_idxs - start_idxs + 1;

start_times_bjt = t_grid_plot(start_idxs);
end_times_bjt   = t_grid_plot(end_idxs);
dur_seconds = seconds(end_times_bjt - start_times_bjt);

T_Index = table(ids, start_idxs, end_idxs, dur_samples, ...
    'VariableNames', {'GestureID', 'Start_Index', 'End_Index', 'Duration_Points'});

T_Time = table(ids, start_times_bjt, end_times_bjt, dur_seconds, ...
    'VariableNames', {'GestureID', 'Start_Time_BJT', 'End_Time_BJT', 'Duration_Sec'});

fprintf('\n=== 表格 1: 采样点信息 ===\n');
disp(T_Index);

fprintf('\n=== 表格 2: 北京时间信息 (UTC+8 -20s修正) ===\n');
disp(T_Time);