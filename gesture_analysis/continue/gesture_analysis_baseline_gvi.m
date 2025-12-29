% =========================================================================
% gesture_analysis_baseline_gvi.m (函数版)
% 功能: 手势分析 Step 1 - 基准线清洗、GVI计算与分段 (v5.4 Points Param)
% 描述:
%   该函数作为手势分析流水线的第一步，负责从原始观测数据中提取干净的信号特征。
%   它集成了基准线锁定、智能底噪清洗、SG滤波去噪以及基于GVI的动作区间检测。
%   [更新] 参数配置已改为直接使用"采样点数"，不再使用"秒"。
%
% [调用格式]:
%   [obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data);
%
% [输入参数]:
%   1. obs_data (struct数组): 原始观测数据。
%
% [返回值说明]:
%   1. obs_clean (struct): 清洗后的观测数据 (SNR 已平滑)。
%   2. step1_res (struct): 打包好的中间结果 (含 segments, volatility, t_grid 等)。
% =========================================================================

function [obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data)

fprintf('--> [Step 1] 启动基准线清洗与GVI分段 (Points Mode)...\n');

%% 1. 初始化与参数设置
obs_clean = obs_data; % 复制一份用于回写

% --- 核心参数 (点数制) ---
% 基础设定: 采样率固定 25Hz (即 1点 = 0.04秒)

% [Baseline] 基准线算法参数
PARA.diff_lag_N         = 5;     % 趋势窗口 (点数)
PARA.noise_cutoff_db    = 1;     % 噪声阈值 (dB)
PARA.spike_th_db        = 2;     % 毛刺阈值 (dB)
PARA.spike_max_duration = 3;     % 毛刺最大宽度 (点数)

% [GVI & Segment] 分段算法参数 (已修改为点数)
PARA.smooth_window_pts  = 25;    % [GVI] 平滑窗口 (点数): 38点 ≈ 1.5秒
PARA.gvi_threshold      = 2;     % [GVI] 激活阈值 (dB)
PARA.merge_gap_pts      = 10;    % [Segment] 合并间隙 (点数): 10点 = 0.4秒
PARA.min_duration_pts   = 5;     % [Segment] 最小动作时长 (点数): 5点 = 0.2秒

% [System] 采样率锁定
PARA.sampling_rate      = 25;    
fprintf('    系统采样率已锁定: %d Hz (dt=0.04s)\n', PARA.sampling_rate);

%% 2. 数据对齐
% 提取有效卫星
all_sat_ids = {};
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

% 构建标准时间网格 (严格按照 25Hz)
raw_times = [obs_data.time];
t_start = min(raw_times);
t_end   = max(raw_times);
t_grid  = (t_start : seconds(1/PARA.sampling_rate) : t_end)'; 

num_samples = length(t_grid);
num_sats = length(valid_sats);

% 提取并插值 SNR 矩阵
cn0_matrix = NaN(num_samples, num_sats);

for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    
    % 智能查找 SNR 字段
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
    
    % 提取原始数据
    s_times = []; s_cn0 = [];
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id).snr, target_snr_code)
            val = obs_data(k).data.(sat_id).snr.(target_snr_code);
            if ~isnan(val) && val > 10
                s_times = [s_times; obs_data(k).time]; s_cn0 = [s_cn0; val];
            end
        end
    end
    
    % 线性插值到标准 25Hz 网格
    if length(s_times) > 20
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
    end
end

%% 3. 核心算法: 基准线清洗 (Baseline Cleaning)
fprintf('--> [预处理] 执行基准线清洗 (Mode Locking)...\n');

for s = 1:num_sats
    raw_col = cn0_matrix(:, s); col = raw_col;
    valid_data = raw_col(~isnan(raw_col)); 
    if isempty(valid_data), continue; end
    
    % 计算众数基准线
    baseline = mode(round(valid_data)); 
    
    for t = 1:num_samples
        curr_val = raw_col(t); if isnan(curr_val), continue; end
        
        % 毛刺剔除逻辑
        if abs(curr_val - baseline) > PARA.spike_th_db
             is_spike = true;
             for k=1:PARA.spike_max_duration
                 if t+k<=num_samples && abs(raw_col(t+k)-baseline)<=PARA.noise_cutoff_db
                     is_spike=true; break; 
                 end
             end
             % 暂时保留显著跳变
        end
        
        % 底噪清洗: 强制归零微小震荡
        if abs(curr_val - baseline) < PARA.noise_cutoff_db
            col(t) = baseline; 
        end
    end
    cn0_matrix(:, s) = col;
end

%% 4. 计算能量与分段
% 计算波动能量 (Volatility)
% [修正] 直接使用点数参数 PARA.smooth_window_pts
cn0_smooth = movmean(cn0_matrix, PARA.smooth_window_pts, 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth); 

% 计算 GVI
gvi_curve = sum(volatility_matrix, 2, 'omitnan');
gvi_curve_clean = movmean(gvi_curve, 5); 

is_active = gvi_curve_clean > PARA.gvi_threshold;

% Merge Gap (缝合间隙)
% [修正] 直接使用点数参数 PARA.merge_gap_pts
min_gap_idx = PARA.merge_gap_pts;
padded_active = [1; is_active; 1];
gap_starts = find(diff(padded_active) == -1); 
gap_ends   = find(diff(padded_active) == 1) - 1; 

for i = 1:length(gap_starts)
    gap_len = gap_ends(i) - gap_starts(i) + 1;
    if gap_len < min_gap_idx && gap_starts(i) > 1
        is_active(gap_starts(i):gap_ends(i)-1) = 1; 
    end
end

% Extract Segments (提取片段)
edges = diff([0; is_active; 0]);
s_indices = find(edges == 1); 
e_indices = find(edges == -1) - 1;

segments = struct('id', {}, 'start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {}, 'peak_idx', {});
num_segs = 0;
% [修正] 直接使用点数参数 PARA.min_duration_pts
min_dur_idx = PARA.min_duration_pts;

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
fprintf('✅ [Step 1] 完成: 识别到 %d 个手势片段。\n', num_segs);

%% 5. 结果打包
step1_res.segments = segments;
step1_res.volatility_matrix = volatility_matrix; 
step1_res.t_grid = t_grid;
step1_res.valid_sats = valid_sats;
step1_res.PARA = PARA; % 包含修改后的点数参数

end