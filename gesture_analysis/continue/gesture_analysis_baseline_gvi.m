% =========================================================================
% gesture_analysis_baseline_gvi.m (函数版)
% 功能: 手势分析 Step 1 - 基准线清洗、GVI特征提取与动作分段 (v7.2 Deviation Mode)
% 描述:
%   [重要更新] 更改了 Volatility (波动能量) 的计算方式。
%   不再使用 "当前值 - 滑动平均" (容易丢失长手势中间部分)，
%   改为 "当前值 - 全局基准线" (abs(Data - Baseline))。
%   这确保了只要信号处于遮挡状态，能量就持续存在，不会归零。
%
% [调用格式]:
%   [obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data);
%
% [输入参数]:
%   1. obs_data (struct数组): 
%      原始观测数据，由 parse_rinex_obs 读取，包含 time 和 data 字段。
%
% [返回值说明]:
%   1. obs_clean (struct): 
%      原始结构的副本 (保持结构体兼容性)。
%   2. step1_res (struct): 
%      核心结果包，包含:
%      - .cn0_clean_matrix: [关键] 经过算法清洗后的绝对 SNR 矩阵。
%      - .volatility_matrix: 波动能量矩阵 (基于绝对偏差计算)。
%      - .segments: 动作分段信息。
%      - .t_grid: 统一的时间轴 (datetime)。
%      - .valid_sats: 有效卫星列表。
%      - .PARA: 参数集 (含采样率)。
% =========================================================================

function [obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data)

fprintf('--> [Step 1] 启动基准线清洗与GVI分段 (Deviation Mode v7.2)...\n');

%% 1. 初始化与参数设置
obs_clean = obs_data; 

% --- 核心参数 ---
PARA.sampling_rate      = 25;    % [系统] 锁定 25Hz (dt=0.04s)

% [Baseline] 清洗参数
PARA.diff_lag_N         = 5;     % 趋势窗口 (N点)
PARA.noise_cutoff_db    = 1.0;   % 噪声/平稳阈值 (dB)
PARA.spike_th_db        = 1.0;   % 毛刺检测阈值 (dB)
PARA.spike_max_duration = 5;     % 毛刺最大持续点数

% [GVI & Segment] 参数
% --- 事件检测参数配置 ---
PARA.smooth_window_pts  = 25;    % [平滑窗口] 滑动平均滤波的窗口点数 (注: 在偏差模式下仅用于GVI曲线平滑以辅助检测，原始数据特征提取不受此影响)
PARA.gvi_threshold      = 2;     % [激活阈值] GVI信号判定为"动作"的触发门限值 (超过此值视为潜在动作区域)
PARA.merge_gap_pts      = 10;    % [合并容差] 两个动作片段间的最大允许间隔点数 (若间隔小于此值，则将两段合并视为同一动作，防止动作破碎)
PARA.min_duration_pts   = 3;     % [最小力度] 有效动作的最小持续点数 (持续时间短于此值的脉冲将被视为噪声过滤掉)

fprintf('    算法模式: 绝对偏差 (Abs Deviation from Baseline)\n');

%% 2. 数据对齐 (插值)
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
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

% 2.2 构建网格
raw_times = [obs_data.time];
t_grid  = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
num_samples = length(t_grid);
num_sats = length(valid_sats);

% 2.3 插值原始 SNR
cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sid = valid_sats{s_idx};
    target_code = '';
    for k = 1:min(50, length(obs_data)) 
        if isfield(obs_data(k).data, sid) && isfield(obs_data(k).data.(sid), 'snr')
            fds = fieldnames(obs_data(k).data.(sid).snr);
            if ~isempty(fds), target_code = fds{1}; break; end
        end
    end
    if isempty(target_code), continue; end
    
    s_times = []; s_vals = [];
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sid) && isfield(obs_data(k).data.(sid).snr, target_code)
            val = obs_data(k).data.(sid).snr.(target_code);
            if ~isnan(val) && val > 10
                s_times = [s_times; obs_data(k).time]; s_vals = [s_vals; val];
            end
        end
    end
    if length(s_times) > 5
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_vals(u_idx), t_grid, 'linear', NaN);
    end
end

%% 3. 基准线清洗 (验证算法逻辑)
fprintf('--> [预处理] 执行基准线清洗...\n');

N = PARA.diff_lag_N;
NoiseTh = PARA.noise_cutoff_db;
SpikeTh = PARA.spike_th_db;
SpikeDur = PARA.spike_max_duration;

% 记录每颗卫星的 Baseline 值
sat_baselines = NaN(1, num_sats);

for s = 1:num_sats
    raw_col = cn0_matrix(:, s);
    valid_mask = ~isnan(raw_col);
    if sum(valid_mask) < 10, continue; end
    
    % 计算全局基准线
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
    cn0_matrix(:, s) = col; % 更新为清洗后的数据
end

%% 4. [核心修改] GVI 计算 (基于绝对偏差)
% 算法: volatility = abs(current - global_baseline)
fprintf('--> [特征提取] 计算绝对偏差能量 (Deviation Energy)...\n');

volatility_matrix = zeros(size(cn0_matrix));

for s = 1:num_sats
    if isnan(sat_baselines(s)), continue; end
    
    % 取出清洗后的数据
    clean_col = cn0_matrix(:, s);
    base_val  = sat_baselines(s);
    
    % 计算偏差能量
    dev_col = abs(clean_col - base_val);
    
    % 双重保险: 再次过滤微小噪声
    dev_col(dev_col <= PARA.noise_cutoff_db) = 0;
    
    volatility_matrix(:, s) = dev_col;
end

% 计算 GVI 总和
gvi_curve = sum(volatility_matrix, 2, 'omitnan');
% 平滑 GVI 曲线
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

% Segments
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

%% 5. 结果打包
step1_res.segments = segments;
step1_res.volatility_matrix = volatility_matrix; 
step1_res.cn0_clean_matrix  = cn0_matrix;  
step1_res.t_grid = t_grid;
step1_res.valid_sats = valid_sats;
step1_res.PARA = PARA;

end